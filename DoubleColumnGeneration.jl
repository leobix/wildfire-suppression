include("CommonStructs.jl")
include("Subproblems.jl")
include("BranchingRules.jl")

using JuMP


@kwdef mutable struct DualWarmStart
	# Linking duals from a previous solve that we can recycle as a warm start.
	# Re-using dual information is especially helpful when exploring nearby nodes
	# of the branch-and-price tree where the optimal duals change slowly.
	linking_values::Matrix{Float64}
	# Strategy flag that can be used by callers to signal how the warm start
	# was generated (currently just informative, but kept for future tuning).
	const strategy::String = "global"
	# Small perturbation that can be applied to avoid ties in the duals.
	epsilon::Float64 = 0.001
end

function adapt_linking_duals(warm_start_matrix::Matrix{Float64}, num_fires::Int, num_time_periods::Int)
	# Resize the incoming warm-start matrix if the current instance has a different
	# number of fires or time periods.  Missing entries default to zero.
	result = zeros(num_fires, num_time_periods)
	min_f = min(size(warm_start_matrix, 1), num_fires)
	min_t = min(size(warm_start_matrix, 2), num_time_periods)
	result[1:min_f, 1:min_t] .= warm_start_matrix[1:min_f, 1:min_t]
	return result
end

"""
Perform double column generation to solve the (relaxed) restricted master problem `rmp`
by iteratively adding columns to model feasible fire suppression plans and crew routes.
This function optimizes suppression plans and routes simultaneously to satisfy
supply-demand constraints and minimize costs.

# Arguments

  - `rmp::RestrictedMasterProblem` (modified): The restricted master problem instance representing the optimization problem.
  - `crew_routes::CrewRouteData` (modified): Data structure containing information (cost, assignments) about all generated crew routes.
  - `fire_plans::FirePlanData` (modified): Data structure containing information (cost, demands) about all generated fire plans.
  - `cut_data::CutData` (modified): Data structure containing information about cuts, including lookups to incorporate coeffs into master problem and subproblems.
  - `crew_subproblems::Vector{TimeSpaceNetwork}`: Vector of data structures containing all static info about crew subproblems.
  - `fire_subproblems::Vector{TimeSpaceNetwork}`: Vector of data structures containing all static info about fire subproblems.
  - `crew_branching_rules::Vector{CrewAssignmentBranchingRule}`: Vector of branching rules that indicate whether crew j suppresses fire g at time t.
  - `fire_branching_rules::Vector{FireDemandBranchingRule}`: Vector of branching rules that indicate whether fire g demands <=d crews or >d crews at time t.
  - `global_fire_allotment_branching_rules::Vector{GlobalFireAllotmentBranchingRule}`: Vector of branching rules that may place a cap on demand across all fires and times
  - `upper_bound::Float64`: Global upper bound given by the best feasible solution we have to the integer problem.
  - `timing::Bool`: Flag for whether we will track and return detailed timings
  - `improving_column_abs_tolerance::Float64`: Absolute improvement needed (using reduced cost bound) to add a column (default: 1e-10).
  - `local_gap_rel_tolerance::Float64`: Relative tolerance (using Lagrangian bound) needed to accept solution as optimal (default: 1e-9).
"""
function double_column_generation!!!!(
	rmp::RestrictedMasterProblem,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	cut_data::CutData,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	crew_branching_rules::Vector{CrewAssignmentBranchingRule},
	fire_branching_rules::Vector{FireDemandBranchingRule},
	global_fire_allotment_branching_rules::Vector{GlobalFireAllotmentBranchingRule},
	fires_to_ignore::Vector{Int64};
	upper_bound::Float64,
	timing::Bool,
	time_limit::Float64 = Inf,
	improving_column_abs_tolerance::Float64 = 1e-10,
	local_gap_rel_tolerance::Float64 = 1e-5,
	final_snapshot_only::Bool = false,
	dual_warm_start::Union{Nothing, DualWarmStart} = nothing)

	# initialize timing dictionary so that callers can understand where the time goes
	details = Dict{String, Float64}()
	details["master_problem"] = 0.0
	details["fire_subproblems"] = 0.0
	details["crew_subproblems"] = 0.0
	t = time()

	# gather global information about the dimensionality of the problem
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	# If the master problem has not been solved yet, create an artificial dual
	# point.  The values do not need to be feasible; they simply guide the
	# pricing problems toward useful early columns.
	if rmp.termination_status == MOI.OPTIMIZE_NOT_CALLED
		if isnothing(dual_warm_start)
			# initialize with an (infeasible) dual solution that will suppress minimally
			fire_duals = zeros(num_fires) .+ Inf
			crew_duals = zeros(num_crews)
			linking_duals = zeros(num_fires, num_time_periods) .+ 1e30
			cut_duals = normalized_rhs.(rmp.gub_cover_cuts) .* 0
			global_fire_allot_duals = normalized_rhs.(rmp.fire_allotment_branches) .* 0
		else
			fire_duals = zeros(num_fires)
			crew_duals = zeros(num_crews)
			linking_duals = adapt_linking_duals(dual_warm_start.linking_values, num_fires, num_time_periods)
			cut_duals = normalized_rhs.(rmp.gub_cover_cuts) .* 0
			global_fire_allot_duals = normalized_rhs.(rmp.fire_allotment_branches) .* 0
		end
	else
		# Otherwise, recycle the duals from the last master problem solve. This
		# is the standard approach in column generation because the duals define
		# the reduced-cost pricing problems.
		fire_duals = dual.(rmp.plan_per_fire)
		crew_duals = dual.(rmp.route_per_crew)
		linking_duals = dual.(rmp.supply_demand_linking)
		cut_duals = dual.(rmp.gub_cover_cuts)
		global_fire_allot_duals = dual.(rmp.fire_allotment_branches)
	end

	# Track the incumbent upper bound separately so we can easily scale duals
	# when infeasibilities are encountered.
	ub = copy(upper_bound)

	# initialize column generation loop
	continue_iterating::Bool = true
	iteration = 0

	# add in dummy plans for the fires_to_ignore
	# We enforce that these fires always have at least one trivial plan so that
	# the linking constraints remain well-defined even though we do not want to
	# explicitly price them in this call.
	for fire in fires_to_ignore
		add_column_to_master_problem!!(
			rmp,
			cut_data,
			fire_plans,
			global_fire_allotment_branching_rules,
			fire,
			1,
		)
	end

	# Main column-generation loop: price crew routes, then fire plans, then
	# re-optimize the restricted master problem until no improving columns exist.
	while continue_iterating

		iteration += 1
		reduced_cost_sum = 0

		if timing
			t = time()
		end

		crew_objectives = zeros(Float64, num_crews)
		crew_arcs_used = [Int[] for crew ∈ 1:num_crews]

		# Solve every crew's resource-allocation pricing problem in parallel.
		# Each crew contributes a single route/column back to the master problem.
		#
		# NOTE: JuMP/Gurobi handles the master problem, so there is no need for
		# locks here—the subproblems are fully independent.
		Threads.@threads for crew in 1:num_crews

			# Reset all per-iteration data on the network so we start from a
			# clean slate before dual adjustments are applied.
			crew_subproblems[crew].prohibited_arcs .&= false
			for i ∈ eachindex(crew_subproblems[crew].modified_arc_costs)
				crew_subproblems[crew].modified_arc_costs[i] = crew_subproblems[crew].arc_costs[i]
			end
			# Incorporate supply/demand dual prices and branching decisions into
			# the arc costs so the shortest-path solution returns the reduced cost.
			adjust_crew_arc_costs!!(
				crew_subproblems[crew].modified_arc_costs,
				crew_subproblems[crew].prohibited_arcs,
				crew,
				linking_duals,
				crew_subproblems[crew].supply_demand_dual_arc_lookup,
				crew_branching_rules,
			)

			# adjust the arc costs for the cuts
			cut_adjust_arc_costs!(
				crew_subproblems[crew].modified_arc_costs,
				cut_data.crew_sp_lookup[crew],
				cut_duals,
			)

			# solve the subproblem
			objective, arcs_used = crew_dp_subproblem(
				crew_subproblems[crew].wide_arcs,
				crew_subproblems[crew].modified_arc_costs,
				crew_subproblems[crew].prohibited_arcs,
				crew_subproblems[crew].state_in_arcs,
			)

			# Adjust the objective for the cuts (we gave - coeff if allotment
			# not broken, give + coeff here) so that the reduced cost matches
			# the dual constraints that involve this route.
			for (ix, cut) in cut_data.cut_dict
				if cut.crew_coeffs[crew] > 1e-20
					objective += (cut.crew_coeffs[crew] * cut_duals[ix])
				end
			end

			crew_objectives[crew] = objective
			crew_arcs_used[crew] = arcs_used
		end
		

		for crew ∈ 1:num_crews

			objective = crew_objectives[crew]
			arcs_used = crew_arcs_used[crew]

			# if there is an improving route
			# (Reduced cost < 0 for minimization, so objective < crew dual.)
			if objective < crew_duals[crew] - improving_column_abs_tolerance

				reduced_cost_sum += (objective - crew_duals[crew])

				# get the real cost, unadjusted for duals
				# (The state machine only considers reduced cost; the master
				# problem objective needs the original physical cost.)
				cost = sum(crew_subproblems[crew].arc_costs[arcs_used])

				# get the indicator matrix of fires fought at each time
				fires_fought = get_fires_fought(
					crew_subproblems[crew].wide_arcs,
					arcs_used,
					(num_fires, num_time_periods),
				)

				@debug "crew route" crew fires_fought

				# add the route to the routes
				new_route_ix =
					add_column_to_route_data!(crew_routes, crew, cost, fires_fought, arcs_used)

				# update the master problem
				add_column_to_master_problem!!(
					rmp,
					cut_data,
					crew_routes,
					crew,
					new_route_ix,
				)
			end
		end

		if timing
			details["crew_subproblems"] += time() - t
			t = time()
		end

		# Solve fire pricing problems next.  They share the same structure as
		# the crew problems but demand different data and dual adjustments.
		fire_objectives = zeros(Float64, num_fires)
		fire_arcs_used = [Int[] for fire ∈ 1:num_fires]

		# for each fire
		# Fires that should be ignored are skipped entirely since they already
		# have dummy plans in the master problem.
		Threads.@threads for fire in [g for g ∈ 1:num_fires if g ∉ fires_to_ignore]

			# generate the local costs of the arcs
			fire_subproblems[fire].prohibited_arcs .&= false
			for i ∈ eachindex(fire_subproblems[fire].modified_arc_costs)
				fire_subproblems[fire].modified_arc_costs[i] = fire_subproblems[fire].arc_costs[i]
			end
			adjust_fire_arc_costs!!(
				fire_subproblems[fire].modified_arc_costs,
				fire_subproblems[fire].prohibited_arcs,
				fire,
				fire_subproblems[fire].supply_demand_dual_arc_lookup,
				fire_subproblems[fire].long_arcs,
				linking_duals[fire, :],
				[rule for rule ∈ fire_branching_rules if rule.fire_ix == fire],
			)

			# adjust the arc costs for the cuts
			cut_adjust_arc_costs!(
				fire_subproblems[fire].modified_arc_costs,
				cut_data.fire_sp_lookup[fire],
				cut_duals,
			)

			# for each branching rule
			for ix in eachindex(global_fire_allotment_branching_rules)

				rule = global_fire_allotment_branching_rules[ix]

				# if this is a <= 0 rule, we have more prohibited arcs
				if ~rule.geq_flag
					for arc_ix ∈ rule.fire_sp_arc_lookup[fire]
						fire_subproblems[fire].prohibited_arcs[arc_ix] = true
					end
				end

				# we need to do a proper dual adjustment
				# The branching rules restrict cumulative demand across several
				# fires; we use an extra set of dual multipliers to communicate
				# them down to each fire subproblem.
				adjust_fire_sp_arc_costs!(
					fire_subproblems[fire].modified_arc_costs,
					rule,
					fire,
					global_fire_allot_duals[ix],
				)
			end

			# solve the subproblem
			objective, arcs_used = fire_dp_subproblem(
				fire_subproblems[fire].wide_arcs,
				fire_subproblems[fire].modified_arc_costs,
				fire_subproblems[fire].prohibited_arcs,
				fire_subproblems[fire].state_in_arcs,
			)

			fire_objectives[fire] = objective
			fire_arcs_used[fire] = arcs_used
		end
		


		for fire ∈ [fire for fire ∈ 1:num_fires if fire ∉ fires_to_ignore]

			objective = fire_objectives[fire]
			arcs_used = fire_arcs_used[fire]

			# if there is an improving plan
			if objective < fire_duals[fire] - improving_column_abs_tolerance

				reduced_cost_sum += (objective - fire_duals[fire])

				# get the real cost, unadjusted for duals
				cost = 0.0
				if final_snapshot_only
					# choose last nonzero (pre-extinguish) end-of-day area as cost
					# This mode is used when only the final fire size matters.
					best_arc = 0
					best_t = -1
					for a in arcs_used
						to_t = fire_subproblems[fire].long_arcs[a, FM.TIME_TO]
						c = fire_subproblems[fire].arc_costs[a]
						if to_t <= num_time_periods + 1 && (to_t - 1) <= num_time_periods && c > 1e-12
							if (to_t - 1) > best_t
								best_t = to_t - 1
								best_arc = a
							end
						end
					end
					if best_arc != 0
						cost = fire_subproblems[fire].arc_costs[best_arc]
					else
						cost = sum(fire_subproblems[fire].arc_costs[arcs_used])
					end
				else
					cost = sum(fire_subproblems[fire].arc_costs[arcs_used])
				end

				# get the vector of crew demands at each time
				# Each plan is summarized by how many crews it would like at
				# every period; the master problem uses these as column data.
				crew_demands = get_crew_demands(
					fire_subproblems[fire].wide_arcs,
					arcs_used,
					num_time_periods,
				)

				# add the plan to the plans
				new_plan_ix =
					add_column_to_plan_data!(fire_plans, fire, cost, crew_demands, arcs_used)

				@debug "fire plan" fire crew_demands
				# update the master problem
				add_column_to_master_problem!!(
					rmp,
					cut_data,
					fire_plans,
					global_fire_allotment_branching_rules,
					fire,
					new_plan_ix,
				)
			end
		end



		if timing
			details["fire_subproblems"] += time() - t
		end

		@debug "total reduced cost" reduced_cost_sum ub reduced_cost_sum / ub local_gap_rel_tolerance
		# Continue iterating if at least one improving column was added or if
		# the reduced-cost bound is still too weak relative to the incumbent.
		continue_iterating =
			((iteration == 1) || (-reduced_cost_sum / ub > local_gap_rel_tolerance)) && (time_limit > time() - t)

		# if we have not found columns with enough reduced cost
		if ~continue_iterating

			# because of discretization, we do not actually have a guarantee that
			# the deferral variables are 0. set them to 0 here and keep iterating.
			# the deferral variables still help a lot with convergence.
			# alternative is to accept the solution with deferrals, implicitly
			# improving discretization... too complicated but should be a better solution
			# The intent is to stabilize the master problem early on yet return
			# a solution that exactly respects the physical supply constraints.
			if maximum(value.(rmp.deferred_num_crews)) > 1e-5
				for g ∈ 1:num_fires
					for t ∈ 1:num_time_periods
						fix(rmp.deferred_num_crews[g, t], 0, force = true)
					end
				end
				continue_iterating = true
                            @debug "remove deferrals" iteration
			end
		end

        if continue_iterating

			# TODO dual warm start passed in here
			if timing
				t = time()
			end
			# Resolve the restricted master problem with the newly added columns.
			# The solution provides updated dual prices that feed back into the
			# next round of pricing problems.
            optimize!(rmp.model)
			if timing
				details["master_problem"] += (time() - t)
			end

			if termination_status(rmp.model) != MOI.OPTIMAL
                            @debug "non optimal termination status" termination_status(rmp.model)
			end


			## TODO FIX THIS LOGIC AND INFEASIBLE LOGIC
			rmp.termination_status = MOI.ITERATION_LIMIT

			# grab dual values (or farkas vector if infeasible)
			fire_duals = dual.(rmp.plan_per_fire)
			crew_duals = dual.(rmp.route_per_crew)
			linking_duals = dual.(rmp.supply_demand_linking)
			cut_duals = dual.(rmp.gub_cover_cuts)
			global_fire_allot_duals = dual.(rmp.fire_allotment_branches)

			# if the master problem is infeasible
			if (termination_status(rmp.model) == MOI.INFEASIBLE) |
			   (termination_status(rmp.model) == MOI.INFEASIBLE_OR_UNBOUNDED)

				# log this
				@debug "RMP is infeasible, using dual certificate to get dual values" iteration fire_duals crew_duals linking_duals cut_duals dual_status(
					rmp.model,
				)

				# set status
				rmp.termination_status = MOI.INFEASIBLE

				# scale dual values to a feasible dual solution with cost "upper_bound"
				# (linking_duals omitted because 0 RHS)
				# When the continuous relaxation is infeasible we extract a Farkas
				# certificate.  Scaling keeps the certificate comparable to the
				# incumbent upper bound so that the reduced-cost bound remains useful.
				dual_costs = 0
				for ix in eachindex(fire_duals)
					dual_costs +=
						(normalized_rhs(rmp.plan_per_fire[ix]) * fire_duals[ix])
				end

				for ix in eachindex(crew_duals)
					dual_costs +=
						(normalized_rhs(rmp.route_per_crew[ix]) * crew_duals[ix])
				end

				for ix in eachindex(cut_duals)
					dual_costs +=
						(normalized_rhs(rmp.gub_cover_cuts[ix]) * cut_duals[ix])
				end

				for ix in eachindex(global_fire_allot_duals)
					dual_costs +=
						(
							normalized_rhs(rmp.fire_allotment_branches[ix]) *
							global_fire_allot_duals[ix]
						)
				end

				scale = 1.0
				if isfinite(upper_bound) && isfinite(dual_costs) && (abs(dual_costs) > 1e-9)
					scale = upper_bound / dual_costs
					if !isfinite(scale)
						@debug "Dual scaling produced non-finite scale, skipping" upper_bound dual_costs
						scale = 1.0
					end
				else
					@debug "Skipping dual scaling" upper_bound dual_costs
				end
				fire_duals = fire_duals .* scale
				crew_duals = crew_duals .* scale
				linking_duals = linking_duals .* scale
				cut_duals = cut_duals .* scale
				global_fire_allot_duals = global_fire_allot_duals .* scale
			else
				ub = objective_value(rmp.model)
				lb = ub + reduced_cost_sum
				if lb > upper_bound
					@debug "prune by bound" lb upper_bound iteration
					# The reduced-cost lower bound exceeds the best feasible
					# solution, so we can stop exploring this node altogether.
					continue_iterating = false
					rmp.termination_status = MOI.OBJECTIVE_LIMIT
				end
			end

			supp = get_fire_incumbent_weighted_average(
				rmp,
				fire_plans,
				num_fires,
				num_time_periods,
			)
			# Tracking the incumbent suppression pattern helps diagnose slow
			# convergence in logs without dumping whole vectors of duals.
			@debug "progress" iteration linking_duals supp


		# if no new column added, we have proof of optimality
		else
			# re-optimze for JuMP reasons (access attrs) just in case we added a column 
			# but then stopped due to too small reduced cost improvement
			# (Without this call JuMP may prevent us from querying objective/dual info.)
			if timing
				t = time()
			end
			optimize!(rmp.model)

			if (termination_status(rmp.model) == MOI.INFEASIBLE) |
				(termination_status(rmp.model) == MOI.INFEASIBLE_OR_UNBOUNDED)
 
				 # log this
				 @debug "prune by infeasibility"
				 rmp.termination_status = MOI.INFEASIBLE

			else

				rmp.termination_status = MOI.LOCALLY_SOLVED
				@debug "end DCG" iteration termination_status(rmp.model) objective_value(
					rmp.model,
				)

			end

			if timing
				details["master_problem"] += (time() - t)
			end

		end
	end

	details["iteration"] = Float64(iteration)
	return details
end



function define_restricted_master_problem(
	gurobi_env,
	crew_route_data::CrewRouteData,
	crew_avail_ixs::Vector{Vector{Int64}},
	fire_plan_data::FirePlanData,
	fire_avail_ixs::Vector{Vector{Int64}},
	cut_data::CutData,
	fire_allotment_branching_rules::Vector{GlobalFireAllotmentBranchingRule},
	deferral_stabilization::Bool,
	fires_to_ignore::Vector{Int64};
	dual_warm_start::Union{Nothing, DualWarmStart} = nothing,
)
    @debug "Define restricted master problem" fires_to_ignore

	# get dimensions
	# These determine the size of all matrices that follow; the problem can
	# have many crews/fires, so everything is sparse/irregular in practice.
	num_crews, _, num_fires, num_time_periods = size(crew_route_data.fires_fought)

	# inititalze JuMP model
	# We rely on Gurobi for performance and therefore set the tolerances tight
	# enough to make reduced-cost reasoning reliable.
	m = direct_model(Gurobi.Optimizer(gurobi_env))
	set_optimizer_attribute(m, "OutputFlag", 0) # put this first so others don't print
	set_optimizer_attribute(m, "OptimalityTol", 1e-9)
	set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
	set_optimizer_attribute(m, "InfUnbdInfo", 1)
	# set_optimizer_attribute(m, "OutputFlag", 1)


	# decision variables for crew routes and fire plans
	# Each variable represents selecting one pre-generated column.  Column
	# generation will grow the sets `crew_avail_ixs` and `fire_avail_ixs` over time.
	@variable(m, route[c = 1:num_crews, r ∈ crew_avail_ixs[c]] >= 0)
	@variable(m, plan[g = 1:num_fires, p ∈ fire_avail_ixs[g]] >= 0)

	# add deferral stabilization variables
	# These slack variables temporarily allow mismatches between supply and
	# demand to keep the master problem well-behaved until sufficient columns
	# have been generated.  We later force them to zero.
	@variable(
		m,
		deferred_num_crews[
			g = 1:num_fires,
			t = 0:num_time_periods,
		] >= 0
	)
	for g ∈ 1:num_fires
		fix(deferred_num_crews[g, 0], 0, force = true)
	end

	# fix deferral stabilization variables to 0 if not stabilizing
	if ~deferral_stabilization
		for g ∈ 1:num_fires
			for t ∈ 1:num_time_periods
				fix(deferred_num_crews[g, t], 0, force = true)
			end
		end
	end

	# dual warm starts currently handled in column generation; no additional
	# model variables are required here.

	# constraints that you must choose a plan per crew and per fire
	# Each crew selects exactly one route (convex combination) and each fire
	# must have at least one plan in the mix.
	@constraint(m, route_per_crew[c = 1:num_crews],
		sum(route[c, r] for r ∈ crew_avail_ixs[c]) == 1)
	@constraint(m, plan_per_fire[g = 1:num_fires],
		sum(plan[g, p] for p ∈ fire_avail_ixs[g]) >= 1)

	## constraints for cuts

	# get proper coefficients of columns in each cut
	# The cut dictionaries refer to columns by (fire, plan) and (crew, route)
	# indices.  Here we translate them into the indices that JuMP expects and
	# cache them to avoid rebuilding the sparse matrix on every solve.
	cut_ixs = keys(cut_data.cut_dict)

	fire_plan_ixs = Dict()
	fire_plan_coeffs = Dict()
	fire_lookup = cut_data.fire_mp_lookup
	for cut_ix in cut_ixs
		cut_plan_ixs = []
		cut_plan_coeffs = []
		for (fire_ix, coeff) in fire_lookup[cut_ix]
			if fire_ix in eachindex(plan)
				push!(cut_plan_ixs, fire_ix)
				push!(cut_plan_coeffs, -coeff)
			end
		end
		fire_plan_ixs[cut_ix] = cut_plan_ixs
		fire_plan_coeffs[cut_ix] = cut_plan_coeffs
	end

	crew_route_ixs = Dict()
	crew_route_coeffs = Dict()
	crew_lookup = cut_data.crew_mp_lookup
	for cut_ix in cut_ixs
		cut_route_ixs = []
		cut_route_coeffs = []
		for (crew_ix, coeff) in crew_lookup[cut_ix]
			if crew_ix in eachindex(route)
				push!(cut_route_ixs, crew_ix)
				push!(cut_route_coeffs, -coeff)
			end
		end
		crew_route_ixs[cut_ix] = cut_route_ixs
		crew_route_coeffs[cut_ix] = cut_route_coeffs
	end

	# need it to default to SparseAxisArray when empty, maybe there is a better way
	# Each cover cut enforces that if a collection of fires demands crews at a
	# given time, at least one crew route capable of serving them must be chosen.
	@constraint(
		m,
		gub_cover_cuts[
			t = 1:num_time_periods,
			u = 1:10000;
			(t, u) ∈ keys(cut_data.cut_dict),
		],
		sum(
			crew_route_coeffs[t, u][i] * route[crew_route_ixs[t, u][i]] for
			i ∈ eachindex(crew_route_coeffs[t, u])
		) +
		sum(
			fire_plan_coeffs[t, u][i] * plan[fire_plan_ixs[t, u][i]] for
			i ∈ eachindex(fire_plan_coeffs[t, u])
		) >=
		-cut_data.cut_dict[(t, u)].rhs
	)

	# container for fire allotment branching rules
	# The branching rules resemble aggregate knapsack constraints that we
	# populate later when applying branch decisions.
	@constraint(
		m,
		fire_allotment_branches[eachindex(fire_allotment_branching_rules)],
		0 >= 0
	)

	# set constraint RHS and plan coeffs
	for ix in eachindex(fire_allotment_branching_rules)
		rule = fire_allotment_branching_rules[ix]

		# (following >= convention)
		# Branching rules can be >= or <=; we convert them into the >= form
		# expected by JuMP by flipping both sides when needed.
		sign = 2 * Int(rule.geq_flag) - 1

		# possible RHS values are <= 0, >= 1
		set_normalized_rhs(fire_allotment_branches[ix], Int(rule.geq_flag))

		# get coeff for any fire plan exceeding the allotment
		for (fire_plan_ix, excess) ∈ rule.mp_lookup
			if fire_plan_ix ∈ eachindex(plan)
				set_normalized_coefficient(
					fire_allotment_branches[ix],
					plan[fire_plan_ix],
					sign * excess,
				)
			end
		end

	end

	# linking constraint
	# These constraints are the heart of the master problem: they ensure that
	# the sum of crews allocated to a fire at a time (from crew columns) meets
	# the demand requested by the chosen fire plans.
	@constraint(m, linking[g = 1:num_fires, t = 1:num_time_periods],

		# crews at fire
		sum(
			route[c, r] * crew_route_data.fires_fought[c, r, g, t]
			for c ∈ 1:num_crews, r ∈ crew_avail_ixs[c]
		) +
		deferred_num_crews[g, t-1] - deferred_num_crews[g, t]
		>=

		# crews suppressing
		sum(
			plan[g, p] * fire_plan_data.crews_present[g, p, t]
			for p ∈ fire_avail_ixs[g]
		))


	@objective(m, Min,

		# route costs
		sum(
			route[c, r] * crew_route_data.route_costs[c, r]
			for c ∈ 1:num_crews, r ∈ crew_avail_ixs[c]
		)
		+

		# suppression plan costs
		sum(
			plan[g, p] * fire_plan_data.plan_costs[g, p]
			for g ∈ 1:num_fires, p ∈ fire_avail_ixs[g]
		)

		- 

		# ignore fires that are not started yet
		# (Their plans may be present for bookkeeping, but we do not want them
		# to bias the objective until they are activated elsewhere.)
		sum(
			plan[g, p] * fire_plan_data.plan_costs[g, p] 
			for g ∈ fires_to_ignore, p ∈ fire_avail_ixs[g]
		)
	)

	return RestrictedMasterProblem(
		m,
		crew_avail_ixs,
		fire_avail_ixs,
		route,
		plan,
		deferred_num_crews,
		route_per_crew,
		plan_per_fire,
		linking,
		gub_cover_cuts,
		fire_allotment_branches,
		MOI.OPTIMIZE_NOT_CALLED,
	)

end

"""
Add a new suppression plan to the FirePlanData structure for a specific fire.

# Arguments

  - `plan_data::FirePlanData` (modified): The FirePlanData structure to which the column will be added.
  - `fire::Int64`: The index of the fire for which the plan is being added.
  - `cost::Float64`: The cost associated with the plan for the given fire.
  - `crew_demands::Vector{Int64}`: The crews required for the plan at each time.

# Returns

  - `ix::Int64`: The index of the newly added column in the plan data for the specified fire.

# Description

  - Increments the count of plans for the specified fire.
  - Appends the plan cost to the plan costs matrix for the specified fire.
  - Appends the crew demands to the crews present matrix for the specified fire.
"""
function add_column_to_plan_data!(
	plan_data::FirePlanData,
	fire::Int64,
	cost::Float64,
	crew_demands::Vector{Int64},
	arcs_used::Vector{Int64},
)
	# add 1 to number of plans for this fire, store the index
	plan_data.plans_per_fire[fire] += 1
	ix = plan_data.plans_per_fire[fire]

	# append the route cost
	plan_data.plan_costs[fire, ix] = cost

	# append the fires fought
	plan_data.crews_present[fire, ix, :] = crew_demands

	# append the arcs used
	# Storing the actual dynamic-programming arcs allows reconstruction of
	# the plan later for visualization or branching.
	plan_data.arcs_used[fire, ix] = arcs_used

	return ix

end

"""
Add a new route to the CrewRouteData structure for a specific crew.

# Arguments

  - `route_data::CrewRouteData` (modified): The CrewRouteData structure to which the column will be added.
  - `crew::Int64`: The index of the crew for which the route data is being added.
  - `cost::Float64`: The cost associated with the route for the given crew.
  - `fires_fought::BitArray{2}`: A 2-dimensional BitArray representing the fires fought by the crew at each time for the given route.

# Returns

  - `ix::Int64`: The index of the newly added column in the route data for the specified crew.

# Description

  - Increments the count of routes for the specified crew.
  - Appends the route cost to the route costs matrix for the specified crew.
  - Appends the fires fought to the fires_fought array for the specified crew.
"""
function add_column_to_route_data!(
	route_data::CrewRouteData,
	crew::Int64,
	cost::Float64,
	fires_fought::BitArray{2},
	arcs_used::Vector{Int64},
)

	# add 1 to number of routes for this crew, store the index
	route_data.routes_per_crew[crew] += 1
	ix = route_data.routes_per_crew[crew]

	# append the route cost
	route_data.route_costs[crew, ix] = cost

	# append the fires fought
	route_data.fires_fought[crew, ix, :, :] = fires_fought

	# append the arcs used
	# Keeping the arcs lets us re-solve or branch on individual transitions.
	route_data.arcs_used[crew, ix] = arcs_used

	return ix

end

"""
Add a new column to the RestrictedMasterProblem representing a crew route.

# Arguments

  - `rmp::RestrictedMasterProblem` (modified): The RestrictedMasterProblem object to which the column will be added.
  - `cut_data::CutData` (modified): The CutData object containing information about the cuts.
  - `crew_routes::CrewRouteData`: The CrewRouteData object containing information about crew routes.
  - `crew::Int64`: The index of the crew for which the route is being added.
  - `ix::Int64`: The index (for this crew) of the route being added.

# Description

  - Defines a variable representing the new route in the RestrictedMasterProblem.
  - Updates the coefficient of this route in the objective function of the RestrictedMasterProblem.
  - Updates coefficient of this route in all relevant constraints.
"""
function add_column_to_master_problem!!(
	rmp::RestrictedMasterProblem,
	cut_data::CutData,
	crew_routes::CrewRouteData,
	crew::Int64,
	ix::Int64,
)

	# define variable
	rmp.routes[crew, ix] =
		@variable(rmp.model, base_name = "route[$crew,$ix]", lower_bound = 0)

	# update index lookup
	push!(rmp.crew_column_ixs[crew], ix)

	# update coefficient in objective
	set_objective_coefficient(
		rmp.model,
		rmp.routes[crew, ix],
		crew_routes.route_costs[crew, ix],
	)

	## update coefficient in constraints

	# route per crew
	set_normalized_coefficient(rmp.route_per_crew[crew], rmp.routes[crew, ix], 1)

	# supply demand linking
	# This is a 4-D array assignment: (fire, time) pairs get +1 when the crew
	# attends the fire during that time in the route.
	set_normalized_coefficient.(
		rmp.supply_demand_linking,
		rmp.routes[crew, ix],
		crew_routes.fires_fought[crew, ix, :, :],
	)

	# cuts

	# for each cut in the cut data
	for (cut_ix, cut) ∈ cut_data.cut_dict

		# if this crew is involved in the cut
		if cut.crew_coeffs[crew] > 1e-20

			# get the fires not suppressed at the given time in order to enter cut
			fires = [i for i in keys(cut.fire_coeffs)]

			# if these fires are not suppressed at the given time
			if maximum(crew_routes.fires_fought[crew, ix, fires, cut.time_ix]) == 0

				# set the coefficient of this route in this cut to -1 (negated because >=)
				set_normalized_coefficient(
					rmp.gub_cover_cuts[cut_ix],
					rmp.routes[crew, ix],
					-1,
				)

				# add the plan to the cut mp lookup
				# (Used when dynamically updating the cut RHS later.)
				cut_data.crew_mp_lookup[cut_ix][(crew, ix)] = 1

			end

		end
	end

end

"""
Add a new column to the RestrictedMasterProblem representing a fire plan.

# Arguments

  - `rmp::RestrictedMasterProblem` (modified): The RestrictedMasterProblem object to which the column will be added.
  - `cut_data::CutData` (modified): The CutData object containing information about the cuts.
  - `fire_plans::FirePlanData`: The FirePlanData object containing information about fire plans.
  - `fire_allotment_branching_rules::Vector{GlobalFireAllotmentBranchingRule}`: The vector of GlobalFireAllotmentBranchingRule containing information about global fire allotment branching rules.
  - `fire::Int64`: The index of the fire for which the plan is being added.
  - `ix::Int64`: The index (for this fire)of the plan being added.

# Description

  - Defines a variable representing the new plan in the RestrictedMasterProblem.
  - Updates the coefficient in the objective function of the RestrictedMasterProblem.
  - Updates coefficient of this plan in all relevant constraints.
"""
function add_column_to_master_problem!!(
	rmp::RestrictedMasterProblem,
	cut_data::CutData,
	fire_plans::FirePlanData,
	fire_allotment_branching_rules::Vector{GlobalFireAllotmentBranchingRule},
	fire::Int64,
	ix::Int64,
)

	# define variable
	rmp.plans[fire, ix] =
		@variable(rmp.model, base_name = "plan[$fire,$ix]", lower_bound = 0)

	# update index lookup
	push!(rmp.fire_column_ixs[fire], ix)

	# update coefficient in objective
	set_objective_coefficient(
		rmp.model,
		rmp.plans[fire, ix],
		fire_plans.plan_costs[fire, ix],
	)

	# update coefficient in constraints

	# plan per fire
	set_normalized_coefficient(rmp.plan_per_fire[fire], rmp.plans[fire, ix], 1)

	# supply demand linking
	# Fire plans request crews, so they enter the linking constraints with a
	# negative sign to offset the positive crew contributions.
	set_normalized_coefficient.(
		rmp.supply_demand_linking[fire, :],
		rmp.plans[fire, ix],
		-fire_plans.crews_present[fire, ix, :],
	)

	# cuts

	# for each cut in the cut data
	for (cut_ix, cut) ∈ cut_data.cut_dict

		# if this fire is involved in the cut
		if fire ∈ keys(cut.fire_coeffs)

			# update mp lookup and see if this plan has >0 coeff in the cut
			# (The helper inspects the plan structure and determines whether it
			# violates the coverage minimum encoded by the cut.)
			plan_in_cut = update_cut_fire_mp_lookup!(
				cut_data.fire_mp_lookup[cut_ix],
				cut,
				fire_plans,
				fire,
				ix,
			)

			if plan_in_cut
				# set the coefficient of this plan in this cut to the coefficient of the cut (negated because >=)
				set_normalized_coefficient(
					rmp.gub_cover_cuts[cut_ix],
					rmp.plans[fire, ix],
					-cut_data.fire_mp_lookup[cut_ix][(fire, ix)],
				)
			end

		end
	end

	# global fire allotment branching rule, add new plan to mp_lookup
	# Branch-and-price adds aggregate branching constraints outside the normal
	# fire linking constraints.  Whenever we introduce a new plan we need to
	# inform those constraints so they can price correctly.
	for rule_ix in eachindex(rmp.fire_allotment_branches)
		rule = fire_allotment_branching_rules[rule_ix]
		sign = 2 * Int(rule.geq_flag) - 1

		# if fire plan exceeds the allotment
		if (fire, ix) ∈ keys(rule.mp_lookup)
			set_normalized_coefficient(
				rmp.fire_allotment_branches[rule_ix],
				rmp.plans[(fire, ix)],
				sign * rule.mp_lookup[(fire, ix)],
			)
		end
	end
end

function get_fire_incumbent_weighted_average(
	rmp::RestrictedMasterProblem,
	fire_plans::FirePlanData,
	num_fires::Int,
	num_time_periods::Int,
)

	# Weighted average of all fire plans that currently have positive weight in
	# the master problem solution.  The result is the expected number of crews
	# requested by each fire at every time.
	fire_allotment = zeros(num_fires, num_time_periods)
	for ix in eachindex(rmp.plans)
		if value(rmp.plans[ix]) > 0
			fire_allotment[ix[1], :] +=
				fire_plans.crews_present[ix..., :] * value(rmp.plans[ix])
		end
	end

	return fire_allotment

end

function get_crew_incumbent_weighted_average(
	rmp::RestrictedMasterProblem,
	crew_routes::CrewRouteData,
)

	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	# Similar to `get_fire_incumbent_weighted_average` but preserves the crew
	# dimension so we can see which team is expected to be at which fire.
	crew_allotment = zeros(Float64, num_crews, num_fires, num_time_periods)
	for crew in 1:num_crews
		for col in rmp.crew_column_ixs[crew]
			weight = value(rmp.routes[(crew, col)])
			if weight > 0
				crew_allotment[crew, :, :] +=
					crew_routes.fires_fought[crew, col, :, :] * weight
			end
		end
	end

	return crew_allotment
end

function get_fire_and_crew_incumbent_weighted_average(
	rmp::RestrictedMasterProblem,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
)

	# get problem dimensions
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	# Convenience wrapper that returns both perspectives at once so callers do
	# not need to recompute the expensive weighted averages twice.
	fire_allotment = get_fire_incumbent_weighted_average(
		rmp,
		fire_plans,
		num_fires,
		num_time_periods,
	)
	crew_allotment = get_crew_incumbent_weighted_average(rmp, crew_routes)
	return fire_allotment, crew_allotment
end

function get_cost_due_to_fires_and_crews(
	solved_rmp::RestrictedMasterProblem,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
)

	# get problem dimensions
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	# get the cost due to fires
	# Iterate over all plan variables because sparsity patterns can change as
	# new columns are added; relying on `fire_column_ixs` would miss new
	# entries that JuMP already created.
	fire_cost = 0
	for ix in eachindex(solved_rmp.plans)
		if value(solved_rmp.plans[ix]) > 0
			f = ix[1]
			plan = ix[2]
			fire_cost += value(solved_rmp.plans[ix]) * fire_plans.plan_costs[f, plan]
		end
	end

	# get the cost due to crews
	# Same idea as above, but over the crew route columns.
	crew_cost = 0
	for ix in eachindex(solved_rmp.routes)
		if value(solved_rmp.routes[ix]) > 0
			c = ix[1]
			route = ix[2]
			crew_cost += value(solved_rmp.routes[ix]) * crew_routes.route_costs[c, route]
		end
	end
			

	return fire_cost, crew_cost
end

function get_fire_and_crew_arcs_used(
	solved_rmp::RestrictedMasterProblem,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
)

	# get problem dimensions
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	# get the arcs used by fires
	# Each arc list contains the dynamic-programming path corresponding to the
	# integer solution.  Having the arcs makes it easier to export maps later.
	fire_arcs_used = Vector{Vector{Int64}}(undef, num_fires)
	for f in 1:num_fires
		fire_arcs_used[f] = Int64[]
	end
	for ix in eachindex(solved_rmp.plans)
		if value(solved_rmp.plans[ix]) > 0.99
			f = ix[1]
			plan = ix[2]
			# Guard against uninitialized plan columns in FirePlanData
			if (plan <= fire_plans.plans_per_fire[f]) && isassigned(fire_plans.arcs_used, f, plan)
				fire_arcs_used[f] = fire_plans.arcs_used[f, plan]
			else
				fire_arcs_used[f] = Int64[]
			end
		end
	end

	# get the arcs used by crews
	# Same story for crews; retain an empty vector for crews not selected.
	crew_arcs_used = Vector{Vector{Int64}}(undef, num_crews)
	for c in 1:num_crews
		crew_arcs_used[c] = Int64[]
	end
	for ix in eachindex(solved_rmp.routes)
		if value(solved_rmp.routes[ix]) > 0.99
			c = ix[1]
			route = ix[2]
			# Guard against uninitialized route columns in CrewRouteData
			if (route <= crew_routes.routes_per_crew[c]) && isassigned(crew_routes.arcs_used, c, route)
				crew_arcs_used[c] = crew_routes.arcs_used[c, route]
			else
				crew_arcs_used[c] = Int64[]
			end
		end
	end

	return fire_arcs_used, crew_arcs_used
end
