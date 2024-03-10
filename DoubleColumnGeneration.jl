include("CommonStructs.jl")
include("Subproblems.jl")
include("BranchingRules.jl")

using JuMP


@kwdef mutable struct DualWarmStart

	linking_values::Matrix{Float64}
	const strategy::String = "global"
	epsilon::Float64 = 0.001

end

function define_restricted_master_problem(
	gurobi_env,
	crew_route_data::CrewRouteData,
	crew_avail_ixs::Vector{Vector{Int64}},
	fire_plan_data::FirePlanData,
	fire_avail_ixs::Vector{Vector{Int64}},
	cut_data::GUBCoverCutData,
	fire_allotment_branching_rules::Vector{GlobalFireAllotmentBranchingRule},
	dual_warm_start::Union{Nothing, DualWarmStart} = nothing,
)

	# get dimensions
	num_crews, _, num_fires, num_time_periods = size(crew_route_data.fires_fought)

	# inititalze JuMP model
	m = Model(() -> Gurobi.Optimizer(gurobi_env))
	set_optimizer_attribute(m, "OptimalityTol", 1e-9)
	set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
	set_optimizer_attribute(m, "OutputFlag", 0)
	set_optimizer_attribute(m, "InfUnbdInfo", 1)


	# decision variables for crew routes and fire plans
	@variable(m, route[c = 1:num_crews, r ∈ crew_avail_ixs[c]] >= 0)
	@variable(m, plan[g = 1:num_fires, p ∈ fire_avail_ixs[g]] >= 0)

	if ~isnothing(dual_warm_start)

		error("Not implemented")
		# dual stabilization variables
		@variable(m, delta_plus[g = 1:num_fires, t = 1:num_time_periods] >= 0)
		@variable(m, delta_minus[g = 1:num_fires, t = 1:num_time_periods] >= 0)

	end

	# constraints that you must choose a plan per crew and per fire
	@constraint(m, route_per_crew[c = 1:num_crews],
		sum(route[c, r] for r ∈ crew_avail_ixs[c]) == 1)
	@constraint(m, plan_per_fire[g = 1:num_fires],
		sum(plan[g, p] for p ∈ fire_avail_ixs[g]) >= 1)

	# container for GUB cover cuts
	# need it to default to SparseAxisArray when empty, maybe there is a better way
	@constraint(
		m,
		gub_cover_cuts[
			t = 1:num_time_periods,
			u = 1:1000;
			(t, u) ∈ keys(cut_data.cut_dict),
		],
		0 >= 0
	)

	# set cut RHS (negating because of >= convention)
	for ix in keys(cut_data.cut_dict)
		set_normalized_rhs(gub_cover_cuts[ix], -cut_data.cut_dict[ix].rhs)
	end

	# get proper coefficients of columns in gub_cover_cuts
	fire_lookup = cut_data.fire_mp_lookup
	for cut_ix in eachindex(fire_lookup)
		for (fire_ix, coeff) in fire_lookup[cut_ix]
			if fire_ix in eachindex(plan)
				set_normalized_coefficient(
					gub_cover_cuts[cut_ix],
					plan[fire_ix],
					-coeff,
				)
			end
		end
	end
	crew_lookup = cut_data.crew_mp_lookup
	for cut_ix in eachindex(crew_lookup)
		for (crew_ix, coeff) in crew_lookup[cut_ix]
			if crew_ix in eachindex(route)
				set_normalized_coefficient(
					gub_cover_cuts[cut_ix],
					route[crew_ix],
					-coeff,
				)
			end
		end
	end


	# container for fire allotment branching rules
	@constraint(
		m,
		fire_allotment_branches[eachindex(fire_allotment_branching_rules)],
		0 >= 0
	)

	# set constraint RHS and plan coeffs
	for ix in eachindex(fire_allotment_branching_rules)
		rule = fire_allotment_branching_rules[ix]

		# (following >= convention)
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
	if isnothing(dual_warm_start)

		# for each fire and time period
		@constraint(m, linking[g = 1:num_fires, t = 1:num_time_periods],

			# crews at fire
			sum(
				route[c, r] * crew_route_data.fires_fought[c, r, g, t]
				for c ∈ 1:num_crews, r ∈ crew_avail_ixs[c]
			)
			>=

			# crews suppressing
			sum(
				plan[g, p] * fire_plan_data.crews_present[g, p, t]
				for p ∈ fire_avail_ixs[g]
			))

	elseif dual_warm_start.strategy == "global"

		error("Not implemented")

		# get expected dual value ratios
		ratios = dual_warm_start.linking_values
		ratios = ratios / sum(ratios)

		@constraint(m, linking[g = 1:num_fires, t = 1:num_time_periods],

			# crews at fire
			sum(
				route[c, r] * crew_route_data.fires_fought[c, r, g, t]
				for c ∈ 1:num_crews, r ∈ crew_avail_ixs[c]
			)
			+

			# perturbation
			delta_plus[g, t] - delta_minus[g, t] -
			sum(ratios .* delta_plus) + sum(ratios .* delta_minus)
			>=

			# crews suppressing
			sum(
				plan[g, p] * fire_plan_data.crews_present[g, p, t]
				for p ∈ fire_avail_ixs[g]
			))

		# this constrant neutralizes the perturbation, will be presolved away if RHS is 0
		# but raising the RHS slightly above 0 allows the perturbation
		@constraint(m, perturb[g = 1:num_fires, t = 1:num_time_periods],
			delta_plus[g, t] + delta_minus[g, t] <= 0)
	else
		error("Dual stabilization type not implemented")
	end


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
	)

	return RestrictedMasterProblem(
		m,
		route,
		plan,
		route_per_crew,
		plan_per_fire,
		linking,
		gub_cover_cuts,
		fire_allotment_branches,
		MOI.OPTIMIZE_NOT_CALLED,
	)

end


function add_column_to_plan_data!(
	plan_data::FirePlanData,
	fire::Int64,
	cost::Float64,
	crew_demands::Vector{Int64},
)
	# add 1 to number of plans for this fire, store the index
	plan_data.plans_per_fire[fire] += 1
	ix = plan_data.plans_per_fire[fire]

	# append the route cost
	plan_data.plan_costs[fire, ix] = cost

	# append the fires fought
	plan_data.crews_present[fire, ix, :] = crew_demands

	return ix

end

function add_column_to_route_data!(
	route_data::CrewRouteData,
	crew::Int64,
	cost::Float64,
	fires_fought::BitArray{2},
)
	# add 1 to number of routes for this crew, store the index
	route_data.routes_per_crew[crew] += 1
	ix = route_data.routes_per_crew[crew]

	# append the route cost
	route_data.route_costs[crew, ix] = cost

	# append the fires fought
	route_data.fires_fought[crew, ix, :, :] = fires_fought

	return ix

end

function add_column_to_master_problem!!(
	rmp::RestrictedMasterProblem,
	cut_data::GUBCoverCutData,
	crew_routes::CrewRouteData,
	crew::Int64,
	ix::Int64,
)

	# define variable
	rmp.routes[crew, ix] =
		@variable(rmp.model, base_name = "route[$crew,$ix]", lower_bound = 0)

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
	set_normalized_coefficient.(
		rmp.supply_demand_linking,
		rmp.routes[crew, ix],
		crew_routes.fires_fought[crew, ix, :, :],
	)

	# cuts

	# for each cut in the cut data
	for (cut_ix, cut) ∈ cut_data.cut_dict

		# if this crew is involved in the cut
		if crew ∈ cut.inactive_crews

			# get the fires not suppressed at the given time in order to enter cut
			fires = [i for i in keys(cut.fire_lower_bounds)]

			# if these fires are not suppressed at the given time
			if maximum(crew_routes.fires_fought[crew, ix, fires, cut.time_ix]) == 0

				# set the coefficient of this route in this cut to -1 (negated because >=)
				set_normalized_coefficient(
					rmp.gub_cover_cuts[cut_ix],
					rmp.routes[crew, ix],
					-1,
				)

				# add the plan to the cut mp lookup
				cut_data.crew_mp_lookup[cut_ix][(crew, ix)] = 1

			end

		end
	end

end

function add_column_to_master_problem!!(
	rmp::RestrictedMasterProblem,
	cut_data::GUBCoverCutData,
	fire_plans::FirePlanData,
	fire_allotment_branching_rules::Vector{GlobalFireAllotmentBranchingRule},
	fire::Int64,
	ix::Int64,
)

	# define variable
	rmp.plans[fire, ix] =
		@variable(rmp.model, base_name = "plan[$fire,$ix]", lower_bound = 0)

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
	set_normalized_coefficient.(
		rmp.supply_demand_linking[fire, :],
		rmp.plans[fire, ix],
		-fire_plans.crews_present[fire, ix, :],
	)

	# cuts

	# for each cut in the cut data
	for (cut_ix, cut) ∈ cut_data.cut_dict

		# if this fire is involved in the cut
		if fire ∈ keys(cut.fire_lower_bounds)

			# update mp lookup and see if this plan has >0 coeff in the cut
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


function get_fire_and_crew_incumbent_weighted_average(
	rmp::RestrictedMasterProblem,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
)

	# get problem dimensions
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	fire_allotment = zeros(num_fires, num_time_periods)
	for ix in eachindex(rmp.plans)
		if value(rmp.plans[ix]) > 0
			fire_allotment[ix[1], :] +=
				fire_plans.crews_present[ix..., :] * value(rmp.plans[ix])
		end
	end

	crew_allotment = zeros(num_crews, num_fires, num_time_periods)
	for ix in eachindex(rmp.routes)
		if value(rmp.routes[ix]) > 0
			crew_allotment[ix[1], :, :] +=
				crew_routes.fires_fought[ix..., :, :] * value(rmp.routes[ix])
		end
	end

	fire_allotment, crew_allotment
end

function double_column_generation!(
	rmp::RestrictedMasterProblem,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	crew_branching_rules::Vector{CrewSupplyBranchingRule},
	fire_branching_rules::Vector{FireDemandBranchingRule},
	global_fire_allotment_branching_rules::Vector{GlobalFireAllotmentBranchingRule},
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	cut_data::GUBCoverCutData;
	upper_bound::Float64,
	improving_column_abs_tolerance::Float64 = 1e-4,
	local_gap_rel_tolerance::Float64 = 1e-9)

	# gather global information
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	if rmp.termination_status == MOI.OPTIMIZE_NOT_CALLED
		# initialize with an (infeasible) dual solution that will suppress minimally
		fire_duals = zeros(num_fires) .+ Inf
		crew_duals = zeros(num_crews)
		linking_duals = zeros(num_fires, num_time_periods) .+ 1e30
		cut_duals = normalized_rhs.(rmp.gub_cover_cuts) .* 0
		global_fire_allot_duals = normalized_rhs.(rmp.fire_allotment_branches) .* 0
	else
		fire_duals = dual.(rmp.plan_per_fire)
		crew_duals = dual.(rmp.route_per_crew)
		linking_duals = dual.(rmp.supply_demand_linking)
		cut_duals = dual.(rmp.gub_cover_cuts)
		global_fire_allot_duals = dual.(rmp.fire_allotment_branches)
	end

	ub = copy(upper_bound)

	# initialize column generation loop
	continue_iterating::Bool = true
	iteration = 0

	while (continue_iterating & (iteration < 200))
		if iteration == 199
			@info "iteration limit hit" iteration
		end

		iteration += 1
		reduced_cost_sum = 0

		# Not for any good reason, the crew subproblems all access the 
		# same set of arcs in matrix form, and each runs its subproblem 
		# on a subset of the arcs. This means that dual-adjusting arc 
		# costs happens once only. In contrast, in the fire subproblems
		# this happens inside the loop for each fire. 
		# (TODO: see which is faster in Julia)

		subproblem = crew_subproblems[1]

		# generate the local costs of the arcs
		# TODO NEEDS refactor on prohibited_arcs
		rel_costs, prohibited_arcs = get_adjusted_crew_arc_costs(
			subproblem.long_arcs,
			linking_duals,
			crew_branching_rules,
		)

		arc_costs = rel_costs .+ subproblem.arc_costs

		# for each crew
		for crew in 1:num_crews

			# adjust the arc costs for the cuts
			cut_adjusted_arc_costs = cut_adjust_arc_costs(
				arc_costs,
				cut_data.crew_sp_lookup[crew],
				cut_duals,
			)

			# extract the subproblem
			subproblem = crew_subproblems[crew]

			# TODO grab the prohibited arcs belonging to this crew only 
			crew_prohibited_arcs = prohibited_arcs

			# solve the subproblem
			objective, arcs_used = crew_dp_subproblem(
				subproblem.wide_arcs,
				cut_adjusted_arc_costs,
				crew_prohibited_arcs,
				subproblem.state_in_arcs,
			)

			# adjust the objective for the cuts (we gave -1 if allotment not broken, give +1 here)
			for (ix, cut) in cut_data.cut_dict
				if crew ∈ cut.inactive_crews
					objective += cut_duals[ix]
				end
			end

			# if there is an improving route
			if objective < crew_duals[crew] - improving_column_abs_tolerance

				reduced_cost_sum += (objective - crew_duals[crew])

				# get the real cost, unadjusted for duals
				cost = sum(subproblem.arc_costs[arcs_used])

				# get the indicator matrix of fires fought at each time
				fires_fought = get_fires_fought(
					subproblem.wide_arcs,
					arcs_used,
					(num_fires, num_time_periods),
				)

				# add the route to the routes
				new_route_ix =
					add_column_to_route_data!(crew_routes, crew, cost, fires_fought)

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

		# for each fire
		for fire in 1:num_fires

			# extract the subproblem
			subproblem = fire_subproblems[fire]

			# generate the local costs of the arcs
			rel_costs, prohibited_arcs = get_adjusted_fire_arc_costs(
				subproblem.long_arcs,
				linking_duals[fire, :],
				[rule for rule ∈ fire_branching_rules if rule.fire_ix == fire],
			)

			arc_costs = rel_costs .+ subproblem.arc_costs

			# adjust the arc costs for the cuts
			cut_adjusted_arc_costs = cut_adjust_arc_costs(
				arc_costs,
				cut_data.fire_sp_lookup[fire],
				cut_duals,
			)

			# for each branching rule
			for ix in eachindex(global_fire_allotment_branching_rules)

				rule = global_fire_allotment_branching_rules[ix]

				# if this is a <= 0 rule, we have more prohibited arcs
				if ~rule.geq_flag
					for arc_ix ∈ rule.fire_sp_arc_lookup[fire]
						prohibited_arcs[arc_ix] = true
					end
				end

				# we need to do a proper dual adjustment
				cut_adjusted_arc_costs = adjust_fire_sp_arc_costs(
					rule,
					fire,
					cut_adjusted_arc_costs,
					global_fire_allot_duals[ix],
				)
			end

			# solve the subproblem
			objective, arcs_used = fire_dp_subproblem(
				subproblem.wide_arcs,
				cut_adjusted_arc_costs,
				prohibited_arcs,
				subproblem.state_in_arcs,
			)

			# if there is an improving plan
			if objective < fire_duals[fire] - improving_column_abs_tolerance

				reduced_cost_sum += (objective - fire_duals[fire])

				# get the real cost, unadjusted for duals
				cost = sum(subproblem.arc_costs[arcs_used])

				# get the vector of crew demands at each time
				crew_demands = get_crew_demands(
					subproblem.wide_arcs,
					arcs_used,
					num_time_periods,
				)

				# add the plan to the plans
				new_plan_ix =
					add_column_to_plan_data!(fire_plans, fire, cost, crew_demands)

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

		continue_iterating =
			(iteration == 1) || (-reduced_cost_sum / ub > local_gap_rel_tolerance)
		if continue_iterating

			# TODO dual warm start passed in here
			optimize!(rmp.model)

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
				@info "RMP is infeasible, using dual certificate to get dual values" iteration fire_duals crew_duals linking_duals cut_duals dual_status(
					rmp.model,
				)

				# set status
				rmp.termination_status = MOI.INFEASIBLE

				# scale dual values to a feasible dual solution with cost "upper_bound"
				# (linking_duals omitted because 0 RHS)
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

				scale = upper_bound / dual_costs
				fire_duals = fire_duals .* scale
				crew_duals = crew_duals .* scale
				linking_duals = linking_duals .* scale
				cut_duals = cut_duals .* scale
			else
				ub = objective_value(rmp.model)
				lb = ub + reduced_cost_sum
				if lb > upper_bound
					@info "prune by bound" lb upper_bound iteration
					continue_iterating = false
					rmp.termination_status = MOI.OBJECTIVE_LIMIT
				end
			end

			# if no new column added, we have proof of optimality
		else
			rmp.termination_status = MOI.LOCALLY_SOLVED
			@info "RMP stats with no more columns found" iteration objective_value(
				rmp.model,
			)
			fire_allots, crew_allots = get_fire_and_crew_incumbent_weighted_average(
				rmp,
				crew_routes,
				fire_plans,
			)
		end
	end
end

