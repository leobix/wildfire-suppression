include("CommonStructs.jl")
include("BranchingRules.jl")
include("GUBKnapsackCoverCuts.jl")

using Gurobi, Statistics

mutable struct MutGRBsvec
    len::Cint
    ind::Ptr{Cint}
    val::Ptr{Cdouble}
end

mutable struct BranchAndBoundNode

	const ix::Int64
	const parent::Union{Nothing, BranchAndBoundNode}
	const new_crew_branching_rules::Vector{CrewSupplyBranchingRule}
	const new_fire_branching_rules::Vector{FireDemandBranchingRule}
	const new_global_fire_allotment_branching_rules::Vector{
		GlobalFireAllotmentBranchingRule,
	}
	const cut_data::GUBCutData
	children::Vector{BranchAndBoundNode}
	l_bound::Float64
	u_bound::Float64
	master_problem::Union{Nothing, RestrictedMasterProblem}
	heuristic_found_master_problem::Union{Nothing, RestrictedMasterProblem}
	feasible::Union{Nothing, Bool}
end

function BranchAndBoundNode(
	;
	ix::Int64,
	parent::Union{Nothing, BranchAndBoundNode},
	cut_data::GUBCutData,
	new_crew_branching_rules::Vector{CrewSupplyBranchingRule} = CrewSupplyBranchingRule[],
	new_fire_branching_rules::Vector{FireDemandBranchingRule} = FireDemandBranchingRule[],
	new_global_fire_allotment_branching_rules::Vector{
		GlobalFireAllotmentBranchingRule,
	} = GlobalFireAllotmentBranchingRule[],
	children::Vector{BranchAndBoundNode} = BranchAndBoundNode[],
	l_bound::Float64 = -Inf,
	u_bound::Float64 = Inf,
	master_problem::Union{Nothing, RestrictedMasterProblem} = nothing,
	heuristic_found_master_problem::Union{Nothing, RestrictedMasterProblem} = nothing,
	feasible::Union{Nothing, Bool} = nothing)

	return BranchAndBoundNode(
		ix,
		parent,
		new_crew_branching_rules,
		new_fire_branching_rules,
		new_global_fire_allotment_branching_rules,
		cut_data,
		children,
		l_bound,
		u_bound,
		master_problem,
		heuristic_found_master_problem,
		feasible,
	)
end

###############

function max_variance_natural_variable(
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	route_values::JuMP.Containers.SparseAxisArray,
	plan_values::JuMP.Containers.SparseAxisArray,
)
	# gather global information
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	# calculate variance of A_{cgt}, crew c suppressing fire g at time t
	crew_means = zeros(Float64, (num_crews, num_fires, num_time_periods))
	for ix ∈ eachindex(route_values)
		crew = ix[1]
		route = ix[2]
		crew_means[crew, :, :] +=
			route_values[ix] * crew_routes.fires_fought[crew, route, :, :]
	end
	crew_variances = crew_means .* (1 .- crew_means)

	# calculate variance of B_{gt}, demand at fire g at time t
	fire_means = zeros(Float64, (num_fires, num_time_periods))
	fire_sq_means = zeros(Float64, (num_fires, num_time_periods))
	for ix ∈ eachindex(plan_values)
		fire = ix[1]
		plan = ix[2]
		fire_means[fire, :] +=
			plan_values[ix] * fire_plans.crews_present[fire, plan, :]
		fire_sq_means[fire, :] +=
			plan_values[ix] * (fire_plans.crews_present[fire, plan, :] .^ 2)
	end
	fire_variances = fire_sq_means - (fire_means .^ 2)
	@debug "Variances" fire_variances # crew_variances
	# get the max variance for each natural variable type
	crew_max_var, crew_max_ix = findmax(crew_variances)
	fire_max_var, fire_max_ix = findmax(fire_variances)

	# return the info needed to create the branching rule
	if fire_max_var > crew_max_var
		return "fire", fire_max_ix, fire_max_var, fire_means[fire_max_ix]
	else
		return "crew", crew_max_ix, crew_max_var, crew_means[crew_max_ix]
	end

end


function apply_branching_rule(
	crew_avail_ixs::Vector{Vector{Int64}},
	crew_routes::CrewRouteData,
	branching_rule::CrewSupplyBranchingRule)

	output = [Int64[] for fire in 1:size(crew_avail_ixs)[1]]
	for crew in 1:size(crew_avail_ixs)[1]
		ixs = copy(crew_avail_ixs[crew])
		if branching_rule.crew_ix == crew
			ixs = [
				i for i in ixs if satisfies_branching_rule(
					branching_rule,
					crew_routes.fires_fought[branching_rule.crew_ix, i, :, :],
				)
			]
		end
		output[crew] = ixs
		ixs = copy(crew_avail_ixs[crew])
		if branching_rule.crew_ix == crew
			ixs = [
				i for i in ixs if satisfies_branching_rule(
					branching_rule,
					crew_routes.fires_fought[branching_rule.crew_ix, i, :, :],
				)
			]
		end
		output[crew] = ixs
	end

	return output
end

function apply_branching_rule(
	fire_avail_ixs::Vector{Vector{Int64}},
	fire_plans::FirePlanData,
	branching_rule::FireDemandBranchingRule,
)

	output = [Int64[] for fire in 1:size(fire_avail_ixs)[1]]
	for fire in 1:size(fire_avail_ixs)[1]
		ixs = copy(fire_avail_ixs[fire])
		if branching_rule.fire_ix == fire
			ixs = [
				i for i in ixs if satisfies_branching_rule(
					branching_rule,
					fire_plans.crews_present[branching_rule.fire_ix, i, :],
				)
			]
		end
		output[fire] = ixs
		ixs = copy(fire_avail_ixs[fire])
		if branching_rule.fire_ix == fire
			ixs = [
				i for i in ixs if satisfies_branching_rule(
					branching_rule,
					fire_plans.crews_present[branching_rule.fire_ix, i, :],
				)
			]
		end
		output[fire] = ixs
	end


	return output

end

function price_and_cut!!(
	rmp::RestrictedMasterProblem,
	cut_data::GUBCutData,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	gub_cut_limit_per_time::Int64,
	crew_rules::Vector{CrewSupplyBranchingRule},
	fire_rules::Vector{FireDemandBranchingRule},
	global_fire_allotment_rules::Vector{GlobalFireAllotmentBranchingRule},
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData;
	upper_bound::Float64 = 1e20,
	gub_cover_cuts = true,
	general_gub_cuts = true,
	single_fire_cuts = false
	)

	loop_ix = 1
	loop_max = 15

	most_recent_obj = 0

	while loop_ix < loop_max
		# run DCG, adding columns as needed
		double_column_generation!(
			rmp,
			crew_subproblems,
			fire_subproblems,
			crew_rules,
			fire_rules,
			global_fire_allotment_rules,
			crew_routes,
			fire_plans,
			cut_data,
			upper_bound = upper_bound,
		)
		if rmp.termination_status == MOI.OBJECTIVE_LIMIT
			@debug "no more cuts needed"
			break
		end

		current_obj = objective_value(rmp.model)
		@info "in cut generation" loop_ix current_obj
		if most_recent_obj / current_obj > 1 - 1e-8
			@info "halting cut generation early, too small improvement" loop_ix
			break
		end
		most_recent_obj = current_obj

		# add cuts
		num_cuts = find_and_incorporate_knapsack_gub_cuts!!(
			cut_data,
			gub_cut_limit_per_time,
			rmp,
			crew_routes,
			fire_plans,
			crew_subproblems,
			fire_subproblems,
			gub_cover_cuts = gub_cover_cuts,
			general_gub_cuts = general_gub_cuts,
			single_fire_cuts = single_fire_cuts
		)
		@debug "after cuts" num_cuts objective_value(rmp.model)
		if (num_cuts == 0)
			@debug "No cuts found, breaking" loop_ix
			break
		end

		if num_cuts == 0
			@debug "No cuts found, breaking" loop_ix
			break

		end

		loop_ix += 1

	end

	# if we got completely through all the loop of cut generation
	if loop_ix == loop_max
		@info "There may be more cuts; one last DCG to have provable lower bound" loop_ix
		double_column_generation!(
			rmp,
			crew_subproblems,
			fire_subproblems,
			crew_rules,
			fire_rules,
			global_fire_allotment_rules,
			crew_routes,
			fire_plans,
			cut_data,
			upper_bound = upper_bound,
		)
	end

	@debug "After price-and-cut" objective_value(rmp.model)
	@debug "Afer price-and-cut" length(eachindex(rmp.routes)) length(
		eachindex(rmp.plans),
	) length([i for i in eachindex(rmp.plans) if reduced_cost(rmp.plans[i]) < 1e-2]) length([
		i for i in eachindex(rmp.routes) if reduced_cost(rmp.routes[i]) < 1e-2
	]) dual.(rmp.supply_demand_linking) dual.(
		rmp.fire_allotment_branches
	) rmp.fire_allotment_branches

end

function find_integer_solution(
	solved_rmp::RestrictedMasterProblem,
	upper_bound::Float64,
	fire_column_limit::Int64,
	crew_column_limit::Int64,
	time_limit::Union{Int, Float64};
	warm_start_routes = nothing,
	warm_start_plans = nothing,
)
	# get reduced costs
	rc_plans = reduced_cost.(solved_rmp.plans)
	rc_routes = reduced_cost.(solved_rmp.routes)

	# get relaxed objective value
	lp_objective = objective_value(solved_rmp.model)

	if fire_column_limit > length(rc_plans)
		@debug "fewer fire plans than limit" length(rc_plans) fire_column_limit
		fire_column_limit = length(rc_plans)
	end

	if crew_column_limit > length(rc_routes)
		@debug "fewer crew routes than limit" length(rc_routes) crew_column_limit
		crew_column_limit = length(rc_routes)
	end


	## find variables we will be allowed to use

	# get keys in vector form
	plan_keys = [i for i in eachindex(rc_plans)]

	# get teduced costs in vector form
	plan_values = [rc_plans[ix] for ix in plan_keys]

	# get the indices of the first "fire_column_limit" indices of plan_values
	sorted_plan_ixs = Int.(plan_values .* 0)
	used_plan_ixs =
		[
			i for
			i in partialsortperm!(sorted_plan_ixs, plan_values, 1:fire_column_limit)
		]

	# get the max reduced cost selected
	max_plan_rc = plan_values[sorted_plan_ixs[fire_column_limit]]

	# if this reduced cost is too high for the upper bound
	if lp_objective + max_plan_rc > upper_bound + 1e-5

		# delete elements of used_plan_ixs
		to_delete = Int64[]
		for (i, ix) in enumerate(used_plan_ixs)
			if lp_objective + plan_values[ix] > upper_bound + 1e-5
				push!(to_delete, i)
			end
		end
		@debug "removing fire plans due to reduced-cost bound" length(to_delete)
		deleteat!(used_plan_ixs, to_delete)
	end

	# get the keys corresponding to these indices
	used_plan_keys = plan_keys[used_plan_ixs]
	unused_plan_keys = [i for i in plan_keys if i ∉ used_plan_keys]


	route_keys = [i for i in eachindex(rc_routes)]
	route_values = [rc_routes[ix] for ix in route_keys]
	sorted_route_ixs = Int.(route_values .* 0)
	used_route_ixs =
		[
			i for i in
			partialsortperm!(sorted_route_ixs, route_values, 1:crew_column_limit)
		]
	max_route_rc = route_values[sorted_route_ixs[crew_column_limit]]

	# if this reduced cost is too high for the upper bound
	if lp_objective + max_route_rc > upper_bound + 1e-5

		# delete elements of used_route_ixs
		to_delete = Int64[]
		for (i, ix) in enumerate(used_route_ixs)
			if lp_objective + route_values[ix] > upper_bound + 1e-5
				push!(to_delete, i)
			end
		end
		@debug "removing crew routes due to reduced-cost bound" length(to_delete)

		deleteat!(used_route_ixs, to_delete)
	end

	used_route_keys = route_keys[used_route_ixs]
	unused_route_keys = [i for i in route_keys if i ∉ used_route_keys]

	# fix unused variables to 0, set binary variables
	for i ∈ unused_route_keys
		if ~isnothing(warm_start_routes) && (i ∈ eachindex(warm_start_routes)) &&
		   (warm_start_routes[i] > 0.99)
			@warn "Route pruned by reduced cost bound but part of warm start solution, should not happen if warm start solution is best UB"
		end
		fix(solved_rmp.routes[i], 0, force = true)
	end
	for i ∈ unused_plan_keys
		if ~isnothing(warm_start_plans) && (i ∈ eachindex(warm_start_plans)) && (warm_start_plans[i] > 0.99)
			@warn "Plan pruned by reduced cost bound but part of warm start solution, should not happen if warm start solution is best UB"
		end
		fix(solved_rmp.plans[i], 0, force = true)
	end
	for i ∈ used_route_keys
		set_binary(solved_rmp.routes[i])
	end
	for i ∈ used_plan_keys
		set_binary(solved_rmp.plans[i])
	end

	if ~isnothing(warm_start_plans)
		@info "pushing solution to IP"
		plans_ws = 0
		for i ∈ eachindex(solved_rmp.plans)
			if (i ∈ eachindex(warm_start_plans)) && (warm_start_plans[i] > 0.99)
				plans_ws += 1
				set_start_value(solved_rmp.plans[i], 1)
			else
				set_start_value(solved_rmp.plans[i], 0)
			end
		end
		@debug "ws_plans" plans_ws
		routes_ws = 0
		for i ∈ eachindex(solved_rmp.routes)
			if (i ∈ eachindex(warm_start_routes)) && (warm_start_routes[i] > 0.99)
				routes_ws += 1
				set_start_value(solved_rmp.routes[i], 1)
			else
				set_start_value(solved_rmp.routes[i], 0)
			end
		end
		@debug "ws_routes" routes_ws

	end

	# set time TimeLimit and MIPFocus
	set_optimizer_attribute(solved_rmp.model, "TimeLimit", time_limit)
	set_optimizer_attribute(solved_rmp.model, "OutputFlag", 1)

	# solve
	optimize!(solved_rmp.model)

	# get values to return
	obj = Inf
	routes_used = nothing
	plans_used = nothing
	if result_count(solved_rmp.model) >= 1
		obj = objective_value(solved_rmp.model)
		routes_used = Int.(round.(value.(solved_rmp.routes)))
		plans_used = Int.(round.(value.(solved_rmp.plans)))
	end
	obj_bound = objective_bound(solved_rmp.model)


	return obj, obj_bound, routes_used, plans_used

end

function heuristic_upper_bound!!(
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	explored_bb_node::BranchAndBoundNode,
	num_iters::Int,
	soft_time_limit::Float64,
	kill_if_no_improvement_rounds::Int64,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	gub_cut_limit_per_time::Int64,
	gurobi_env;
	routes_best_sol = nothing,
	plans_best_sol = nothing
)
	start_time = time()
	@info "Finding heuristic upper bound" explored_bb_node.ix

	# gather global information
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	cut_data = deepcopy(explored_bb_node.cut_data)

	# keep all columns from explored node
	crew_ixs =
		[
			[i[1] for i in eachindex(explored_bb_node.master_problem.routes[j, :])]
			for j ∈ 1:num_crews
		]
	fire_ixs =
		[
			[i[1] for i in eachindex(explored_bb_node.master_problem.plans[j, :])]
			for j ∈ 1:num_fires
		]

	# add in columns from best solution
	if ~isnothing(routes_best_sol)
		for (crew, route) ∈ eachindex(routes_best_sol)
			if value(routes_best_sol[(crew, route)]) > 0.99
				if route ∉ crew_ixs[crew]
					push!(crew_ixs[crew], route) 
					@info "pushing route best sol"
				end
			end
		end

		for (fire, plan) ∈ eachindex(plans_best_sol)
			if value(plans_best_sol[(fire, plan)]) > 0.99
				if plan ∉ fire_ixs[fire]
					push!(fire_ixs[fire], plan) 
					@info "pushing plan best sol"
				end
			end
		end
	end

	# grab any fire and crew branching rules
	crew_rules = CrewSupplyBranchingRule[]
	fire_rules = FireDemandBranchingRule[]
	cur_node = explored_bb_node
	while ~isnothing(cur_node)

		crew_rules = vcat(cur_node.new_crew_branching_rules, crew_rules)
		fire_rules = vcat(cur_node.new_fire_branching_rules, fire_rules)
		cur_node = cur_node.parent
	end

	global_rules = GlobalFireAllotmentBranchingRule[]
	ub = Inf
	ub_rmp = nothing
	routes = routes_best_sol
	plans = plans_best_sol
	rounds_since_improvement = 0


	fractional_weighted_average_solution, _ =
		get_fire_and_crew_incumbent_weighted_average(explored_bb_node.master_problem,
			crew_routes,
			fire_plans,
		)

	current_allotment = Int.(ceil.(fractional_weighted_average_solution))


	for iter ∈ 1:num_iters

		cur_time = time() - start_time
		if cur_time > soft_time_limit
			@info "Breaking heuristic round" iter-1 cur_time soft_time_limit
			break
		end


		# get rid of past cuts
		# do we want cuts? lose guarantee of feasibility in next step because we may find a cut later 
		# cut_data = GUBCutData(num_crews, num_fires, num_time_periods)

		branching_rule =
			GlobalFireAllotmentBranchingRule(current_allotment,
				false,
				fire_plans,
				fire_subproblems,
			)
		global_rules = [branching_rule]

		for fire ∈ 1:num_fires
			to_delete = Int64[]
			for ix ∈ eachindex(fire_ixs[fire])
				plan_ix = fire_ixs[fire][ix]
				if (fire, plan_ix) ∈ keys(branching_rule.mp_lookup)
					push!(to_delete, ix)
				end
			end
			deleteat!(fire_ixs[fire], to_delete)
		end

		@debug "entering heuristic round" branching_rule.allotment_matrix
		rmp = define_restricted_master_problem(
			gurobi_env,
			crew_routes,
			crew_ixs,
			fire_plans,
			fire_ixs,
			cut_data,
			global_rules,
		)


		# TODO consider cut management
		t = @elapsed price_and_cut!!(
			rmp,
			cut_data,
			crew_subproblems,
			fire_subproblems,
			gub_cut_limit_per_time,
			crew_rules,
			fire_rules,
			global_rules,
			crew_routes,
			fire_plans,
			upper_bound = ub)
		@info "Price and cut time (heuristic)" t 

		# add in columns from best feasible solution so far
		if ~isnothing(routes)

			# this should usually be unnecessary since we don't often branch on crews
			for (crew, route) ∈ eachindex(routes)
				if (routes[(crew, route)] > 0.99) && ((crew, route) ∉ eachindex(rmp.routes))
					add_column_to_master_problem!!(rmp, cut_data, crew_routes, crew, route)
					@warn "Adding crew column to heuristic from best feasible solution, why was it dropped?"
				end
			end

			for (fire, plan) ∈ eachindex(plans)
				if (plans[(fire, plan)] > 0.99) && ((fire, plan) ∉ eachindex(rmp.plans))
					add_column_to_master_problem!!(rmp, cut_data, fire_plans, global_rules, fire, plan) 
				end
			end
		end

		# remove capacity constraint to allow best feasible solution so far
		set_normalized_rhs.(rmp.fire_allotment_branches, -10000)
		optimize!(rmp.model)

		crew_ixs = [[i[1] for i in eachindex(rmp.routes[j, :])] for j ∈ 1:num_crews]
		fire_ixs = [[i[1] for i in eachindex(rmp.plans[j, :])] for j ∈ 1:num_fires]


		lb_node =
			rmp.termination_status == MOI.LOCALLY_SOLVED ?
			objective_value(rmp.model) : -Inf

		if lb_node / ub > 1 - 1e-9
			@debug "pruning heuristic round by bound, bumping allotments"
			current_allotment = current_allotment .+ 1
			continue
		end

		# TODO add global branching rule as a valid inequality to help solver find cuts

		t = @elapsed obj, obj_bound, routes, plans =
			find_integer_solution(
				rmp,
				ub,
				120000,
				400000,
				40.0,
				warm_start_plans = plans,
				warm_start_routes = routes,
			)
		
		rounds_since_improvement += 1
		if obj < ub - 1e-9
			ub = obj
			ub_rmp = rmp
			rounds_since_improvement = 0
		elseif obj == Inf && ub < Inf
			@warn "Failure of warm start solution"
		end
		@info "found sol" t obj obj_bound
		if rounds_since_improvement >= kill_if_no_improvement_rounds
			@info "Too long since improvement in heuristic, killing early" rounds_since_improvement
			break
		end

		fire_allots, _ =
			get_fire_and_crew_incumbent_weighted_average(ub_rmp,
				crew_routes,
				fire_plans,
			)

		@debug "solution bounds" t lb_node ub obj obj_bound fire_allots

		current_allotment = current_allotment .+ 1
	end
	@info "Found heuristic upper bound" ub

	return ub, ub_rmp, routes, plans
end


function explore_node!!(
	branch_and_bound_node::BranchAndBoundNode,
	all_nodes::Vector{BranchAndBoundNode},
	current_global_upper_bound::Float64,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	gub_cut_limit_per_time::Int64,
	warm_start_strategy::Union{String, Nothing},
	gurobi_env;
	rel_tol = 1e-9,
	restore_integrality = false,
	gub_cover_cuts=true,
	general_gub_cuts=true,
	single_fire_cuts=false)

	@info "Exploring node" branch_and_bound_node.ix

	# gather global information
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)
	cut_data = branch_and_bound_node.cut_data
	## get the columns with which to initialize restricted master problem

	# if we are at the root node, there are no columns yet
	if isnothing(branch_and_bound_node.parent)
		crew_ixs = [Int[] for i ∈ 1:num_crews]
		fire_ixs = [Int[] for i ∈ 1:num_fires]


	else
		# if we are not at the root node, there are a lot of options here, but
		# for now take all the columns generated by the parent RMP that satisfy 
		# the new branching rule. Maybe could improve performance by culling
		# columns based on reduced costs or absence in basis, but first try at this
		# showed worse performance

		cull = false
		parent_rmp = branch_and_bound_node.parent.master_problem
		if ~cull
			crew_ixs =
				[[i[1] for i in eachindex(parent_rmp.routes[j, :])] for j ∈ 1:num_crews]
			fire_ixs =
				[[i[1] for i in eachindex(parent_rmp.plans[j, :])] for j ∈ 1:num_fires]
		else
			# TODO if we want to cull, need some better solution of timing of restore_integrality
			# in parent node
			@warn "If we restored integrality in parent node, only getting these columns"
			crew_ixs =
			[[i[1] for i in eachindex(parent_rmp.routes[j, :]) if value(parent_rmp.routes[j, i...]) > 1e-4] for j ∈ 1:num_crews]
			fire_ixs =
			[[i[1] for i in eachindex(parent_rmp.plans[j, :]) if value(parent_rmp.plans[j, i...])  > 1e-4] for j ∈ 1:num_fires]
		end
		@debug "num ix" length(fire_ixs[1])

		for rule in branch_and_bound_node.new_crew_branching_rules
			crew_ixs = apply_branching_rule(crew_ixs, crew_routes, rule)
		end
		for rule in branch_and_bound_node.new_fire_branching_rules
			fire_ixs = apply_branching_rule(fire_ixs, fire_plans, rule)
		end

		# remove plans excluded by left branching rules
		for rule in branch_and_bound_node.new_global_fire_allotment_branching_rules
			if ~rule.geq_flag
				for fire in 1:num_fires
					fire_ixs[fire] = filter!(
						x -> (fire, x) ∉ keys(rule.mp_lookup),
						fire_ixs[fire],
					)
				end
			end
		end
	end


	# get the branching rules
	crew_rules = CrewSupplyBranchingRule[]
	fire_rules = FireDemandBranchingRule[]
	global_rules = GlobalFireAllotmentBranchingRule[]

	cur_node = branch_and_bound_node
	while ~isnothing(cur_node)

		crew_rules = vcat(cur_node.new_crew_branching_rules, crew_rules)
		fire_rules = vcat(cur_node.new_fire_branching_rules, fire_rules)
		global_rules =
			vcat(cur_node.new_global_fire_allotment_branching_rules, global_rules)
		cur_node = cur_node.parent

	end

	# define the restricted master problem
	## TODO how do we handle existing cuts
	## currently carrying them all
	
	t = @elapsed rmp = define_restricted_master_problem(
		gurobi_env,
		crew_routes,
		crew_ixs,
		fire_plans,
		fire_ixs,
		cut_data,
		global_rules,
	)
	@info "Define rmp time (b-and-b)" t

	t = @elapsed price_and_cut!!(
		rmp,
		cut_data,
		crew_subproblems,
		fire_subproblems,
		gub_cut_limit_per_time,
		crew_rules,
		fire_rules,
		global_rules,
		crew_routes,
		fire_plans,
		upper_bound = current_global_upper_bound,
		gub_cover_cuts=gub_cover_cuts,
		general_gub_cuts=general_gub_cuts,
		single_fire_cuts=single_fire_cuts)

	@info "Price and cut time (b-and-b)" t
	@debug "after price and cut" objective_value(rmp.model) crew_routes.routes_per_crew fire_plans.plans_per_fire

	# extract some data
	binding_cuts = [
		i for i in eachindex(rmp.gub_cover_cuts) if
		dual(rmp.gub_cover_cuts[i]) > 1e-4
	]
	@debug "cut_duals" dual.(rmp.gub_cover_cuts)
	used_plans = [
		i for i in eachindex(rmp.plans) if
		value(rmp.plans[i]) > 1e-4
	]
	used_routes = [
		i for i in eachindex(rmp.routes) if
		value(rmp.routes[i]) > 1e-4
	]

	# update the rmp 
	branch_and_bound_node.master_problem = rmp

	all_fire_allots, all_crew_allots = extract_usages(crew_routes, fire_plans, rmp)
	@info "usages" all_fire_allots all_crew_allots

	# update the branch-and-bound node to be feasible or not
	if rmp.termination_status == MOI.INFEASIBLE
		branch_and_bound_node.feasible = false
		branch_and_bound_node.l_bound = Inf
		@debug "infeasible node" node_ix
	else
		if rmp.termination_status ∉ [MOI.LOCALLY_SOLVED, MOI.OBJECTIVE_LIMIT]
			@debug "Node not feasible or optimal, CG did not terminate?"
			error("Node not feasible or optimal, CG did not terminate?")
		else
			branch_and_bound_node.l_bound = objective_value(rmp.model)
		end
		branch_and_bound_node.feasible = true

		plan_values = value.(rmp.plans)
		route_values = value.(rmp.routes)

		if restore_integrality
			error("Doing this here loses rmp duals")
			# t = @elapsed obj, obj_bound, routes, plans =
			# 	find_integer_solution(
			# 		rmp,
			# 		current_global_upper_bound,
			# 		300,
			# 		1000,
			# 		0.5,
			# 	)
			# branch_and_bound_node.u_bound = obj
			# @debug "after restoring integrality" t obj obj_bound
		else
			tolerance = 1e-8
			integer =
				all((plan_values .< tolerance) .|| (plan_values .> 1 - tolerance)) &&
				all((route_values .< tolerance) .|| (route_values .> 1 - tolerance))
			if integer
				branch_and_bound_node.u_bound = branch_and_bound_node.l_bound
			end
		end
	end

	# if we cannot prune
	if (branch_and_bound_node.u_bound - branch_and_bound_node.l_bound > rel_tol) &
	   branch_and_bound_node.feasible &
	   (branch_and_bound_node.l_bound / current_global_upper_bound < 1 - rel_tol)

		# TODO think about branching rules

		@debug "duals" value.(rmp.supply_demand_linking)
		# decide the next branching rules
		branch_type, branch_ix, var_variance, var_mean =
			max_variance_natural_variable(
				crew_routes,
				fire_plans,
				route_values,
				plan_values,
			)

		# restrict to used cuts
		used_cuts = restrict_GUBCutData(cut_data, binding_cuts)
		@debug "cuts" used_cuts.cut_dict

		@assert var_variance > 0 "Cannot branch on variable with no variance, should already be integral"

		# create two new nodes with branching rules
		if branch_type == "fire"
			left_branching_rule = FireDemandBranchingRule(
				Tuple(branch_ix)...,
				Int(floor(var_mean)),
				"less_than_or_equal",
			)
			right_branching_rule = FireDemandBranchingRule(
				Tuple(branch_ix)...,
				Int(floor(var_mean)) + 1,
				"greater_than_or_equal",
			)
			left_child = BranchAndBoundNode(
				ix = size(all_nodes)[1] + 1,
				parent = branch_and_bound_node,
				new_fire_branching_rules = [left_branching_rule],
				cut_data = used_cuts,
			)
			right_child = BranchAndBoundNode(
				ix = size(all_nodes)[1] + 2,
				parent = branch_and_bound_node,
				new_fire_branching_rules = [right_branching_rule],
				cut_data = used_cuts,
			)
			push!(all_nodes, left_child)
			push!(all_nodes, right_child)
			branch_and_bound_node.children = [left_child, right_child]


		else
			left_branching_rule = CrewSupplyBranchingRule(
				Tuple(branch_ix)...,
				false,
			)
			right_branching_rule = CrewSupplyBranchingRule(
				Tuple(branch_ix)...,
				true,
			)
			left_child = BranchAndBoundNode(
				ix = size(all_nodes)[1] + 1,
				parent = branch_and_bound_node,
				new_crew_branching_rules = [left_branching_rule],
				cut_data = used_cuts,
			)
			right_child = BranchAndBoundNode(
				ix = size(all_nodes)[1] + 2,
				parent = branch_and_bound_node,
				new_crew_branching_rules = [right_branching_rule],
				cut_data = used_cuts,
			)
			push!(all_nodes, left_child)
			push!(all_nodes, right_child)
			branch_and_bound_node.children = [left_child, right_child]

		end
		@debug "branching rules" left_branching_rule right_branching_rule
	end

	return used_plans, used_routes, binding_cuts
end