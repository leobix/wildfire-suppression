include("CommonStructs.jl")
include("BranchingRules.jl")
include("GUBKnapsackCoverCuts.jl")

using Gurobi, Statistics


mutable struct BranchAndBoundNode

	const ix::Int64
	const parent::Union{Nothing, BranchAndBoundNode}
	const new_crew_branching_rules::Vector{CrewSupplyBranchingRule}
	const new_fire_branching_rules::Vector{FireDemandBranchingRule}
	const new_global_fire_allotment_branching_rules::Vector{
		GlobalFireAllotmentBranchingRule,
	}
	const cut_data::GUBCoverCutData
	children::Vector{BranchAndBoundNode}
	l_bound::Float64
	u_bound::Float64
	master_problem::Union{Nothing, RestrictedMasterProblem}
	feasible::Union{Nothing, Bool}
end

function BranchAndBoundNode(
	;
	ix::Int64,
	parent::Union{Nothing, BranchAndBoundNode},
	cut_data::GUBCoverCutData,
	new_crew_branching_rules::Vector{CrewSupplyBranchingRule} = CrewSupplyBranchingRule[],
	new_fire_branching_rules::Vector{FireDemandBranchingRule} = FireDemandBranchingRule[],
	new_global_fire_allotment_branching_rules::Vector{
		GlobalFireAllotmentBranchingRule,
	} = GlobalFireAllotmentBranchingRule[],
	children::Vector{BranchAndBoundNode} = BranchAndBoundNode[],
	l_bound::Float64 = -Inf,
	u_bound::Float64 = Inf,
	master_problem::Union{Nothing, RestrictedMasterProblem} = nothing,
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
	@debug "Means" crew_means # crew_means
	@debug "Variances" crew_variances # crew_variances
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
		if plan_values[ix] > 0.0001
			@debug "used plan" ix plan_values[ix] fire_plans.crews_present[
				fire,
				plan,
				:,
			]
		end
	end
	fire_variances = fire_sq_means - (fire_means .^ 2)
	@debug "Means" fire_means # crew_means
	@info "Variances" fire_variances # crew_variances
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
	cut_data::GUBCoverCutData,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	crew_rules::Vector{CrewSupplyBranchingRule},
	fire_rules::Vector{FireDemandBranchingRule},
	global_fire_allotment_rules::Vector{GlobalFireAllotmentBranchingRule},
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData;
	upper_bound::Float64=1e20)

	@debug "at price and cut start" length(rmp.plans[1, :]) length(rmp.plans[2, :]) length(
		rmp.plans[3, :],
	)

	loop_ix = 1
	loop_max = 8
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
			upper_bound = upper_bound
		)
		if rmp.termination_status == MOI.OBJECTIVE_LIMIT
			@info "no more cuts needed"
			break
		end
		@debug "objective after cg" objective_value(rmp.model)
		# add cuts
		num_cuts = find_and_incorporate_knapsack_gub_cuts!!(
			cut_data,
			rmp,
			crew_routes,
			fire_plans,
			crew_subproblems,
			fire_subproblems,
		)
		@info "after cuts" num_cuts objective_value(rmp.model)
		if num_cuts == 0
			@info "No cuts found, breaking" loop_ix
			break
		end

		loop_ix += 1

	end

	# if we got completely through all the loop of cut generation
	if loop_ix == loop_max
		@info "One last DCG to have provable lower bound"
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
			upper_bound = upper_bound
		)
	end

	@info "rmp stats end" objective_value(rmp.model) length(eachindex(rmp.routes)) length(
		eachindex(rmp.plans),
	) length([i for i in eachindex(rmp.plans) if reduced_cost(rmp.plans[i]) < 1e-2]) length([
		i for i in eachindex(rmp.routes) if reduced_cost(rmp.routes[i]) < 1e-2
	]) dual.(rmp.supply_demand_linking) dual.(rmp.gub_cover_cuts) cut_data.cut_dict dual.(
		rmp.fire_allotment_branches
	) rmp.fire_allotment_branches
	@debug "reduced_costs" reduced_cost.(rmp.plans) reduced_cost.(rmp.routes)

end

function find_integer_solution(
	solved_rmp::RestrictedMasterProblem,
	upper_bound::Float64,
	fire_column_limit::Int64,
	crew_column_limit::Int64,
	time_limit::Float64,
)
	# get reduced costs
	rc_plans = reduced_cost.(solved_rmp.plans)
	rc_routes = reduced_cost.(solved_rmp.routes)

	# get relaxed objective value
	lp_objective = objective_value(solved_rmp.model)

	if fire_column_limit > length(rc_plans)
		@info "fewer fire plans than limit" length(rc_plans) fire_column_limit
		fire_column_limit = length(rc_plans)
	end

	if crew_column_limit > length(rc_routes)
		@info "fewer crew routes than limit" length(rc_routes) crew_column_limit
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
	if lp_objective + max_plan_rc > upper_bound

		# delete elements of used_plan_ixs
		to_delete = Int64[]
		for (i, ix) in enumerate(used_plan_ixs)
			if lp_objective + plan_values[ix] > upper_bound
				push!(to_delete, i)
			end
		end
		@info "removing fire plans due to reduced-cost bound" length(to_delete)
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
	if lp_objective + max_plan_rc > upper_bound

		# delete elements of used_plan_ixs
		to_delete = Int64[]
		for (i, ix) in enumerate(used_route_ixs)
			if lp_objective + route_values[ix] > upper_bound
				push!(to_delete, i)
			end
		end
		@info "removing crew routes due to reduced-cost bound" length(to_delete)

		deleteat!(used_route_ixs, to_delete)
	end

	used_route_keys = route_keys[used_route_ixs]
	unused_route_keys = [i for i in route_keys if i ∉ used_route_keys]

	@debug "allowed columns for restore_integrality" length(used_plan_keys) length(
		used_route_keys,
	) maximum(plan_values[used_plan_ixs]) maximum(route_values[used_route_ixs])


	# fix unused variables to 0, set binary variables
	for i ∈ unused_route_keys
		fix(solved_rmp.routes[i], 0, force = true)
	end
	for i ∈ unused_plan_keys
		fix(solved_rmp.plans[i], 0, force = true)
	end
	for i ∈ used_route_keys
		set_binary(solved_rmp.routes[i])
	end
	for i ∈ used_plan_keys
		set_binary(solved_rmp.plans[i])
	end

	# set time TimeLimit
	set_optimizer_attribute(solved_rmp.model, "TimeLimit", time_limit)

	# solve
	optimize!(solved_rmp.model)

	# get values to return
	obj = Inf
	routes_used = nothing
	plans_used = nothing
	if is_solved_and_feasible(solved_rmp.model)
		obj = objective_value(solved_rmp.model)
		routes_used = value.(solved_rmp.routes)
		plans_used = value.(solved_rmp.plans)
	end
	obj_bound = objective_bound(solved_rmp.model)


	return obj, obj_bound, routes_used, plans_used

end

function heuristic_upper_bound!!(
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	num_iters::Int,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	gurobi_env,
)

	# gather global information
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	# initialize side data
	cut_data = GUBCoverCutData(num_crews, num_fires, num_time_periods)
	crew_ixs = [Int[] for i ∈ 1:num_crews]
	fire_ixs = [Int[] for i ∈ 1:num_fires]
	crew_rules = CrewSupplyBranchingRule[]
	fire_rules = FireDemandBranchingRule[]
	global_rules = GlobalFireAllotmentBranchingRule[]
	lb = -Inf
	ub = Inf

	## solve at root node so that we can get fractional weighted-average solution

	rmp = define_restricted_master_problem(
		gurobi_env,
		crew_routes,
		crew_ixs,
		fire_plans,
		fire_ixs,
		cut_data,
		global_rules,
	)

	# TODO smarter cuts that incorporate these upper bounds (crews_available)
	price_and_cut!!(
		rmp,
		cut_data,
		crew_subproblems,
		fire_subproblems,
		crew_rules,
		fire_rules,
		global_rules,
		crew_routes,
		fire_plans,
		upper_bound = ub)

	if rmp.termination_status == MOI.OPTIMAL
		lb = objective_value(rmp.model)
	end

	fractional_weighted_average_solution, _ =
		get_fire_and_crew_incumbent_weighted_average(rmp,
			crew_routes,
			fire_plans,
		)

	current_allotment = Int.(ceil.(fractional_weighted_average_solution))


	for i ∈ 1:num_iters

		branching_rule =
			GlobalFireAllotmentBranchingRule(current_allotment,
				false,
				fire_plans,
				fire_subproblems,
			)
		global_rules = [branching_rule]

		# TODO maybe we only keep active columns, though cut addition not so bad
		crew_ixs =
			[[i[1] for i in eachindex(rmp.routes[j, :])] for j ∈ 1:num_crews]
		fire_ixs =
			[[i[1] for i in eachindex(rmp.plans[j, :])] for j ∈ 1:num_fires]

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


		@info "entering heuristic round" branching_rule.allotment_matrix crew_ixs fire_ixs
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
		price_and_cut!!(
			rmp,
			cut_data,
			crew_subproblems,
			fire_subproblems,
			crew_rules,
			fire_rules,
			global_rules,
			crew_routes,
			fire_plans,
			upper_bound = ub)

		lb_node =
			rmp.termination_status == MOI.OPTIMAL ? objective_value(rmp.model) : Inf

		t = @elapsed obj, obj_bound, routes, plans =
			find_integer_solution(rmp, ub, 1200, 4000, 20.0)
		if obj < ub
			ub = obj
		end

		fire_allots, _ =
			get_fire_and_crew_incumbent_weighted_average(rmp,
				crew_routes,
				fire_plans,
			)

		@info "solution bounds" t lb_node ub obj obj_bound fire_allots

		binding_non_zero_allotments =
			(
				current_allotment .==
				Int.(round.(fire_allots))
			)

		if ~any(binding_non_zero_allotments)
			@info "No upper bounds are binding on the allotment, breaking"
			break
		end

		current_allotment = current_allotment .+ binding_non_zero_allotments
	end
end


function explore_node!!(
	branch_and_bound_node::BranchAndBoundNode,
	all_nodes::Vector{BranchAndBoundNode},
	current_global_upper_bound::Float64,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	warm_start_strategy::Union{String, Nothing},
	gurobi_env)

	# gather global information
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)
	cut_data = branch_and_bound_node.cut_data
	## get the columns with which to initialize restricted master problem

	# if we are at the root node, there are no columns yet
	if isnothing(branch_and_bound_node.parent)
		crew_ixs = [Int[] for i ∈ 1:num_crews]
		fire_ixs = [Int[] for i ∈ 1:num_fires]

		# if we are not at the root node, there are a lot of options here, but
		# for now take all the columns generated by the parent RMP that satisfy 
		# the new branching rule. Probably could improve performance by culling
		# columns based on reduced costs or absence in basis.
	else
		parent_rmp = branch_and_bound_node.parent.master_problem
		crew_ixs =
			[[i[1] for i in eachindex(parent_rmp.routes[j, :])] for j ∈ 1:num_crews]
		fire_ixs =
			[[i[1] for i in eachindex(parent_rmp.plans[j, :])] for j ∈ 1:num_fires]
		@debug "avaliable columns before cull" crew_ixs fire_ixs
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
					@info "fire_ixs" length(fire_ixs[fire])
					fire_ixs[fire] = filter!(
						x -> (fire, x) ∉ keys(rule.mp_lookup),
						fire_ixs[fire],
					)
					@info "fire_ixs" length(fire_ixs[fire])
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
	@debug "all branching rules found to pass to DCG" crew_rules fire_rules global_fire_allotment_branching_rules

	# define the restricted master problem
	## TODO how do we handle existing cuts
	## currently carrying them all
	rmp = define_restricted_master_problem(
		gurobi_env,
		crew_routes,
		crew_ixs,
		fire_plans,
		fire_ixs,
		cut_data,
		global_rules,
	)

	price_and_cut!!(
		rmp,
		cut_data,
		crew_subproblems,
		fire_subproblems,
		crew_rules,
		fire_rules,
		global_rules,
		crew_routes,
		fire_plans, 
		upper_bound = current_global_upper_bound)
	@info "after price and cut" objective_value(rmp.model) crew_routes.routes_per_crew fire_plans.plans_per_fire

	# extract some data
	binding_cuts = [
		i for i in eachindex(rmp.gub_cover_cuts) if
		dual(rmp.gub_cover_cuts[i]) > 1e-4
	]
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
	fire_alloc, crew_alloc = get_fire_and_crew_incumbent_weighted_average(
		rmp,
		crew_routes,
		fire_plans,
	)
	@info "usages" all_fire_allots all_crew_allots fire_alloc crew_alloc

	# update the branch-and-bound node to be feasible or not
	if rmp.termination_status == MOI.INFEASIBLE
		branch_and_bound_node.feasible = false
		branch_and_bound_node.l_bound = Inf
		@info "infeasible node" node_ix
	else
		if rmp.termination_status ∉ [MOI.LOCALLY_SOLVED, MOI.OBJECTIVE_LIMIT]
			@info "Node not feasible or optimal, CG did not terminate?"
			error("Node not feasible or optimal, CG did not terminate?")
		else
			branch_and_bound_node.l_bound = objective_value(rmp.model)
		end
		branch_and_bound_node.feasible = true

		plan_values = value.(rmp.plans)
		route_values = value.(rmp.routes)

		t = @elapsed obj, obj_bound, routes, plans =
			find_integer_solution(rmp, current_global_upper_bound, 300, 1000, 0.5)
		branch_and_bound_node.u_bound = obj
	end

	@info "after restoring integrality" t obj obj_bound

	@debug "used fire plans" [
		(ix, value(rmp.plans[ix]), fire_plans.crews_present[ix..., :]) for
		ix in eachindex(rmp.plans) if value(rmp.plans[ix]) > 0
	]
	@debug "used crew routes" [
		(ix, value(rmp.routes[ix])) for
		ix in eachindex(rmp.routes) if value(rmp.routes[ix]) > 0
	]


	# if we cannot prune
	if (branch_and_bound_node.u_bound - branch_and_bound_node.l_bound > 1e-5) &
	   branch_and_bound_node.feasible &
	   (branch_and_bound_node.l_bound < current_global_upper_bound)

		# TODO think about branching rules
		@info "fire_allots" all_fire_allots fire_alloc
		# cap = zeros(Int, num_fires, num_time_periods)
		# for g ∈ 1:num_fires
		# 	for t ∈ 1:num_time_periods
		# 		l = length(all_fire_allots[g, t])
		# 		if l > 0
		# 			if all_fire_allots[g, t][l][2] > 0
		# 				cap[g, t] = all_fire_allots[g, t][l][1] - 1
		# 			else
		# 				cap[g, t] = num_crews
		# 			end
		# 		else
		# 			cap[g, t] = num_crews
		# 		end
		# 	end
		# end

		# left_branching_rule =
		# 	GlobalFireAllotmentBranchingRule(Int.(ceil.(fire_alloc) .+ 2),
		# 		false,
		# 		fire_plans,
		# 		fire_subproblems,
		# 	)
		# # TODO examine side effects of shared reference
		# right_branching_rule = GlobalFireAllotmentBranchingRule(
		# 	left_branching_rule.allotment_matrix,
		# 	true,
		# 	deepcopy(left_branching_rule.fire_sp_arc_lookup),
		# 	deepcopy(left_branching_rule.mp_lookup),
		# )


		# left_child = BranchAndBoundNode(
		# 	ix = size(all_nodes)[1] + 1,
		# 	parent = branch_and_bound_node,
		# 	new_global_fire_allotment_branching_rules = [left_branching_rule],
		# )
		# right_child = BranchAndBoundNode(
		# 	ix = size(all_nodes)[1] + 2,
		# 	parent = branch_and_bound_node,
		# 	new_global_fire_allotment_branching_rules = [right_branching_rule],
		# )

		# decide the next branching rules
		branch_type, branch_ix, var_variance, var_mean =
			max_variance_natural_variable(
				crew_routes,
				fire_plans,
				route_values,
				plan_values,
			)

		# restrict to used cuts
		used_cuts = restrict_GUBCoverCutData(cut_data, binding_cuts)

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
				cut_data = used_cuts
			)
			right_child = BranchAndBoundNode(
				ix = size(all_nodes)[1] + 2,
				parent = branch_and_bound_node,
				new_fire_branching_rules = [right_branching_rule],
				cut_data = used_cuts
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
				cut_data = used_cuts
			)
			right_child = BranchAndBoundNode(
				ix = size(all_nodes)[1] + 2,
				parent = branch_and_bound_node,
				new_crew_branching_rules = [right_branching_rule],
				cut_data = used_cuts
			)
			push!(all_nodes, left_child)
			push!(all_nodes, right_child)
			branch_and_bound_node.children = [left_child, right_child]

		end
		@info "branching rules" left_branching_rule right_branching_rule
	end

	return used_plans, used_routes, binding_cuts
end




function test_BranchAndBoundNode()

	bb_node = BranchAndBoundNode(ix = 1, parent = nothing)
	println(bb_node)

end

# test_BranchAndBoundNode()
