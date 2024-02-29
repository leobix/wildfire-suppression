include("CommonStructs.jl")
include("BranchingRules.jl")

using Gurobi, Statistics


mutable struct BranchAndBoundNode

	const ix::Int64
	const parent::Union{Nothing, BranchAndBoundNode}
	const new_crew_branching_rules::Vector{CrewSupplyBranchingRule}
	const new_fire_branching_rules::Vector{FireDemandBranchingRule}
	const new_global_fire_allotment_branching_rules::Vector{
		GlobalFireAllotmentBranchingRule,
	}
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
	fire_plans::FirePlanData)

	@debug "at price and cut start" length(rmp.plans[1, :]) length(rmp.plans[2, :]) length(rmp.plans[3, :])

	for loop_ix ∈ 1:8
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
		)
		@info "objective after cg" objective_value(rmp.model)
		# add cuts
		num_cuts = find_and_incorporate_knapsack_gub_cuts!!(
			cut_data,
			rmp,
			crew_routes,
			fire_plans,
			crew_subproblems,
			fire_subproblems,
		)
		@info "objective after cuts" objective_value(rmp.model)
		if num_cuts == 0
			@info "No cuts found, breaking" loop_ix
			break
		end

	end
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
	)

	@info "rmp stats end" dual.(rmp.gub_cover_cuts) dual.(rmp.fire_allotment_branches) rmp.fire_allotment_branches

end
function explore_node!!(
	branch_and_bound_node::BranchAndBoundNode,
	all_nodes::Vector{BranchAndBoundNode},
	current_global_upper_bound::Float64,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	cut_data::GUBCoverCutData,
	warm_start_strategy::Union{String, Nothing},
	gurobi_env)

	# gather global information
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

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
		fire_plans)
	@info "after price and cut" objective_value(rmp.model) crew_routes.routes_per_crew fire_plans.plans_per_fire

	# update the rmp 
	branch_and_bound_node.master_problem = rmp

	@info "used fire plans" [
		(ix, value(rmp.plans[ix]), fire_plans.crews_present[ix..., :], fire_plans.plan_costs[ix...]) for
		ix in eachindex(rmp.plans) if value(rmp.plans[ix]) > 0
	]
	@info "used crew routes" [
		(ix, value(rmp.routes[ix]), crew_routes.route_costs[ix...]) for
		ix in eachindex(rmp.routes) if value(rmp.routes[ix]) > 0
	]

	# update the branch-and-bound node to be feasible or not
	if rmp.termination_status == MOI.INFEASIBLE
		branch_and_bound_node.feasible = false
		branch_and_bound_node.l_bound = Inf
		@info "infeasible node" node_ix
	else
		branch_and_bound_node.feasible = true
		branch_and_bound_node.l_bound = objective_value(rmp.model)

		plan_values = value.(rmp.plans)
		route_values = value.(rmp.routes)
		# # update the branch-and-bound node to be integer or not
		# tolerance = 1e-4
		# integer =
		# 	all((plan_values .< tolerance) .| (plan_values .> 1 - tolerance)) &
		# 	all((route_values .< tolerance) .| (route_values .> 1 - tolerance))
		# branch_and_bound_node.integer = integer

		set_binary.(rmp.routes)
		set_binary.(rmp.plans)
		set_optimizer_attribute(rmp.model, "TimeLimit", 2)
		t = @elapsed optimize!(rmp.model)
		branch_and_bound_node.u_bound = objective_value(rmp.model)
	end

	@info "after restoring integrality" t objective_value(rmp.model) objective_bound(
		rmp.model,
	)
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
		fire_alloc, _ = get_fire_and_crew_incumbent_weighted_average(
			rmp,
			crew_routes,
			fire_plans,
		)
		@info fire_alloc
		left_branching_rule = GlobalFireAllotmentBranchingRule(Int.(ceil.(fire_alloc) .+ 1),
			false,
			fire_plans,
			fire_subproblems,
		)
		# TODO examine side effects of shared reference
		right_branching_rule = GlobalFireAllotmentBranchingRule(left_branching_rule.allotment_matrix, true, left_branching_rule.fire_sp_arc_lookup, left_branching_rule.mp_lookup)


		left_child = BranchAndBoundNode(
			ix = size(all_nodes)[1] + 1,
			parent = branch_and_bound_node,
			new_global_fire_allotment_branching_rules = [left_branching_rule]
		)
		right_child = BranchAndBoundNode(
			ix = size(all_nodes)[1] + 2,
			parent = branch_and_bound_node,
			new_global_fire_allotment_branching_rules = [right_branching_rule]
		)

		# # decide the next branching rules
		# branch_type, branch_ix, var_variance, var_mean =
		# 	max_variance_natural_variable(
		# 		crew_routes,
		# 		fire_plans,
		# 		route_values,
		# 		plan_values,
		# 	)

		# @assert var_variance > 0 "Cannot branch on variable with no variance, should already be integral"

		# # create two new nodes with branching rules
		# if branch_type == "fire"
		# 	left_branching_rule = FireDemandBranchingRule(
		# 		Tuple(branch_ix)...,
		# 		Int(floor(var_mean)),
		# 		"less_than_or_equal",
		# 	)
		# 	right_branching_rule = FireDemandBranchingRule(
		# 		Tuple(branch_ix)...,
		# 		Int(floor(var_mean)) + 1,
		# 		"greater_than_or_equal",
		# 	)
		# 	left_child = BranchAndBoundNode(
		# 		ix = size(all_nodes)[1] + 1,
		# 		parent = branch_and_bound_node,
		# 		new_fire_branching_rules = [left_branching_rule],
		# 	)
		# 	right_child = BranchAndBoundNode(
		# 		ix = size(all_nodes)[1] + 2,
		# 		parent = branch_and_bound_node,
		# 		new_fire_branching_rules = [right_branching_rule],
		# 	)
		# 	push!(all_nodes, left_child)
		# 	push!(all_nodes, right_child)
		# 	branch_and_bound_node.children = [left_child, right_child]


		# else
		# 	left_branching_rule = CrewSupplyBranchingRule(
		# 		Tuple(branch_ix)...,
		# 		false,
		# 	)
		# 	right_branching_rule = CrewSupplyBranchingRule(
		# 		Tuple(branch_ix)...,
		# 		true,
		# 	)
		# 	left_child = BranchAndBoundNode(
		# 		ix = size(all_nodes)[1] + 1,
		# 		parent = branch_and_bound_node,
		# 		new_crew_branching_rules = [left_branching_rule],
		# 	)
		# 	right_child = BranchAndBoundNode(
		# 		ix = size(all_nodes)[1] + 2,
		# 		parent = branch_and_bound_node,
		# 		new_crew_branching_rules = [right_branching_rule],
		# 	)
			push!(all_nodes, left_child)
			push!(all_nodes, right_child)
			branch_and_bound_node.children = [left_child, right_child]

		# end
		@info "branching rules" left_branching_rule.allotment_matrix left_branching_rule.geq_flag right_branching_rule.allotment_matrix right_branching_rule.geq_flag
	end
end




function test_BranchAndBoundNode()

	bb_node = BranchAndBoundNode(ix = 1, parent = nothing)
	println(bb_node)

end

# test_BranchAndBoundNode()
