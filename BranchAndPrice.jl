include("CommonStructs.jl")

using Gurobi
# const GRB_ENV = Gurobi.Env()


## data processing

function preprocess_fire_data()

end

function preprocess_crew_data()

end


### cg

function intialize_column_generation()

end




end

function solve_crew_subproblem(
	ts_network::TimeSpaceNetwork,
	linking_constraint_duals,
	branching_rules,
)

end

function solve_fire_subproblem(
	ts_network::TimeSpaceNetwork,
	linking_constraint_duals,
	branching_rules,
)

end

function choose_natural_variable_for_branching(
	branch_and_bound_node::BranchAndBoundNode,
)

end

function update_master_problem()

end




function explore_node(branch_and_bound_node::BranchAndBoundNode,
	all_nodes::Vector{BranchAndBoundNode},
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
	warm_start_strategy::Union{String, Nothing},
	gurobi_env)

	# get active branching rules by following tree upward

	# get dual values at solution of parent node for dual warm start

	# find available crew and fire plans based on branching rule to warm-start DCG

	# define master problem
	rmp = define_restricted_master_problem(
		gurobi_env,
		crew_routes,
		crew_plan_ixs,
		fire_plans,
		fire_plan_ixs,
		dual_warm_start,
	)

	# run DCG, adding columns as needed
	double_column_generation!(rmp, crew_subproblems, fire_subproblems)

	# update the branch-and-bound node to be feasible or not

	# decide the branching rule


end

function branch_and_price()


end




function test_BranchAndBoundNode()

	bb_node = BranchAndBoundNode(ix = 1, parent_ix = -1)
	println(bb_node)

end

test_BranchAndBoundNode()
