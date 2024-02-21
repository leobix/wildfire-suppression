include("BranchAndPrice.jl")
include("DoubleColumnGeneration.jl")

using JuMP, Gurobi, Profile
const GRB_ENV = Gurobi.Env()

function initialize_data_structures(
	num_fires::Int64,
	num_crews::Int64,
	num_time_periods::Int64,
)
	crew_models = build_crew_models(
		"data/raw/big_fire",
		num_fires,
		num_crews,
		num_time_periods,
	)

	fire_models = build_fire_models(
		"data/raw/big_fire",
		num_fires,
		num_crews,
		num_time_periods,
	)


	crew_routes = CrewRouteData(10000, num_fires, num_crews, num_time_periods)
	fire_plans = FirePlanData(10000, num_fires, num_time_periods)

	return crew_routes, fire_plans, crew_models, fire_models
end


function branch_and_price(num_fires::Int, num_crews::Int, num_time_periods::Int)

    # initialize input data
	crew_routes, fire_plans, crew_models, fire_models =
		initialize_data_structures(3, 10, 14)

    # initialize nodes list with the root node
	nodes = BranchAndBoundNode[]
	first_node = BranchAndBoundNode(ix = 1, parent = nothing)
	push!(nodes, first_node)
	
    # initialize global variables to track in branch-and-bound tree
    node_ix = 1
	ub = Inf
    ub_ix::Int = -1

	## breadth-first search for now, can get smarter/add options
	
    # while there are more nodes to explore
    while node_ix <= size(nodes)[1]
		
        # explore the next node
        explore_node!!(
			nodes[node_ix],
			nodes,
			ub,
			crew_routes,
			fire_plans,
			crew_models,
			fire_models,
			nothing,
			GRB_ENV,
		)

        # if this node has an integer solution, check if we have found 
        # a better solution than the incumbent
        if nodes[node_ix].integer == true
            node_ub = objective_value(nodes[nodes_ix].master_problem.model)
            if node_ub < ub
                ub = node_ub
                ub_ix = ix
            end
        end

        # calculate the best current lower bound by considering all nodes with
        # fully explored children 
        lb = Inf

        # for each node
        for node in nodes

            # if it is an integer solution or has children
            if (node.integer == true) | (size(node.children)[1] > 0)

                # a vaild lower bound is the max of the bound found at this node and 
                # at its children (using the initialization of l_bound to Inf)
                node_lb = node.l_bound
                for child in node.children
                    node_lb = max(node_lb, child.l_bound)
                end

                # if this is less than the best lower found so far, update it
                if node_lb < lb
                    lb = node_lb
                end
            end
        end
            
        # print progress
        println("_____")
        println(node_ix)
        println(lb)
        println(ub)
        println("_____")

        # go to the next node
        node_ix += 1 
        println(length(nodes))
        println("_____")    
        if node_ix > 100
            println("yay.")
            return
        end
	end

end

branch_and_price(3, 10, 14)
