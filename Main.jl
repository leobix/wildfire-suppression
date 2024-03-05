include("BranchAndPrice.jl")

using JuMP, Gurobi, Profile, ArgParse, Logging

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


	crew_routes = CrewRouteData(100000, num_fires, num_crews, num_time_periods)
	fire_plans = FirePlanData(100000, num_fires, num_time_periods)
	cut_data = GUBCoverCutData(num_crews, num_fires, num_time_periods)

	return crew_routes, fire_plans, crew_models, fire_models, cut_data
end

# TODO if this is a bottleneck, can cache lower bounds
# in a new field in branch and bound node
function find_lower_bound(node::BranchAndBoundNode)

	# child.l_bound = -Inf if unexplored
	child_lbs = [find_lower_bound(child) for child in node.children]

	if length(child_lbs) > 0
		return max(minimum(child_lbs), node.l_bound)
	else
		return node.l_bound
	end

end
function branch_and_price(num_fires::Int, num_crews::Int, num_time_periods::Int)

	# initialize input data
	crew_routes, fire_plans, crew_models, fire_models, cut_data =
		initialize_data_structures(num_fires, num_crews, num_time_periods)

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
			cut_data,
			nothing,
			GRB_ENV,
		)

		# if this node has an integer solution, check if we have found 
		# a better solution than the incumbent
        if nodes[node_ix].u_bound < ub
            ub = nodes[node_ix].u_bound
            ub_ix = node_ix
        end

		# calculate the best current lower bound by considering all nodes with
		# fully explored children 
		lb = find_lower_bound(nodes[1])

		# print progress

		@info "current bounds" node_ix lb ub

		# go to the next node
		node_ix += 1
		@info "number of nodes" node_ix length(nodes)
		@info "columns" crew_routes.routes_per_crew fire_plans.plans_per_fire

		if node_ix > 7
			println("halted early.")
            # for g in 1:num_fires
            #     num_plans = fire_plans.plans_per_fire[g]
            #     plans = eachrow(fire_plans.crews_present[g, 1:num_plans, :])
            #     plans = [i for i in plans if sum(i) > 0]
            #     @info plans
            #     @assert allunique(plans)
            # end

            # for c in 1:num_crews
            #     num_routes = crew_routes.routes_per_crew[c]
            #     routes = [crew_routes.fires_fought[c, i] for i in 1:num_routes]
            #     routes = [i for i in routes if sum(i) > 0]
            #     @info routes
            #     @assert allunique(routes)
            # end
            return
		end
	end

end

function get_command_line_args()
	arg_parse_settings = ArgParseSettings()
	@add_arg_table arg_parse_settings begin
		"--debug"
		help = "run in debug mode, exposing all logging that uses @debug macro"
		action = :store_true
	end
	return parse_args(arg_parse_settings)
end

args = get_command_line_args()
io = open("logs.txt", "w")
if args["debug"] == true
	global_logger(ConsoleLogger(io, Logging.Debug, show_limited = false))
else
	global_logger(ConsoleLogger(io, Logging.Info, show_limited = false))
end

branch_and_price(6, 20, 14)
# branch_and_price(9, 30, 14)
close(io)
