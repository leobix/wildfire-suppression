include("BranchAndPrice.jl")

using JuMP, Gurobi, JSON, Profile, ArgParse, Logging, IterTools

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
    cut_data = GUBCutData(num_crews, num_fires, num_time_periods)

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
function branch_and_price(
    num_fires::Int,
    num_crews::Int,
    num_time_periods::Int;
    max_nodes=10000,
    algo_tracking=false,
    gub_cut_limit_per_time=10000,
    cut_loop_max=15,
    relative_improvement_cut_req=1e-8,
    soft_heuristic_time_limit=300.0,
    hard_heuristic_iteration_limit=10,
    heuristic_must_improve_rounds=2,
    heuristic_cadence=10,
    total_time_limit=1800.0,
    bb_node_gub_cover_cuts=true,
    bb_node_general_gub_cuts=true,
    bb_node_single_fire_cuts=false,
    bb_node_decrease_gub_allots=true,
    bb_node_single_fire_lift=true,
    heuristic_gub_cover_cuts=true,
    heuristic_general_gub_cuts=true,
    heuristic_single_fire_cuts=false,
    heuristic_decrease_gub_allots=true,
    heuristic_single_fire_lift=true,
    price_and_cut_file=nothing
)
    start_time = time()

    @info "Initializing data structures"
    # initialize input data
    crew_routes, fire_plans, crew_models, fire_models, cut_data =
        initialize_data_structures(num_fires, num_crews, num_time_periods)


    algo_tracking ?
    (@info "Checkpoint after initializing data structures" time() - start_time) :
    nothing

    explored_nodes = []
    ubs = []
    lbs = []
    columns = []
    times = []
    time_1 = time() - start_time

    # initialize nodes list with the root node
    nodes = BranchAndBoundNode[]
    first_node = BranchAndBoundNode(ix=1, parent=nothing, cut_data=cut_data)
    push!(nodes, first_node)

    # initialize global variables to track in branch-and-bound tree
    ub = Inf
    ub_ix::Int = -1
    routes_best_sol = nothing
    plans_best_sol = nothing

    ## breadth-first search for now, can get smarter/add options

    unexplored = [i for i in eachindex(nodes) if isnothing(nodes[i].master_problem)]
    node_lbs = []
    for i in unexplored
        if (~isnothing(nodes[i].parent))
            push!(node_lbs, nodes[i].parent.l_bound)
        else
            @warn "Node with no parent" i
            push!(node_lbs, -Inf)
        end
    end

    node_explored_count = 0
    # while there are more nodes to explore
    while length(unexplored) > 0

        node_explored_count += 1

        ### breadth first ###
        # node_ix = unexplored[1]

        ### raise lb ###
        node_ix = unexplored[findmin(node_lbs)[2]]

        # explore the next node
        explore_node!!(
            nodes[node_ix],
            nodes,
            ub,
            crew_routes,
            fire_plans,
            crew_models,
            fire_models,
            gub_cut_limit_per_time,
            nothing,
            GRB_ENV,
            restore_integrality=false,
            cut_loop_max=cut_loop_max,
            relative_improvement_cut_req=relative_improvement_cut_req,
            gub_cover_cuts=bb_node_gub_cover_cuts,
            general_gub_cuts=bb_node_general_gub_cuts,
            single_fire_cuts=bb_node_single_fire_cuts,
            decrease_gub_allots=bb_node_decrease_gub_allots,
            single_fire_lift=bb_node_single_fire_lift,
            log_cuts_file=price_and_cut_file,
        )

        if time() - start_time > total_time_limit
            @info "Full time limit reached"
            break
        end

        if (node_explored_count % heuristic_cadence == 1) &&
           (soft_heuristic_time_limit > 0.0)
            heuristic_ub, ub_rmp, routes_best_sol, plans_best_sol = heuristic_upper_bound!!(
                crew_routes,
                fire_plans,
                nodes[node_ix],
                hard_heuristic_iteration_limit,
                soft_heuristic_time_limit,
                heuristic_must_improve_rounds,
                crew_models,
                fire_models,
                gub_cut_limit_per_time,
                GRB_ENV,
                routes_best_sol=routes_best_sol,
                plans_best_sol=plans_best_sol,
                cut_loop_max=cut_loop_max,
                relative_improvement_cut_req=relative_improvement_cut_req,
                gub_cover_cuts=heuristic_gub_cover_cuts,
                general_gub_cuts=heuristic_general_gub_cuts,
                single_fire_cuts=heuristic_single_fire_cuts,
                decrease_gub_allots=heuristic_decrease_gub_allots,
                single_fire_lift=heuristic_single_fire_lift
            )

            if time() - start_time > total_time_limit
                @info "Full time limit reached"
                break
            end

            nodes[node_ix].heuristic_found_master_problem = ub_rmp

            if heuristic_ub < nodes[node_ix].u_bound
                nodes[node_ix].u_bound = heuristic_ub
            end
        end

        ## now that we have done the heuristic, can mess up rmp and lose duals

        # t = @elapsed obj, obj_bound, routes, plans =
        # 	find_integer_solution(
        # 		nodes[node_ix].master_problem,
        # 		ub,
        # 		1200,
        # 		4000,
        # 		10.0,)			
        # @info "restore_integrality" t obj obj_bound
        # if obj < nodes[node_ix].u_bound
        # 	nodes[node_ix].u_bound = obj
        # end


        # if this node has an integer solution, check if we have found 
        # a better solution than the incumbent
        # TODO keep track if it comes from heurisitc or no
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
        @info "number of nodes" node_explored_count length(nodes)
        @info "columns" sum(crew_routes.routes_per_crew) sum(
            fire_plans.plans_per_fire,
        )
        algo_tracking ?
        (@info "Time check" time() - start_time) :
        nothing

        if algo_tracking
            push!(explored_nodes, node_explored_count)
            push!(lbs, lb)
            push!(ubs, ub)
            push!(
                columns,
                sum(crew_routes.routes_per_crew) + sum(fire_plans.plans_per_fire),
            )
            push!(times, time() - start_time)
        end

        unexplored =
            [i for i in eachindex(nodes) if isnothing(nodes[i].master_problem)]
        node_lbs = []
        for i in unexplored
            if (~isnothing(nodes[i].parent))
                push!(node_lbs, nodes[i].parent.l_bound)
            else
                @warn "Node with no parent" i
                push!(node_lbs, -Inf)
            end
        end

        if node_explored_count >= max_nodes
            @info "Hit node limit, breaking" node_explored_count
            break
        end
    end

    return explored_nodes, ubs, lbs, columns, times, time_1

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

io = open("logs_precompile2.txt", "w")
if args["debug"] == true
    global_logger(ConsoleLogger(io, Logging.Debug, show_limited=false))
else
    global_logger(ConsoleLogger(io, Logging.Info, show_limited=false))
end

##########
# ROOT NODE CUT EXPERIMENT #
##########

params = Dict()
params["gub_only"] = Dict(:bb_node_gub_cover_cuts => true,
    :bb_node_general_gub_cuts => false,
    :bb_node_decrease_gub_allots => false,
    :bb_node_single_fire_lift => false
)
params["gub_plus_strengthen"] = Dict(:bb_node_gub_cover_cuts => true,
    :bb_node_general_gub_cuts => false,
    :bb_node_decrease_gub_allots => true,
    :bb_node_single_fire_lift => true
)
params["cglp_only"] = Dict(:bb_node_gub_cover_cuts => false,
    :bb_node_general_gub_cuts => true,
    :bb_node_decrease_gub_allots => false,
    :bb_node_single_fire_lift => false
)
params["everything"] = Dict(:bb_node_gub_cover_cuts => true,
    :bb_node_general_gub_cuts => true,
    :bb_node_decrease_gub_allots => true,
    :bb_node_single_fire_lift => true
)

# precompile
branch_and_price(3, 10, 14, algo_tracking=false)
for (key, param_set) in params
    branch_and_price(3, 10, 14, 
        algo_tracking=true, 
        soft_heuristic_time_limit=0.0, 
        gub_cut_limit_per_time=100000, 
        max_nodes=1, 
        cut_loop_max=100, 
        relative_improvement_cut_req=1e-25, 
        price_and_cut_file="data/experiment_outputs/20240416/" * key * "_cut_progress_precompile.json"; 
        param_set...)
end

# experiment
sizes = [(3, 10, 14), (6, 20, 14), (9, 30, 14), (12, 40, 14), (15, 50, 14)]
for (g, c, t) ∈ sizes
    for (key, param_set) in params
        branch_and_price(g, c, t, 
            algo_tracking=true, 
            soft_heuristic_time_limit=0.0, 
            gub_cut_limit_per_time=100000, 
            max_nodes=1, 
            cut_loop_max=100, 
            relative_improvement_cut_req=1e-25, 
            price_and_cut_file="data/experiment_outputs/20240416/" * key * "_cut_progress_" * string(g) * ".json"; 
            param_set...)
    end
end


# Profile.init()
# @profile branch_and_price(6, 20, 14, algo_tracking=true, soft_heuristic_time_limit=20.0, heuristic_cadence=5, total_time_limit=60.0)
# io2 = open("prof.txt", "w")
# Profile.print(io2, mincount=100)
# close(io2)
close(io)

error("done")
# sizes = [(3, 10, 14)]
sizes = [(3, 10, 14), (6, 20, 14), (9, 30, 14), (12, 40, 14), (15, 50, 14)]
cut_limits = [10000, 0]
heuristic_time_limits = [180.0, 0.0]
cglps = [true, false]
out_dir = "data/experiment_outputs/20240401/"

for (problem_size, cut_limit, heuristic_time_limit, cglp) ∈ product(sizes, cut_limits, heuristic_time_limits, cglps)

    (g, c, t) = problem_size
    heuristic_enabled = heuristic_time_limit > 0.0
    file_name = string(c) * "_" * string(cut_limit) * "_heuristic_" * string(heuristic_enabled) * "_cglp_" * string(cglp)
    local io = open(out_dir * "logs_" * file_name * ".txt", "w")
    if args["debug"] == true
        global_logger(ConsoleLogger(io, Logging.Debug, show_limited=false))
    else
        global_logger(ConsoleLogger(io, Logging.Info, show_limited=false))
    end
    try
        explored_nodes, ubs, lbs, columns, times, time_1 =
            branch_and_price(
                g,
                c,
                t,
                algo_tracking=true,
                gub_cut_limit_per_time=cut_limit,
                heuristic_cadence=5,
                soft_heuristic_time_limit=heuristic_time_limit,
                total_time_limit=1800.0,
                bb_node_general_gub_cuts=cglp
            )
        json_name = file_name * ".json"
        outputs = Dict(
            "explored_nodes" => explored_nodes,
            "upper_bounds" => ubs,
            "lower_bounds" => lbs,
            "times" => times,
            "num_columns" => columns,
            "times" => times,
            "init_time" => time_1,
        )

        open(out_dir * json_name, "w") do f
            JSON.print(f, outputs, 4)
        end
    catch e
        println(e)
    finally
        close(io)
    end
end


