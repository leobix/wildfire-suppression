include("../BranchAndPrice.jl")

using Gurobi, JSON, ArgParse, Logging, IterTools


function get_command_line_args()
    arg_parse_settings = ArgParseSettings()
    @add_arg_table arg_parse_settings begin
        "--debug"
        help = "run in debug mode, exposing all logging that uses @debug macro"
        action = :store_true
        "--directory_output", "-d"
        help = "directory to write outputs, must exist"
        arg_type = String
        default = "data/experiment_outputs/cuts_at_root_node/"
    end
    return parse_args(arg_parse_settings)
end

args = get_command_line_args()


const GRB_ENV = Gurobi.Env()

io = open("logs_cuts_at_root_node.txt", "w")
if args["debug"] == true
    global_logger(ConsoleLogger(io, Logging.Debug, show_limited=false))
else
    global_logger(ConsoleLogger(io, Logging.Info, show_limited=false))
end


params = Dict()
params["gub_only"] = Dict(:bb_node_gub_cover_cuts => true,
    :bb_node_general_gub_cuts => "none",
    :bb_node_decrease_gub_allots => false,
    :bb_node_single_fire_lift => false,
    :cut_search_enumeration_limit => 100000
)
params["gub_plus_strengthen"] = Dict(:bb_node_gub_cover_cuts => true,
    :bb_node_general_gub_cuts => "none",
    :bb_node_decrease_gub_allots => true,
    :bb_node_single_fire_lift => true,
    :cut_search_enumeration_limit => 100000
)
params["cglp_cutting_plane"] = Dict(:bb_node_gub_cover_cuts => false,
    :bb_node_general_gub_cuts => "cutting_plane",
    :bb_node_decrease_gub_allots => false,
    :bb_node_single_fire_lift => false,
    :cut_search_enumeration_limit => 100000
)

params["cglp_enumerate"] = Dict(:bb_node_gub_cover_cuts => false,
    :bb_node_general_gub_cuts => "enumerate",
    :bb_node_decrease_gub_allots => false,
    :bb_node_single_fire_lift => false,
    :cut_search_enumeration_limit => 100000
)

params["cglp_adaptive"] = Dict(:bb_node_gub_cover_cuts => false,
    :bb_node_general_gub_cuts => "adaptive",
    :bb_node_decrease_gub_allots => false,
    :bb_node_single_fire_lift => false,
    :cut_search_enumeration_limit => 10000
)
params["everything"] = Dict(:bb_node_gub_cover_cuts => true,
    :bb_node_general_gub_cuts => "adaptive",
    :bb_node_decrease_gub_allots => true,
    :bb_node_single_fire_lift => true,
    :cut_search_enumeration_limit => 10000
)

# precompile
for (key, param_set) in params
    _, _, _, _, _, _, _, root_sol, root_sol_time = branch_and_price(3, 10, 14, 
        line_per_crew = 20,
        algo_tracking=true, 
        soft_heuristic_time_limit=0.0, 
        price_and_cut_soft_time_limit=1200.0,
        max_nodes=1, 
        cut_loop_max=2, 
        relative_improvement_cut_req=1e-25, 
        root_node_ip = (key == "everything"),
        price_and_cut_file=args["directory_output"] * key * "_cut_progress_precompile.json"; 
        param_set...)
end

# experiment
sizes = [(3, 10, 14, 20), (6, 20, 14, 20), (9, 30, 14, 20), (12, 40, 14, 20), (15, 50, 14, 20), (18, 60, 14, 20), (21, 70, 14, 20)]
output = []
for (g, c, t, l) ∈ sizes
    for (key, param_set) in params
        _, _, _, _, _, _, _, root_sol, root_sol_time = branch_and_price(g, c, t, 
            line_per_crew=l,
            algo_tracking=true, 
            soft_heuristic_time_limit=0.0, 
            price_and_cut_soft_time_limit=1200.0,
            max_nodes=1, 
            cut_loop_max=50, 
            relative_improvement_cut_req=1e-25, 
            root_node_ip = (key == "everything"),
            price_and_cut_file=args["directory_output"] * key * "_cut_progress_" * string(g) * ".json"; 
            param_set...)

        if key == "everything"
            push!(output, (g, c, t, l, root_sol, root_sol_time))
        end
    end
end
open(args["directory_output"] * "cgip.json", "w") do f
    JSON.print(f, output, 4)
end