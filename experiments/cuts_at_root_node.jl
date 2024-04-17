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
for (key, param_set) in params
    branch_and_price(3, 10, 14, 
        algo_tracking=true, 
        soft_heuristic_time_limit=0.0, 
        gub_cut_limit_per_time=100000, 
        max_nodes=1, 
        cut_loop_max=100, 
        relative_improvement_cut_req=1e-25, 
        price_and_cut_file=args["directory_output"] * key * "_cut_progress_precompile.json"; 
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
            price_and_cut_file=args["directory_output"] * key * "_cut_progress_" * string(g) * ".json"; 
            param_set...)
    end
end