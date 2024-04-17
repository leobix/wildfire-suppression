include("BranchAndPrice.jl")

using JuMP, Gurobi, JSON, Profile, ArgParse, Logging, IterTools

const GRB_ENV = Gurobi.Env()


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


