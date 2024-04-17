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
        default = "data/experiment_outputs/branch_price_and_cut/"
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

cut_limits = [100000, 0]
heuristic_time_limits = [180.0, 0.0]
cglps = [true, false]

# precompile
sizes = [(3, 10, 14)]
out_dir = args["directory_output"]

for (problem_size, cut_limit, heuristic_time_limit, cglp) ∈ product(sizes, cut_limits, heuristic_time_limits, cglps)

    (g, c, t) = problem_size
    heuristic_enabled = heuristic_time_limit > 0.0
    file_name = "precompile" * "_" * string(cut_limit) * "_heuristic_" * string(heuristic_enabled) * "_cglp_" * string(cglp)
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
                total_time_limit=15.0,
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
        @error "Failed branch and price" e
    finally
        close(io)
    end
end

# experiment
sizes = [(3, 10, 14), (6, 20, 14), (9, 30, 14), (12, 40, 14), (15, 50, 14)]
sizes = [(3, 10, 14), (6, 20, 14)]
out_dir = args["directory_output"]

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

