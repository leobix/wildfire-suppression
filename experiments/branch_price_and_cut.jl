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

io = open("logs_branch_price_and_cut.txt", "w")
if args["debug"] == true
    global_logger(ConsoleLogger(io, Logging.Debug, show_limited=false))
else
    global_logger(ConsoleLogger(io, Logging.Info, show_limited=false))
end

function run_experiment(out_dir, sizes, cuts, branching_rules, heuristic_time_limits; total_time_limit, precompile)

    for (cut, b_rule, heuristic_time_limit, problem_size) ∈ product(cuts, branching_rules, heuristic_time_limits, sizes)

        (g, c, t, l, crew_speed) = problem_size
        heuristic_enabled = heuristic_time_limit > 0.0
        s = string(c) * "+" * string(l) * "+" * string(crew_speed)
        if precompile
            s = "precompile"
        end
        file_name = s * "_cut_" * string(cut) * "_branch_strategy_" * string(b_rule) * "_heuristic_" * string(heuristic_enabled)
        logfile = out_dir * "logs_" * file_name * ".txt"

        # check if we tried this file, skip if logs exist
        if (~precompile) && (isfile(logfile))
            println("Skipping ", file_name)
            continue
        end

        local io = open(logfile, "w")
        if args["debug"] == true
            global_logger(ConsoleLogger(io, Logging.Debug, show_limited=false))
        else
            global_logger(ConsoleLogger(io, Logging.Info, show_limited=false))
        end
            
        explored_nodes, ubs, lbs, columns, heuristic_times, times, time_1 =
            branch_and_price(
                g,
                c,
                t,
                line_per_crew=l,
                travel_speed = crew_speed,
                algo_tracking=true,
                soft_heuristic_time_limit=heuristic_time_limit,
                total_time_limit=total_time_limit,
                branching_strategy=b_rule,
                bb_node_general_gub_cuts=cut ? "adaptive" : "none",
                heuristic_general_gub_cuts=cut ? "adaptive" : "none"
            )
        json_name = file_name * ".json"
        outputs = Dict(
            "explored_nodes" => explored_nodes,
            "upper_bounds" => ubs,
            "lower_bounds" => lbs,
            "times" => times,
            "heuristic_times" => heuristic_times,
            "num_columns" => columns,
            "init_time" => time_1,
        )

        open(out_dir * json_name, "w") do f
            JSON.print(f, outputs, 4)
        end
        # return
    end
end

out_dir = args["directory_output"]
cuts = [true, false]
branching_rules = ["most_fractional", "max_variance", "linking_dual_max_variance"]
heuristic_time_limits = [0.0, 60.0]

# precompile
sizes = [(3, 10, 14, 20, 640.0)]
run_experiment(out_dir, sizes, cuts, branching_rules, heuristic_time_limits, precompile=true, total_time_limit=5.0)

# value of cuts, branching, heuristic experiment
cuts = [true, false]
branching_rules = ["most_fractional", "max_variance", "linking_dual_max_variance"]
heuristic_time_limits = [0.0, 60.0]
sizes = [(3, 10, 14, 20, 640.0), (6, 20, 14, 20, 640.0), (9, 30, 14, 20, 640.0)]
run_experiment(out_dir, sizes, cuts, branching_rules, heuristic_time_limits, precompile=false, total_time_limit=1200.0)


# sensitivity experiment
branching_rules = ["linking_dual_max_variance"]
heuristic_time_limits = [60.0]
cuts = [true]
sizes_1 = [(6, 20, 14, i, j) for i ∈ [16, 18, 20, 22, 24] for j ∈ [640.0, 240.0, Inf]] 
sizes_2 = [(9, 30, 14, i, j) for i ∈ [16, 18, 20, 22, 24] for j ∈ [640.0, 240.0, Inf]] 
sizes = vcat(sizes_1, sizes_2)
run_experiment(out_dir, sizes, cuts, branching_rules, heuristic_time_limits, precompile=false, total_time_limit=1200.0)

# scalability experiment
branching_rules = ["linking_dual_max_variance"]
heuristic_time_limits = [60.0]
cuts = [true]
sizes = [(3, 10, 14, 20, 640.0), (6, 20, 14, 20, 640.0), (9, 30, 14, 20, 640.0), (12, 40, 14, 20, 640.0), (15, 50, 14, 20, 640.0), (18, 60, 14, 20, 640.0), (21, 70, 14, 20, 640.0)]
run_experiment(out_dir, sizes, cuts, branching_rules, heuristic_time_limits, precompile=false, total_time_limit=1200.0)
