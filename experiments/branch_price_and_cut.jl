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

        (g, c, t) = problem_size
        heuristic_enabled = heuristic_time_limit > 0.0
        s = string(c)
        if precompile
            s = "precompile"
        end
        file_name = s * "_cut_" * string(cut) * "_branch_strategy_" * string(b_rule) * "_heuristic_" * string(heuristic_enabled)
        local io = open(out_dir * "logs_" * file_name * ".txt", "w")
        if args["debug"] == true
            global_logger(ConsoleLogger(io, Logging.Debug, show_limited=false))
        else
            global_logger(ConsoleLogger(io, Logging.Info, show_limited=false))
        end
        explored_nodes, ubs, lbs, columns, times, time_1 =
            branch_and_price(
                g,
                c,
                t,
                algo_tracking=true,
                soft_heuristic_time_limit=heuristic_time_limit,
                total_time_limit=total_time_limit,
                branching_strategy=b_rule,
                bb_node_gub_cover_cuts=cut,
                bb_node_general_gub_cuts=cut,
                bb_node_decrease_gub_allots=cut,
                bb_node_single_fire_lift=cut,
                heuristic_gub_cover_cuts=cut,
                heuristic_general_gub_cuts=cut,
                heuristic_decrease_gub_allots=cut,
                heuristic_single_fire_lift=cut,
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
    end
end

out_dir = args["directory_output"]
cuts = [false]
branching_rules = ["linking_dual_max_variance"]
heuristic_time_limits = [180.0]

out_dir = args["directory_output"]
cuts = [true, false]
branching_rules = ["linking_dual_max_variance", "max_variance"]
heuristic_time_limits = [180.0, 0.0]

# precompile
sizes = [(3, 10, 14), (6, 20, 14), (9, 30, 14)]
run_experiment(out_dir, sizes, cuts, branching_rules, heuristic_time_limits, precompile=true, total_time_limit=5.0)

# experiment
sizes = [(9, 30, 14)]
run_experiment(out_dir, sizes, cuts, branching_rules, heuristic_time_limits, precompile=false, total_time_limit=1800.0)
