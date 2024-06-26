include("../BranchAndPrice.jl")

using JuMP, Gurobi, JSON, Profile, ArgParse, Logging, IterTools

const GRB_ENV = Gurobi.Env()


function get_command_line_args()
	arg_parse_settings = ArgParseSettings()
	@add_arg_table arg_parse_settings begin
		"""--debug"""
		help = "run in debug mode, exposing all logging that uses @debug macro"
		action = :store_true
		"--directory_output", "-d"
		help = "directory to write outputs, must exist"
		arg_type = String
		default = "data/experiment_outputs/cglp_lazy_constraints/"
	end
	return parse_args(arg_parse_settings)
end


args = get_command_line_args()

io = open("logs_lazy_constraints.txt", "w")
if args["debug"] == true
	global_logger(ConsoleLogger(io, Logging.Debug, show_limited = false))
else
	global_logger(ConsoleLogger(io, Logging.Info, show_limited = false))
end

cglp_strategies = ["cutting_plane", "adaptive", "enumerate"]

# precompile
sizes = [(6, 20, 14)]
loop_max = 2

for (g, c, t) ∈ sizes
	for cglp_strat ∈ cglp_strategies

		branch_and_price(
			g,
			c,
			t,
			algo_tracking = true,
			bb_node_gub_cover_cuts = false,
			bb_node_general_gub_cuts = cglp_strat,
			max_nodes = 1,
			cut_loop_max = loop_max,
			total_time_limit = 1200.0,
			soft_heuristic_time_limit = 0.0,
			price_and_cut_file = args["directory_output"] * cglp_strat * "_" *
								 "precompile.json";)

	end
end


sizes = [(3, 10, 14), (6, 20, 14), (9, 30, 14), (12, 40, 14), (15, 50, 14)]
# sizes = [(3, 10, 14), (6, 20, 14)]

loop_max = 50

for (g, c, t) ∈ sizes

	for cglp_strat ∈ cglp_strategies

		branch_and_price(
			g,
			c,
			t,
			algo_tracking = true,
			bb_node_gub_cover_cuts = false,
			bb_node_general_gub_cuts = cglp_strat,
			max_nodes = 1,
			cut_loop_max = loop_max,
			total_time_limit = 1200.0,
			soft_heuristic_time_limit = 0.0,
			price_and_cut_file = args["directory_output"] * cglp_strat * "_" *
								 string(c) * ".json";)

	end
end
