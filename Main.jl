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

io = open("logs_precompile_3.txt", "w")
if args["debug"] == true
	global_logger(ConsoleLogger(io, Logging.Debug, show_limited = false))
else
	global_logger(ConsoleLogger(io, Logging.Info, show_limited = false))
end

# precompile
branch_and_price(3, 10, 14, algo_tracking = false)

branch_and_price(
	6,
	20,
	14,
	algo_tracking = true,
	# price_and_cut_file = "timings_ref2.json",
	bb_node_gub_cover_cuts=true,
	total_time_limit = 600.0,
	soft_heuristic_time_limit = 0.0,
)



