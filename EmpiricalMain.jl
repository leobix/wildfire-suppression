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

io = open("logs_precompile_5.txt", "w")
if args["debug"] == true
	global_logger(ConsoleLogger(io, Logging.Debug, show_limited = false))
else
	global_logger(ConsoleLogger(io, Logging.Info, show_limited = false))
end

# precompile
io = open("logs_60.txt", "w")
if args["debug"] == true
	global_logger(ConsoleLogger(io, Logging.Debug, show_limited = false))
else
	global_logger(ConsoleLogger(io, Logging.Info, show_limited = false))
end
explored_nodes, ubs, lbs, columns, heuristic_times, times, time_1, root_node_ip_sol, root_node_ip_sol_time, fire_arcs_used, crew_arcs_used = branch_and_price(20, 60, 14, from_empirical = true, travel_speed = 40.0 * 6.0)
@info "test" fire_arcs_used, crew_arcs_used
# io = open("logs_40.txt", "w")
# if args["debug"] == true
# 	global_logger(ConsoleLogger(io, Logging.Debug, show_limited = false))
# else
# 	global_logger(ConsoleLogger(io, Logging.Info, show_limited = false))
# end
# branch_and_price(8, 40, 14, from_empirical = true, travel_speed = 40.0 * 16.0)
# io = open("logs_60.txt", "w")
# if args["debug"] == true
# 	global_logger(ConsoleLogger(io, Logging.Debug, show_limited = false))
# else
# 	global_logger(ConsoleLogger(io, Logging.Info, show_limited = false))
# end
# branch_and_price(8, 60, 14, from_empirical = true, travel_speed = 40.0 * 16.0)
