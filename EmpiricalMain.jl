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
num_fires = 10
num_crews = 30
num_time_periods = 14
travel_speed = 40.0 * 6.0
GC.gc()
for t in 0:2

	crew_routes = nothing
	fire_plans = nothing
	crew_models = nothing
	fire_models = nothing
	cut_data = nothing

	result = branch_and_price(num_fires,
		num_crews,
		num_time_periods,
		from_empirical = true, 
		travel_speed = travel_speed,
		crew_routes = crew_routes,
		fire_plans = fire_plans,
		crew_models = crew_models,
		fire_models = fire_models,
		cut_data = cut_data
		)
		# Unpack as many variables as branch_and_price returns, e.g.:
	explored_nodes, ubs, lbs, columns, heuristic_times, times, time_1, root_node_ip_sol, root_node_ip_sol_time, fire_arcs_used, crew_arcs_used = result
	@info "test" fire_arcs_used, crew_arcs_used
	crew_routes, fire_plans, crew_models, fire_models, cut_data = initialize_data_structures(num_fires, num_crews, num_time_periods, 20, travel_speed, from_empirical = true)
	for g in 1:num_fires
		@info "before modify_in_arcs_and_out_arcs!" fire_models[g].state_in_arcs
		modify_in_arcs_and_out_arcs!(fire_models[g], t+1, fire_arcs_used[g])
		@info "after modify_in_arcs_and_out_arcs!" fire_models[g].state_in_arcs
	end
	for j in 1:num_crews
		modify_in_arcs_and_out_arcs!(crew_models[j], t+1, crew_arcs_used[j])
	end


end

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
