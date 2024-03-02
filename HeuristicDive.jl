include("BranchAndPrice.jl")

using JuMP, Gurobi, Profile, ArgParse, Logging

const GRB_ENV = Gurobi.Env()

function find_heuristic_upper_bound(
	num_fires::Int64,
	num_crews::Int64,
	num_time_periods::Int64,
)
	crew_models = build_crew_models(
		"data/raw/big_fire",
		num_fires,
		num_crews,
		num_time_periods,
	)

	fire_models = build_fire_models(
		"data/raw/big_fire",
		num_fires,
		num_crews,
		num_time_periods,
	)


	crew_routes = CrewRouteData(100000, num_fires, num_crews, num_time_periods)
	fire_plans = FirePlanData(100000, num_fires, num_time_periods)

    heuristic_upper_bound!!(
        crew_routes,
        fire_plans,
        2,
        crew_models,
        fire_models,
        GRB_ENV
    )
end

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
io = open("logs.txt", "w")
io2 = open("prof.txt", "w")
if args["debug"] == true
	global_logger(ConsoleLogger(io, Logging.Debug, show_limited = false))
else
	global_logger(ConsoleLogger(io, Logging.Info, show_limited = false))
end

find_heuristic_upper_bound(3, 10, 14)

Profile.init()
@profile find_heuristic_upper_bound(6, 20, 14)
Profile.print(io2)

close(io)
close(io2)