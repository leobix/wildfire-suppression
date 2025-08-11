include("../BranchAndPrice.jl")

using JSON, ArgParse

function get_command_line_args()
	arg_parse_settings = ArgParseSettings()
	@add_arg_table arg_parse_settings begin
		"""--debug"""
		help = "run in debug mode, exposing all logging that uses @debug macro"
		action = :store_true
		"--directory_output", "-d"
		help = "directory to write outputs, must exist"
		arg_type = String
		default = "data/experiment_outputs/"
	end
	return parse_args(arg_parse_settings)
end

sizes = [
	(3, 10, 14, 20, 640.0),
	(6, 20, 14, 20, 240.0),
    (6, 20, 14, 20, 640.0),
    (6, 20, 14, 20, Inf),
	(9, 30, 14, 20, 240.0),
    (9, 30, 14, 20, 640.0),
    (9, 30, 14, 20, Inf),
	(12, 40, 14, 20, 640.0),
	(15, 50, 14, 20, 640.0),
	(18, 60, 14, 20, 640.0),
	(21, 70, 14, 20, 640.0),
]

l = []
for (num_fires, num_crews, num_time_periods, line_per_crew, travel_speed) in sizes
        _, _, crew_models, fire_models, _, _ =
                initialize_data_structures(
                        num_fires,
                        num_crews,
                        num_time_periods,
                        line_per_crew,
                        travel_speed,
                )
    GC.gc()
	d = Dict{Any, Any}("crews" => num_crews, "fires" => num_fires)
	num_crew_arcs = [size(crew_models[i].wide_arcs)[2] for i in 1:num_crews]
	d["crew_arcs"] = sum(num_crew_arcs)
	num_fire_arcs = [size(fire_models[i].wide_arcs)[2] for i in 1:num_fires]
	d["fire_arcs"] = sum(num_fire_arcs)
    d["fire_states"] = sum([size(fire_models[i].state_in_arcs)[1] for i in 1:num_fires])
    d["travel_speed"] = travel_speed
	all_arc_times = []
	for i in 1:num_crews
		arcs = crew_models[i].long_arcs
		travelling_arcs = findall(
			(arcs[:, CM.LOC_FROM] .!= arcs[:, CM.LOC_TO]) .||
			(arcs[:, CM.FROM_TYPE] .!= arcs[:, CM.TO_TYPE]),
		)
		arc_times =
			arcs[travelling_arcs, CM.TIME_TO] .- arcs[travelling_arcs, CM.TIME_FROM]
		all_arc_times = vcat(all_arc_times, arc_times)
	end
	arc_times_hist = Dict(k => count(==(k), all_arc_times) for k in unique(all_arc_times))
    d["travel_times"] = arc_times_hist
    push!(l, d)
    println(d)
	
end

args = get_command_line_args()
open(args["directory_output"] * "problem_sizes.json", "w") do f
	JSON.print(f, l, 4)
end
