include("BranchAndPrice.jl")

using JuMP, Gurobi, JSON, Profile, ArgParse, Logging, IterTools, CSV, DataFrames
import Logging: min_enabled_level, shouldlog, handle_message

struct DualLogger <: AbstractLogger
        loggers::NTuple{2,AbstractLogger}
end

min_enabled_level(l::DualLogger) = min(min_enabled_level(l.loggers[1]), min_enabled_level(l.loggers[2]))
shouldlog(l::DualLogger, level, _module, group, id) =
        shouldlog(l.loggers[1], level, _module, group, id) ||
        shouldlog(l.loggers[2], level, _module, group, id)
handle_message(l::DualLogger, level, message, _module, group, id, file, line; kwargs...) = begin
        handle_message(l.loggers[1], level, message, _module, group, id, file, line; kwargs...)
        handle_message(l.loggers[2], level, message, _module, group, id, file, line; kwargs...)
end

const GRB_ENV = Gurobi.Env()


const GACC_ABBR = Dict(
        "AK" => "Alaska",
        "EA" => "Eastern",
        "GB" => "Great Basin",
        "NC" => "Northern California",
        "CA-N" => "Northern California",
        "NR" => "Northern Rockies",
        "NW" => "Northwest",
        "RM" => "Rocky Mountain",
        "SA" => "Southern",
        "SC" => "Southern California",
        "CA-S" => "Southern California",
        "SW" => "Southwest",
)

const ALL_GACCS = unique(collect(values(GACC_ABBR)))

function parse_gaccs(str::String)
        s = uppercase(strip(str))
        if s == "ALL"
                return ALL_GACCS
        elseif s == "ALL_NO_AK"
                return filter(!=("Alaska"), ALL_GACCS)
        else
                abbrs = split(s, ",")
                return [haskey(GACC_ABBR, a) ? GACC_ABBR[a] : error("Unknown GACC abbreviation: $a") for a in abbrs]
        end
end

function parse_fires_by_gacc(str::String)
        s = strip(str)
        if isempty(s)
                return Dict{String,Vector{Int64}}()
        end
        result = Dict{String,Vector{Int64}}()
        groups = split(s, ';')
        for g in groups
                isempty(strip(g)) && continue
                parts = split(g, ':')
                length(parts) == 2 || error("Invalid GACC/fire mapping: $g")
                gacc_abbr = uppercase(strip(parts[1]))
                gacc = haskey(GACC_ABBR, gacc_abbr) ? GACC_ABBR[gacc_abbr] : error("Unknown GACC abbreviation: $gacc_abbr")
                fire_ids = [parse(Int, f) for f in split(strip(parts[2]), ',') if !isempty(strip(f))]
                result[gacc] = fire_ids
        end
        return result
end

function count_selected_fires(fire_gaccs::Vector{String}, fires_by_gacc::Dict{String,Vector{Int64}})
        fire_folder = "data/empirical_fire_models/raw/arc_arrays"
        selected_fires = CSV.read(joinpath(fire_folder, "selected_fires.csv"), DataFrame)
        if !isempty(fires_by_gacc)
                mask = falses(nrow(selected_fires))
                for (gacc, fires) in fires_by_gacc
                        mask .|= (selected_fires[:, "GACC"] .== gacc) .& in.(selected_fires[:, "FIRE_EVENT_ID"], Ref(fires))
                end
                selected_fires = selected_fires[mask, :]
        else
                selected_fires = selected_fires[in.(selected_fires[:, "GACC"], Ref(fire_gaccs)), :]
        end
        return length(unique(selected_fires[:, "FIRE_EVENT_ID"]))
end

function get_command_line_args()
        arg_parse_settings = ArgParseSettings()
        @add_arg_table arg_parse_settings begin
                "--debug"
                help = "run in debug mode, exposing all logging that uses @debug macro"
                action = :store_true
                "--crew-gaccs"
                help = "Allowed crew GACCs: 'all', 'all_no_ak', or comma-separated list of abbreviations (e.g. GB,NW)"
                default = "GB"
                "--fire-gaccs"
                help = "Allowed fire GACCs: 'all', 'all_no_ak', or comma-separated list of abbreviations (e.g. GB,NW). Defaults to crew GACCs if omitted"
                default = ""
                "--firefighters-per-crew"
                help = "Number of firefighters per crew"
                arg_type = Int
                default = 70
                "--personnel-per-crew"
                help = "Personnel count representing one crew when initializing from empirical data"
                arg_type = Int
                default = 20
                "--fires"
                help = "Mapping of GACC abbreviations to comma-separated FIRE_EVENT_IDs, separated by semicolons (e.g. GB:1,2;NW:3)"
                default = ""
                "--time-limit"
                help = "Time limit in seconds for the branch-and-price algorithm"
                arg_type = Float64
                default = 1800.0
                "--output-folder"
                help = "Directory to store output files"
                default = "data/output"
        end
        return parse_args(arg_parse_settings)
end


args = get_command_line_args()
crew_gaccs = parse_gaccs(args["crew-gaccs"])
fire_gaccs = isempty(args["fire-gaccs"]) ? crew_gaccs : parse_gaccs(args["fire-gaccs"])
firefighters_per_crew = args["firefighters-per-crew"]
personnel_per_crew = args["personnel-per-crew"]
fires_by_gacc = parse_fires_by_gacc(args["fires"])
time_limit = args["time-limit"]
output_folder = args["output-folder"]
mkpath(output_folder)

# send logs to both console and file so users can see initialization details
log_file = open("logs_$(Int(time_limit)).txt", "w")
if args["debug"] == true
        console_logger = ConsoleLogger(stdout, Logging.Debug, show_limited = false)
        file_logger = ConsoleLogger(log_file, Logging.Debug, show_limited = false)
else
        console_logger = ConsoleLogger(stdout, Logging.Info, show_limited = false)
        file_logger = ConsoleLogger(log_file, Logging.Info, show_limited = false)
end

global_logger(DualLogger((console_logger, file_logger)))

@info "Arguments" args
@info "Crew GACCs" crew_gaccs
@info "Fire GACCs" fire_gaccs
@info "Firefighters per crew" firefighters_per_crew
@info "Personnel per crew" personnel_per_crew
@info "Total time limit" time_limit

num_fires = count_selected_fires(fire_gaccs, fires_by_gacc)
num_crews = 0

num_time_periods = 14
travel_speed = 40.0 * 6.0
GC.gc()

crew_routes, fire_plans, crew_models, fire_models, cut_data, init_info = initialize_data_structures(
        num_fires,
        num_crews,
        num_time_periods,
        firefighters_per_crew,
        travel_speed,
        from_empirical = true,
        crew_gaccs = crew_gaccs,
        fire_gaccs = fire_gaccs,
        firefighters_per_crew = firefighters_per_crew,
        initial_firefighters_per_crew = personnel_per_crew,
        fires_by_gacc = fires_by_gacc,
)

num_crews = length(crew_models)
num_fires = length(fire_models)

@info "Total crews" num_crews
@info "Total fires" num_fires
@info "Fire selection criterion" init_info.selection
@info "Fires included in optimization" init_info.fire_ids
@info "Initial crew assignments" init_info.crew_assignments
for (ix, fire_id) in enumerate(init_info.fire_ids)
        start_day = init_info.start_days[ix]
        crew_list = findall(==(ix), init_info.crew_assignments)
        @info "Using fire" ix id=fire_id start_day=start_day initial_crews=crew_list
end
unassigned_crews = findall(==( -1 ), init_info.crew_assignments)
if !isempty(unassigned_crews)
        @info "Crews initially without fire assignment" unassigned_crews
end
for j in 1:num_crews
	no_fire_anticipation!(crew_models[j], [fsp.start_time_period for fsp in fire_models])
end

# add a dummy plan with cost 0 and no crew demands
for fire in 1:num_fires
	@info "adding dummy plan for fire" fire
	new_plan_ix =
		add_column_to_plan_data!(fire_plans, fire, 0.0, zeros(Int64, num_time_periods), Int[])
	if new_plan_ix == -1
		@error "failed to add dummy plan for fire" fire
		error()
	end
	@info "added dummy plan for fire" fire "with index" new_plan_ix
end

for t in 0:num_time_periods

	global crew_routes, fire_plans, crew_models, fire_models, cut_data

	crew_routes = CrewRouteData(Int(floor(6 * 1e6 / num_crews)), num_fires, num_crews, num_time_periods)
	fire_plans = FirePlanData(Int(floor(6 * 1e6  / num_crews)), num_fires, num_time_periods)
	cut_data = CutData(num_crews, num_fires, num_time_periods)

        result = branch_and_price(num_fires,
                num_crews,
                num_time_periods,
                current_time = t,
                from_empirical = true,
                gaccs = crew_gaccs,
                fire_gaccs = fire_gaccs,
                travel_speed = travel_speed,
                firefighters_per_crew = firefighters_per_crew,
                initial_firefighters_per_crew = personnel_per_crew,
                fires_by_gacc = fires_by_gacc,
                crew_routes = crew_routes,
                fire_plans = fire_plans,
                crew_models = crew_models,
                fire_models = fire_models,
                cut_data = cut_data,
                total_time_limit = time_limit,
                output_folder = output_folder,
                )
                # Unpack as many variables as branch_and_price returns, e.g.:
        explored_nodes, ubs, lbs, columns, heuristic_times, times, time_1, root_node_ip_sol, root_node_ip_sol_time, fire_arcs_used, crew_arcs_used = result
        @info "final arcs used" fire_arcs_used, crew_arcs_used

        if fire_arcs_used === nothing || crew_arcs_used === nothing
                @warn "No arc information returned from branch_and_price; stopping early"
                break
        end

        for g in 1:num_fires
                @info "before modify_in_arcs_and_out_arcs!" fire_models[g].state_in_arcs fire_models[g].state_out_arcs fire_arcs_used[g]
                if !isnothing(fire_models[g].start_time_period) && fire_models[g].start_time_period > t
                        @info "fire model start time period is greater than current time, skipping modify_in_arcs_and_out_arcs!" g
                        continue
                end
		modify_in_arcs_and_out_arcs!(fire_models[g], t+1, fire_arcs_used[g], FM.TIME_FROM)
		@info "after modify_in_arcs_and_out_arcs!" fire_models[g].state_in_arcs fire_models[g].state_out_arcs
	end
	for j in 1:num_crews
		@info "before modify_in_arcs_and_out_arcs!" crew_models[j].state_in_arcs crew_models[j].state_out_arcs crew_arcs_used[j]
		modify_in_arcs_and_out_arcs!(crew_models[j], t+1, crew_arcs_used[j], CM.TIME_FROM)
		@info "after modify_in_arcs_and_out_arcs!" crew_models[j].state_in_arcs crew_models[j].state_out_arcs
	end

	# now extract the arc data and costs from the fire_arcs_used and crew_arcs_used and the models
	fire_arcs = Vector{Matrix{Int64}}(undef, num_fires)
	fire_arc_costs = Vector{Vector{Float64}}(undef, num_fires)
	crew_arcs = Vector{Matrix{Int64}}(undef, num_crews)
	crew_arc_costs = Vector{Vector{Float64}}(undef, num_crews)
	for g in 1:num_fires
		fire_arcs[g] = fire_models[g].wide_arcs[:, reverse(fire_arcs_used[g])]
		fire_arc_costs[g] = fire_models[g].arc_costs[reverse(fire_arcs_used[g])]
	end
	for j in 1:num_crews
		crew_arcs[j] = crew_models[j].wide_arcs[:, reverse(crew_arcs_used[j])]
		crew_arc_costs[j] = crew_models[j].arc_costs[reverse(crew_arcs_used[j])]
	end

        # write these to files
        for g in 1:num_fires
                open(joinpath(output_folder, "fire_arcs_$(g)_$(t).json"), "w") do io
                        JSON.print(io, fire_arcs[g])
                end
                open(joinpath(output_folder, "fire_arc_costs_$(g)_$(t).json"), "w") do io
                        JSON.print(io, fire_arc_costs[g])
                end
        end
        for j in 1:num_crews
                open(joinpath(output_folder, "crew_arcs_$(j)_$(t).json"), "w") do io
                        JSON.print(io, crew_arcs[j])
                end
                open(joinpath(output_folder, "crew_arc_costs_$(j)_$(t).json"), "w") do io
                        JSON.print(io, crew_arc_costs[j])
                end
        end
end

close(log_file)

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
