include("BranchAndPrice.jl")

using JuMP, Gurobi, JSON, Profile, ArgParse, Logging, IterTools, CSV, DataFrames, Dates
import DataFrames: groupby
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
        "SA" => "Southern Area",
        "SC" => "Southern California",
        "CA-S" => "Southern California",
        "SW" => "Southwest",
)

const ALL_GACCS = unique(collect(values(GACC_ABBR)))

const JULIA_EXT_STATE_CODE = 1

function default_discretization_bins()
        step_sizes = (5.0, 15.0, 30.0, 50.0, 100.0)
        bins = Float64[]
        append!(bins, collect(1.0:step_sizes[0]:100.0))
        append!(bins, collect(100.0:step_sizes[1]:2000.0))
        append!(bins, collect(2000.0:step_sizes[2]:5000.0))
        append!(bins, collect(5000.0:step_sizes[3]:10000.0))
        append!(bins, collect(10000.0:step_sizes[4]:100000.0))
        unique_bins = sort(unique(bins))
        if unique_bins[end] < 100000.0
                push!(unique_bins, 100000.0)
        end
        return unique_bins
end

function decode_packed_state_area(packed_code::Int, bins::Vector{Float64})
        packed_code <= JULIA_EXT_STATE_CODE && return 0.0
        num_bins = length(bins)
        num_bins == 0 && return 0.0
        active = packed_code - (JULIA_EXT_STATE_CODE + 1)
        next_idx = active ÷ num_bins
        next_idx = clamp(next_idx, 0, num_bins - 1)
        return bins[next_idx + 1]
end

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

function count_selected_fires(
        fire_gaccs::Vector{String},
        fires_by_gacc::Dict{String,Vector{Int64}},
        input_folder::String,
)
        selected_fires = CSV.read(joinpath(input_folder, "selected_fires.csv"), DataFrame)
        selected_fires[!, "GACC"] = normalize_gacc.(selected_fires[!, "GACC"])
        fire_gaccs = normalize_gacc.(fire_gaccs)
        if !isempty(fires_by_gacc)
                normalized_fires_by_gacc = Dict{String,Vector{Int64}}()
                for (gacc, fires) in fires_by_gacc
                        normalized_fires_by_gacc[normalize_gacc(gacc)] = fires
                end
                mask = falses(nrow(selected_fires))
                for (gacc, fires) in normalized_fires_by_gacc
                        mask .|= (selected_fires[:, "GACC"] .== gacc) .& in.(selected_fires[:, "FIRE_EVENT_ID"], Ref(fires))
                end
                selected_fires = selected_fires[mask, :]
        else
                selected_fires = selected_fires[in.(selected_fires[:, "GACC"], Ref(fire_gaccs)), :]
        end
        return length(unique(selected_fires[:, "FIRE_EVENT_ID"]))
end

function build_day_one_fire_subset(
        fire_gaccs::Vector{String},
        fires_by_gacc::Dict{String,Vector{Int64}},
        input_folder::String,
)
        selected_fires = CSV.read(joinpath(input_folder, "selected_fires.csv"), DataFrame)
        # normalize column names to handle legacy exports with different casing
        name_lookup = Dict(lowercase(String(col)) => col for col in names(selected_fires))
        start_key = "start_day_of_sim"
        alt_keys = ("sim_start_day_dsfr", "day_since_first_report", "start_day")
        if haskey(name_lookup, start_key)
                start_col = name_lookup[start_key]
        else
                start_col = nothing
                for key in alt_keys
                        if haskey(name_lookup, key)
                                start_col = name_lookup[key]
                                break
                        end
                end
                isnothing(start_col) && error("selected_fires.csv missing 'start_day_of_sim' column required for --day-1-only")
                selected_fires[!, :start_day_of_sim] = copy(selected_fires[!, start_col])
        end
        if start_col !== :start_day_of_sim
                selected_fires[!, :start_day_of_sim] = copy(selected_fires[!, start_col])
        end
        selected_fires[!, :GACC] = normalize_gacc.(selected_fires[!, :GACC])

        if !isempty(fires_by_gacc)
                normalized = Dict{String,Set{Int64}}()
                for (gacc, fire_list) in fires_by_gacc
                        normalized[normalize_gacc(gacc)] = Set(Int64.(fire_list))
                end
                mask = falses(nrow(selected_fires))
                for (gacc, fire_ids) in normalized
                        mask .|= (selected_fires[:, :GACC] .== gacc) .&
                                in.(selected_fires[:, :FIRE_EVENT_ID], Ref(fire_ids))
                end
        else
                normalized_gaccs = normalize_gacc.(fire_gaccs)
                mask = in.(selected_fires[:, :GACC], Ref(normalized_gaccs))
        end

        mask .&= (selected_fires[:, :start_day_of_sim] .== 0)
        filtered = selected_fires[mask, :]

        subset = Dict{String,Vector{Int64}}()
        if nrow(filtered) == 0
                return subset
        end

        for subdf in groupby(filtered, :GACC)
                gacc = subdf[1, :GACC]
                subset[gacc] = collect(unique(Int64.(subdf[!, :FIRE_EVENT_ID])))
        end

        return subset
end

function resolve_input_folder(folder::String)
        # if the provided path already points to a folder with selected_fires.csv, use it
        if isfile(joinpath(folder, "selected_fires.csv"))
                return folder
        end
        # check for an arc_arrays subdirectory
        arc_path = joinpath(folder, "arc_arrays")
        if isfile(joinpath(arc_path, "selected_fires.csv"))
                return arc_path
        end
        # otherwise look relative to data/empirical_fire_models
        candidate = joinpath("data", "empirical_fire_models", folder, "arc_arrays")
        if isfile(joinpath(candidate, "selected_fires.csv"))
                return candidate
        end
        error("Input folder $folder not found or missing selected_fires.csv")
end

function get_command_line_args()
        arg_parse_settings = ArgParseSettings()
        @add_arg_table arg_parse_settings begin
                "--debug"
                help = "run in debug mode, exposing all logging that uses @debug macro"
                action = :store_true
                "--seed-fastest-extinguish"
                help = "seed a plan that reaches extinguishment as early as possible for each active fire"
                action = :store_true
                "--seed-frontloaded"
                help = "seed a plan that uses the maximum feasible crews on the first active day for each fire"
                action = :store_true
                "--seed-max-daily"
                help = "seed a plan that uses the maximum feasible crews on every remaining day (greedy per-day max) for each fire"
                action = :store_true
                "--seed-crew-routes"
                help = "seed crew routes that arrive at early target times so heavy early fire plans are feasible"
                action = :store_true
                "--final-snapshot-only"
                help = "optimize only the end-of-horizon area (final snapshot) for each fire (last nonzero pre-extinguish area)"
                action = :store_true
                "--day-1-only"
                help = "restrict inputs to fires that start on day 1 of the planning window (start_day_of_sim == 0)"
                action = :store_true
                "--crew-costs"
                help = "crew cost mode: 'on' (default) or 'off' to ignore crew travel/rest costs in the objective"
                default = "on"
                "--crew-gaccs"
                help = "Allowed crew GACCs: 'all', 'all_no_ak', or comma-separated list of abbreviations (e.g. GB,NW)"
                default = "GB"
                "--fire-gaccs"
                help = "Allowed fire GACCs: 'all', 'all_no_ak', or comma-separated list of abbreviations (e.g. GB,NW). Defaults to crew GACCs if omitted"
                default = ""
                "--firefighters-per-crew"
                help = "Number of firefighters per crew"
                arg_type = Int
                default = 50
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
                "--input-folder"
                help = "Directory containing input files (full path or dataset under data/empirical_fire_models)"
                default = "raw"
                "--output-folder"
                help = "Directory to store output files"
                default = "data/output"
        end
        return parse_args(arg_parse_settings)
end


args = get_command_line_args()
crew_gaccs = parse_gaccs(args["crew-gaccs"])
fire_gaccs = isempty(args["fire-gaccs"]) ? crew_gaccs : parse_gaccs(args["fire-gaccs"])
seed_fastest_extinguish = args["seed-fastest-extinguish"]
seed_frontloaded = args["seed-frontloaded"]
seed_max_daily = args["seed-max-daily"]
seed_crew_routes = args["seed-crew-routes"]
final_snapshot_only = args["final-snapshot-only"]
day_one_only = args["day-1-only"]
crew_costs_mode = lowercase(String(args["crew-costs"]))
zero_crew_costs = (crew_costs_mode in ("off","0","false","no"))
firefighters_per_crew = args["firefighters-per-crew"]
personnel_per_crew = args["personnel-per-crew"]
fires_by_gacc = parse_fires_by_gacc(args["fires"])
time_limit = args["time-limit"]
input_folder = resolve_input_folder(args["input-folder"])
output_folder = args["output-folder"]

if day_one_only
        day_one_subset = build_day_one_fire_subset(fire_gaccs, fires_by_gacc, input_folder)
        isempty(day_one_subset) && error("No fires begin on day 1 after applying --day-1-only filter.")

        normalized_request_order = normalize_gacc.(fire_gaccs)
        filtered_order = [g for g in normalized_request_order if haskey(day_one_subset, g)]
        for gacc in keys(day_one_subset)
                gacc ∈ filtered_order || push!(filtered_order, gacc)
        end

        fire_gaccs = filtered_order
        fires_by_gacc = Dict(g => sort(day_one_subset[g]) for g in keys(day_one_subset))
        total_day_one_fires = sum(length(ids) for ids in values(fires_by_gacc))
        @info "--day-1-only filter active" total_fires=total_day_one_fires gaccs=fire_gaccs
end

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
@info "Seed fastest-extinguish" seed_fastest_extinguish
@info "Seed frontloaded" seed_frontloaded
@info "Seed max-daily" seed_max_daily
@info "Seed crew routes" seed_crew_routes
@info "Final snapshot only" final_snapshot_only
@info "Day 1 only" day_one_only
@info "Crew costs mode" (zero_crew_costs ? "off" : "on")

num_fires = count_selected_fires(fire_gaccs, fires_by_gacc, input_folder)
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
        input_folder = input_folder,
        from_empirical = true,
        crew_gaccs = crew_gaccs,
        fire_gaccs = fire_gaccs,
        firefighters_per_crew = firefighters_per_crew,
        initial_firefighters_per_crew = personnel_per_crew,
        fires_by_gacc = fires_by_gacc,
        sorted_fire_output_folder = output_folder,
        zero_crew_costs = zero_crew_costs,
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

# Track committed arcs (history) so past-day decisions are preserved across re-solves
committed_fire_arcs = [Set{Int64}() for _ in 1:num_fires]
committed_crew_arcs = [Set{Int64}() for _ in 1:num_crews]

let prev_dual_warm_start = nothing
optimizer_day_rollup = Any[]

for t in 0:num_time_periods

    global crew_routes, fire_plans, crew_models, fire_models, cut_data

    crew_routes = CrewRouteData(Int(floor(6 * 1e6 / num_crews)), num_fires, num_crews, num_time_periods)
    fire_plans = FirePlanData(Int(floor(6 * 1e6  / num_crews)), num_fires, num_time_periods)
	cut_data = CutData(num_crews, num_fires, num_time_periods)

    # seed dummy plan/route columns so the restricted master is always feasible
    dummy_plan_cost = 1.0e8
    dummy_route_cost = 1.0e8
    for fire in 1:num_fires
        add_column_to_plan_data!(
            fire_plans,
            fire,
            dummy_plan_cost,
            zeros(Int64, num_time_periods),
            Int[],
        )
    end
    for crew in 1:num_crews
        fires_fought = falses(num_fires, num_time_periods)
        add_column_to_route_data!(
            crew_routes,
            crew,
            dummy_route_cost,
            fires_fought,
            Int[],
        )
    end

    current_day = t + 1
    total_days = num_time_periods + 1
    # determine which fires are active or starting this day for user feedback
    fire_start_periods = [fsp.start_time_period for fsp in fire_models]
    fires_active = [g for g in 1:num_fires if isnothing(fire_start_periods[g]) || fire_start_periods[g] <= current_day]
    fires_starting_today = [g for g in 1:num_fires if fire_start_periods[g] == current_day]

    @info "##### Simulation Day $(current_day) of $(total_days) #####"
    @info "Active fires" length(fires_active) fires_active
    if !isempty(fires_starting_today)
        @info "Fires starting today" fires_starting_today
    end
    @info "Fires within planning horizon" day=current_day window_start=current_day - 1 window_end=(current_day - 1 + num_time_periods - 1) fires=[
        (
            optimizer_index = g,
            fire_id = init_info.fire_ids[g],
            start_day = init_info.start_days[g],
            start_date = begin
                base_label = replace(basename(input_folder), "fire_models_" => "")
                base_date = try
                    Date(base_label, dateformat"yyyy_mm_dd")
                catch
                    nothing
                end
                isnothing(base_date) ? nothing : Dates.format(base_date + Day(init_info.start_days[g]), dateformat"yyyy-mm-dd")
            end
        )
        for g in 1:num_fires
        if init_info.start_days[g] <= (current_day - 1 + num_time_periods - 1)
    ]

        # Helper to compute plan cost per mode
        compute_plan_cost = function(fm, path_chrono::Vector{Int})
            arcs_used = path_chrono
            if final_snapshot_only
                best_arc = 0
                best_t = -1
                for a_ix in arcs_used
                    to_t = fm.long_arcs[a_ix, FM.TIME_TO]
                    c = fm.arc_costs[a_ix]
                    if (to_t - 1) <= num_time_periods && c > 1e-12 && (to_t - 1) > best_t
                        best_t = to_t - 1
                        best_arc = a_ix
                    end
                end
                return best_arc == 0 ? sum(fm.arc_costs[arcs_used]) : fm.arc_costs[best_arc]
            else
                return sum(fm.arc_costs[arcs_used])
            end
        end

        # Seeding: fastest extinguish
        if seed_fastest_extinguish
            for g in fires_active
                fm = fire_models[g]
                states, times = size(fm.state_in_arcs)
                # Earliest reachability to extinguish state (assume state 1)
                prev_arc = fill(0, states, times)
                reachable = falses(states, times)
                ext_state_id = 1
                found_t = nothing
                for tt in 1:times
                    for s in 1:states
                        for arc_ix in fm.state_in_arcs[s, tt]
                            arc = fm.long_arcs[arc_ix, :]
                            tf = arc[FM.TIME_FROM]
                            sf = arc[FM.STATE_FROM]
                            if (tf == 0) || (tf >= 1 && reachable[sf, tf])
                                if prev_arc[s, tt] == 0
                                    prev_arc[s, tt] = arc_ix
                                    reachable[s, tt] = true
                                end
                            end
                        end
                    end
                    if reachable[ext_state_id, tt]
                        found_t = tt
                        break
                    end
                end
                if isnothing(found_t)
                    continue
                end
                # Backtrack to build path
                cur_s = ext_state_id
                cur_t = found_t
                path_rev = Int[]
                while cur_t != 0
                    arc_ix = prev_arc[cur_s, cur_t]
                    if arc_ix == 0
                        break
                    end
                    push!(path_rev, arc_ix)
                    arc = fm.long_arcs[arc_ix, :]
                    cur_s = arc[FM.STATE_FROM]
                    cur_t = arc[FM.TIME_FROM]
                end
                path_chrono = reverse(path_rev)
                if isempty(path_chrono)
                    continue
                end
                # Seeding guard: if fire is at/after its start day, require at least one post-start arc
                start_day = fire_models[g].start_time_period
                if !isnothing(start_day) && current_day >= start_day
                    has_post_start = false
                    for a_ix in path_chrono
                        if fm.long_arcs[a_ix, FM.TIME_FROM] >= start_day + 1
                            has_post_start = true
                            break
                        end
                    end
                    if !has_post_start
                        continue
                    end
                end
                # Build plan data
                cost = compute_plan_cost(fm, path_chrono)
                crew_demands = zeros(Int, num_time_periods)
                for a_ix in path_chrono
                    tf = fm.long_arcs[a_ix, FM.TIME_FROM]
                    if 1 <= tf <= num_time_periods
                        crew_demands[tf] = fm.long_arcs[a_ix, FM.CREWS_PRESENT]
                    end
                end
                add_column_to_plan_data!(fire_plans, g, cost, crew_demands, reverse(path_chrono))
            end
        end

        # Seeding: frontloaded (max crews on first active day, then min-cost continuation)
        if seed_frontloaded
            for g in fires_active
                fm = fire_models[g]
                states, times = size(fm.state_in_arcs)
                # Build reachability to current_day
                prev_arc = fill(0, states, times)
                reachable = falses(states, times)
                for tt in 1:current_day
                    for s in 1:states
                        for arc_ix in fm.state_in_arcs[s, tt]
                            arc = fm.long_arcs[arc_ix, :]
                            tf = arc[FM.TIME_FROM]
                            sf = arc[FM.STATE_FROM]
                            if (tf == 0) || (tf >= 1 && reachable[sf, tf])
                                if prev_arc[s, tt] == 0
                                    prev_arc[s, tt] = arc_ix
                                    reachable[s, tt] = true
                                end
                            end
                        end
                    end
                end
                # Pick max crews at current_day
                start_arc = 0
                best_crews = -1
                best_cost = Inf
                best_from_state = 0
                best_to_state = 0
                if current_day <= size(fm.state_out_arcs,2)
                    for s in 1:states
                        if !reachable[s, current_day]
                            continue
                        end
                        for arc_ix in fm.state_out_arcs[s, current_day]
                            crews_here = fm.long_arcs[arc_ix, FM.CREWS_PRESENT]
                            c = fm.arc_costs[arc_ix]
                            if (crews_here > best_crews) || (crews_here == best_crews && c < best_cost)
                                best_crews = crews_here
                                best_cost = c
                                start_arc = arc_ix
                                best_from_state = s
                                best_to_state = fm.long_arcs[arc_ix, FM.STATE_TO]
                            end
                        end
                    end
                end
                if start_arc == 0
                    continue
                end
                # Backtrack prefix
                cur_s = best_from_state
                cur_t = current_day
                prefix_rev = Int[]
                while cur_t != 0
                    arc_ix = prev_arc[cur_s, cur_t]
                    if arc_ix == 0
                        break
                    end
                    push!(prefix_rev, arc_ix)
                    arc = fm.long_arcs[arc_ix, :]
                    cur_s = arc[FM.STATE_FROM]
                    cur_t = arc[FM.TIME_FROM]
                end
                path_chrono = vcat(reverse(prefix_rev), [start_arc])
                # Greedy min-cost continuation
                cur_state = best_to_state
                cur_time = fm.long_arcs[start_arc, FM.TIME_TO]
                while cur_time <= num_time_periods
                    picked = 0
                    picked_cost = Inf
                    if cur_state >= 1 && cur_state <= size(fm.state_out_arcs,1) && cur_time >= 1 && cur_time <= size(fm.state_out_arcs,2)
                        for arc_ix in fm.state_out_arcs[cur_state, cur_time]
                            c = fm.arc_costs[arc_ix]
                            if c < picked_cost
                                picked = arc_ix
                                picked_cost = c
                            end
                        end
                    end
                    if picked == 0
                        break
                    end
                    push!(path_chrono, picked)
                    cur_state = fm.long_arcs[picked, FM.STATE_TO]
                    cur_time = fm.long_arcs[picked, FM.TIME_TO]
                    if cur_time == 0
                        break
                    end
                end
                if isempty(path_chrono)
                    continue
                end
                # Seeding guard: if fire is at/after its start day, require at least one post-start arc
                start_day = fire_models[g].start_time_period
                if !isnothing(start_day) && current_day >= start_day
                    has_post_start = false
                    for a_ix in path_chrono
                        if fm.long_arcs[a_ix, FM.TIME_FROM] >= start_day + 1
                            has_post_start = true
                            break
                        end
                    end
                    if !has_post_start
                        continue
                    end
                end
                cost = compute_plan_cost(fm, path_chrono)
                crew_demands = zeros(Int, num_time_periods)
                for a_ix in path_chrono
                    tf = fm.long_arcs[a_ix, FM.TIME_FROM]
                    if 1 <= tf <= num_time_periods
                        crew_demands[tf] = fm.long_arcs[a_ix, FM.CREWS_PRESENT]
                    end
                end
                add_column_to_plan_data!(fire_plans, g, cost, crew_demands, reverse(path_chrono))
            end
        end

        # Seeding: max-daily (greedy per-day maximum crews)
        if seed_max_daily
            for g in fires_active
                fm = fire_models[g]
                states, times = size(fm.state_in_arcs)
                # Build reachability to current_day
                prev_arc = fill(0, states, times)
                reachable = falses(states, times)
                for tt in 1:current_day
                    for s in 1:states
                        for arc_ix in fm.state_in_arcs[s, tt]
                            arc = fm.long_arcs[arc_ix, :]
                            tf = arc[FM.TIME_FROM]
                            sf = arc[FM.STATE_FROM]
                            if (tf == 0) || (tf >= 1 && reachable[sf, tf])
                                if prev_arc[s, tt] == 0
                                    prev_arc[s, tt] = arc_ix
                                    reachable[s, tt] = true
                                end
                            end
                        end
                    end
                end
                # Select start arc at current_day with max crews (tie-break by low cost)
                start_arc = 0
                best_crews = -1
                best_cost = Inf
                start_from_state = 0
                start_to_state = 0
                if current_day <= size(fm.state_out_arcs,2)
                    for s in 1:states
                        if !reachable[s, current_day]
                            continue
                        end
                        for arc_ix in fm.state_out_arcs[s, current_day]
                            crews_here = fm.long_arcs[arc_ix, FM.CREWS_PRESENT]
                            c = fm.arc_costs[arc_ix]
                            if (crews_here > best_crews) || (crews_here == best_crews && c < best_cost)
                                best_crews = crews_here
                                best_cost = c
                                start_arc = arc_ix
                                start_from_state = s
                                start_to_state = fm.long_arcs[arc_ix, FM.STATE_TO]
                            end
                        end
                    end
                end
                if start_arc == 0
                    continue
                end
                # Backtrack prefix
                cur_s = start_from_state
                cur_t = current_day
                prefix_rev = Int[]
                while cur_t != 0
                    arc_ix = prev_arc[cur_s, cur_t]
                    if arc_ix == 0
                        break
                    end
                    push!(prefix_rev, arc_ix)
                    arc = fm.long_arcs[arc_ix, :]
                    cur_s = arc[FM.STATE_FROM]
                    cur_t = arc[FM.TIME_FROM]
                end
                # Greedy max crews continuation
                path_chrono = vcat(reverse(prefix_rev), [start_arc])
                cur_state = start_to_state
                cur_time = fm.long_arcs[start_arc, FM.TIME_TO]
                while cur_time <= num_time_periods
                    picked = 0
                    picked_crews = -1
                    picked_cost = Inf
                    if cur_state >= 1 && cur_state <= size(fm.state_out_arcs,1) && cur_time >= 1 && cur_time <= size(fm.state_out_arcs,2)
                        for arc_ix in fm.state_out_arcs[cur_state, cur_time]
                            crews_here = fm.long_arcs[arc_ix, FM.CREWS_PRESENT]
                            c = fm.arc_costs[arc_ix]
                            if (crews_here > picked_crews) || (crews_here == picked_crews && c < picked_cost)
                                picked = arc_ix
                                picked_crews = crews_here
                                picked_cost = c
                            end
                        end
                    end
                    if picked == 0
                        break
                    end
                    push!(path_chrono, picked)
                    cur_state = fm.long_arcs[picked, FM.STATE_TO]
                    cur_time = fm.long_arcs[picked, FM.TIME_TO]
                    if cur_time == 0
                        break
                    end
                end
                if isempty(path_chrono)
                    continue
                end
                # Seeding guard: if fire is at/after its start day, require at least one post-start arc
                start_day = fire_models[g].start_time_period
                if !isnothing(start_day) && current_day >= start_day
                    has_post_start = false
                    for a_ix in path_chrono
                        if fm.long_arcs[a_ix, FM.TIME_FROM] >= start_day + 1
                            has_post_start = true
                            break
                        end
                    end
                    if !has_post_start
                        continue
                    end
                end
                cost = compute_plan_cost(fm, path_chrono)
                crew_demands = zeros(Int, num_time_periods)
                for a_ix in path_chrono
                    tf = fm.long_arcs[a_ix, FM.TIME_FROM]
                    if 1 <= tf <= num_time_periods
                        crew_demands[tf] = fm.long_arcs[a_ix, FM.CREWS_PRESENT]
                    end
                end
                add_column_to_plan_data!(fire_plans, g, cost, crew_demands, reverse(path_chrono))
            end
        end

        # Additional seeding: frontload explicitly at the first positive day (start_day+1)
        # This ensures a heavy plan exists exactly at the first actionable day even when current_day < start_day+1
        if seed_frontloaded
            for g in fires_active
                fm = fire_models[g]
                start_day = fm.start_time_period
                if isnothing(start_day)
                    continue
                end
                target_day = start_day + 1
                if target_day < 1 || target_day > num_time_periods
                    continue
                end
                # Build reachability up to target_day
                states, times = size(fm.state_in_arcs)
                prev_arc = fill(0, states, times)
                reachable = falses(states, times)
                for tt in 1:target_day
                    for s in 1:states
                        for arc_ix in fm.state_in_arcs[s, tt]
                            arc = fm.long_arcs[arc_ix, :]
                            tf = arc[FM.TIME_FROM]
                            sf = arc[FM.STATE_FROM]
                            if (tf == 0) || (tf >= 1 && reachable[sf, tf])
                                if prev_arc[s, tt] == 0
                                    prev_arc[s, tt] = arc_ix
                                    reachable[s, tt] = true
                                end
                            end
                        end
                    end
                end
                # Pick max crews at target_day
                start_arc = 0
                best_crews = -1
                best_cost = Inf
                best_from_state = 0
                best_to_state = 0
                if target_day <= size(fm.state_out_arcs,2)
                    for s in 1:states
                        if !reachable[s, target_day]
                            continue
                        end
                        for arc_ix in fm.state_out_arcs[s, target_day]
                            crews_here = fm.long_arcs[arc_ix, FM.CREWS_PRESENT]
                            c = fm.arc_costs[arc_ix]
                            if (crews_here > best_crews) || (crews_here == best_crews && c < best_cost)
                                best_crews = crews_here
                                best_cost = c
                                start_arc = arc_ix
                                best_from_state = s
                                best_to_state = fm.long_arcs[arc_ix, FM.STATE_TO]
                            end
                        end
                    end
                end
                if start_arc == 0
                    continue
                end
                # Backtrack prefix to build path to target_day
                cur_s = best_from_state
                cur_t = target_day
                prefix_rev = Int[]
                while cur_t != 0
                    arc_ix = prev_arc[cur_s, cur_t]
                    if arc_ix == 0
                        break
                    end
                    push!(prefix_rev, arc_ix)
                    arc = fm.long_arcs[arc_ix, :]
                    cur_s = arc[FM.STATE_FROM]
                    cur_t = arc[FM.TIME_FROM]
                end
                path_chrono = vcat(reverse(prefix_rev), [start_arc])
                # Greedy min-cost continuation from target_day onward
                cur_state = best_to_state
                cur_time = fm.long_arcs[start_arc, FM.TIME_TO]
                while cur_time <= num_time_periods
                    picked = 0
                    picked_cost = Inf
                    if cur_state >= 1 && cur_state <= size(fm.state_out_arcs,1) && cur_time >= 1 && cur_time <= size(fm.state_out_arcs,2)
                        for arc_ix in fm.state_out_arcs[cur_state, cur_time]
                            c = fm.arc_costs[arc_ix]
                            if c < picked_cost
                                picked = arc_ix
                                picked_cost = c
                            end
                        end
                    end
                    if picked == 0
                        break
                    end
                    push!(path_chrono, picked)
                    cur_state = fm.long_arcs[picked, FM.STATE_TO]
                    cur_time = fm.long_arcs[picked, FM.TIME_TO]
                    if cur_time == 0
                        break
                    end
                end
                if isempty(path_chrono)
                    continue
                end
                # Build and add plan
                cost = compute_plan_cost(fm, path_chrono)
                crew_demands = zeros(Int, num_time_periods)
                for a_ix in path_chrono
                    tf = fm.long_arcs[a_ix, FM.TIME_FROM]
                    if 1 <= tf <= num_time_periods
                        crew_demands[tf] = fm.long_arcs[a_ix, FM.CREWS_PRESENT]
                    end
                end
                add_column_to_plan_data!(fire_plans, g, cost, crew_demands, reverse(path_chrono))
            end
        end

        warm_start_to_use = prev_dual_warm_start
        if !isnothing(warm_start_to_use)
            dims = size(warm_start_to_use.linking_values)
            if dims[1] != num_fires || dims[2] != num_time_periods
                @info "Skipping dual warm start due to dimension mismatch" warm_dims = dims current_dims = (num_fires, num_time_periods)
                warm_start_to_use = nothing
            else
                @info "Reusing dual warm start" warm_dims = dims
            end
        end

        # Seed early crew routes to ensure feasible supply at first active days
        # Only if enabled via --seed-crew-routes
        if seed_crew_routes
            target_times = Int[]
            push!(target_times, current_day)
            if current_day + 1 <= num_time_periods
                push!(target_times, current_day + 1)
            end
            for g in fires_active
                for j in 1:num_crews
                    long_arcs = crew_models[j].long_arcs
                    for tt in target_times
                        # find arcs that arrive to fire g at time tt
                        cand_ixs = [
                            i for i in 1:size(long_arcs, 1) if
                            (long_arcs[i, CM.TO_TYPE] == CM.FIRE_CODE) &&
                            (long_arcs[i, CM.LOC_TO] == g) &&
                            (long_arcs[i, CM.TIME_TO] == tt)
                        ]
                        if isempty(cand_ixs)
                            continue
                        end
                        # build a modified cost vector that strongly rewards hitting (g,tt)
                        modified_costs = copy(crew_models[j].arc_costs)
                        bonus = 1.0e6
                        for ix in cand_ixs
                            modified_costs[ix] -= bonus
                        end
                        prohibited = falses(length(modified_costs))
                        obj, arcs_used = crew_dp_subproblem(
                            crew_models[j].wide_arcs,
                            modified_costs,
                            prohibited,
                            crew_models[j].state_in_arcs,
                        )
                        if isempty(arcs_used)
                            continue
                        end
                        # ensure the selected path actually includes one of the target arrival arcs
                        contains_target = any(a -> a in cand_ixs, arcs_used)
                        if !contains_target
                            continue
                        end
                        # Compute true static cost and fires_fought mask for this path
                        true_cost = sum(crew_models[j].arc_costs[arcs_used])
                        fires_fought = get_fires_fought(
                            crew_models[j].wide_arcs,
                            arcs_used,
                            (num_fires, num_time_periods),
                        )
                        add_column_to_route_data!(
                            crew_routes,
                            j,
                            true_cost,
                            fires_fought,
                            arcs_used,
                        )
                    end
                end
            end
        end

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
                input_folder = input_folder,
                crew_routes = crew_routes,
                fire_plans = fire_plans,
                crew_models = crew_models,
                fire_models = fire_models,
                cut_data = cut_data,
                total_time_limit = time_limit,
                output_folder = output_folder,
                dual_warm_start = warm_start_to_use,
                final_snapshot_only = final_snapshot_only,
                )
                # Unpack as many variables as branch_and_price returns, e.g.:
        explored_nodes, ubs, lbs, columns, heuristic_times, times, time_1, root_node_ip_sol, root_node_ip_sol_time, fire_arcs_used, crew_arcs_used, root_dual_warm_start = result
        @debug "final arcs used" fire_arcs_used, crew_arcs_used

        if fire_arcs_used === nothing || crew_arcs_used === nothing
                @warn "No arc information returned from branch_and_price; stopping early"
                break
        end

        prev_dual_warm_start = root_dual_warm_start

        # Commit decisions up to current day (t+1): preserve those arcs in future iterations
        commit_cutoff = t + 1
        for g in 1:num_fires
                for a_ix in fire_arcs_used[g]
                        tf = fire_models[g].long_arcs[a_ix, FM.TIME_FROM]
                        if tf < commit_cutoff
                                push!(committed_fire_arcs[g], a_ix)
                        end
                end
        end
        for j in 1:num_crews
                for a_ix in crew_arcs_used[j]
                        tf = crew_models[j].long_arcs[a_ix, CM.TIME_FROM]
                        if tf < commit_cutoff
                                push!(committed_crew_arcs[j], a_ix)
                        end
                end
        end

        for g in 1:num_fires
                @debug "before modify_in_arcs_and_out_arcs!" fire_models[g].state_in_arcs fire_models[g].state_out_arcs fire_arcs_used[g]
                if !isnothing(fire_models[g].start_time_period) && fire_models[g].start_time_period > t
                        @debug "fire model start time period is greater than current time, skipping modify_in_arcs_and_out_arcs!" g
                        continue
                end
		# Preserve both committed history and the current solution’s arcs
		begin
		    local keep_arcs = unique(vcat(collect(committed_fire_arcs[g]), fire_arcs_used[g]))
		    modify_in_arcs_and_out_arcs!(fire_models[g], t+1, keep_arcs, FM.TIME_FROM)
		end
		@debug "after modify_in_arcs_and_out_arcs!" fire_models[g].state_in_arcs fire_models[g].state_out_arcs
	end
	for j in 1:num_crews
		@debug "before modify_in_arcs_and_out_arcs!" crew_models[j].state_in_arcs crew_models[j].state_out_arcs crew_arcs_used[j]
		begin
		    local keep_arcs = unique(vcat(collect(committed_crew_arcs[j]), crew_arcs_used[j]))
		    modify_in_arcs_and_out_arcs!(crew_models[j], t+1, keep_arcs, CM.TIME_FROM)
		end
		@debug "after modify_in_arcs_and_out_arcs!" crew_models[j].state_in_arcs crew_models[j].state_out_arcs
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

        discretization_bins_used = (:discretization_bins in propertynames(init_info)) ? init_info.discretization_bins : nothing
        bins_for_decoding = discretization_bins_used === nothing ? default_discretization_bins() : collect(Float64.(discretization_bins_used))
        template_fire_entries = (:fire_order in propertynames(init_info)) ? [deepcopy(entry) for entry in init_info.fire_order] : [Dict{String,Any}("optimizer_index" => g) for g in 1:num_fires]
        fire_manifest_entries = copy(template_fire_entries)
        crew_manifest_entries = Vector{Dict{String,Any}}()
        day_summary_entries = String[]

        for g in 1:num_fires
                arcs_filename = "fire_arcs_$(g)_$(t).json"
                arc_costs_filename = "fire_arc_costs_$(g)_$(t).json"

                state_entry = fire_manifest_entries[g]
                state_meta = get(state_entry, "state_metadata", nothing)
                state_lookup = Dict{Int,Dict{String,Any}}()
                state_area_map_sim = Dict{Int, Float64}()
                state_area_map_discrete = Dict{Int, Float64}()
                packed_lookup = Dict{Int, Int}()
                if state_meta isa AbstractVector
                        for sm in state_meta
                                if sm isa Dict && haskey(sm, "state_id")
                                        sid = Int(sm["state_id"])
                                        state_lookup[sid] = sm
                                        if haskey(sm, "area_acres_sim") && sm["area_acres_sim"] !== nothing
                                                state_area_map_sim[sid] = Float64(sm["area_acres_sim"])
                                        end
                                        if haskey(sm, "area_acres_discrete") && sm["area_acres_discrete"] !== nothing
                                                state_area_map_discrete[sid] = Float64(sm["area_acres_discrete"])
                                        end
                                        if haskey(sm, "packed_code") && sm["packed_code"] !== nothing
                                                try
                                                        packed_lookup[sid] = Int(sm["packed_code"])
                                                catch
                                                end
                                        end
                                end
                        end
                end

                fire_arcs_export = fire_arcs[g]
                if !isempty(state_lookup)
                        fire_arcs_export = copy(fire_arcs[g])
                        map_state = function(state_idx::Int)
                                sm = get(state_lookup, state_idx, nothing)
                                if sm isa Dict && haskey(sm, "packed_code")
                                        raw = sm["packed_code"]
                                        if !(raw === nothing || raw === missing)
                                                return Int(raw)
                                        end
                                end
                                return state_idx
                        end
                        fire_arcs_export[FM.STATE_FROM, :] = map(map_state, fire_arcs_export[FM.STATE_FROM, :])
                        fire_arcs_export[FM.STATE_TO, :] = map(map_state, fire_arcs_export[FM.STATE_TO, :])
                end

                open(joinpath(output_folder, arcs_filename), "w") do io
                        JSON.print(io, fire_arcs_export)
                end
                open(joinpath(output_folder, arc_costs_filename), "w") do io
                        JSON.print(io, fire_arc_costs[g])
                end

                state_area_map = state_area_map_sim
                state_area_map_discrete = isempty(state_area_map_discrete) ? Dict{Int,Float64}() : state_area_map_discrete

                num_periods = num_time_periods
                daily_crews = fill(0, num_periods)
                daily_area = Vector{Union{Nothing, Float64}}(undef, num_periods)
                daily_area_discrete = Vector{Union{Nothing, Float64}}(undef, num_periods)
                running_area = 0.0
                running_area_discrete = 0.0
                for i in 1:num_periods
                        daily_area[i] = nothing
                        daily_area_discrete[i] = nothing
                end

                # Fill committed prefix first (days < current_day), then plan suffix (≥ current_day)
                for arc_idx in collect(committed_fire_arcs[g])
                        arc_row = fire_models[g].long_arcs[arc_idx, :]
                        time_from = arc_row[FM.TIME_FROM]
                        time_to = arc_row[FM.TIME_TO]
                        crew_count = arc_row[FM.CREWS_PRESENT]
                        state_to = arc_row[FM.STATE_TO]

                        if 1 ≤ time_from ≤ num_periods && time_from < current_day
                                daily_crews[time_from] = crew_count
                        end

                        period_ix = time_to - 1
                        packed_code = get(packed_lookup, state_to, nothing)
                        area_val = nothing
                        if packed_code !== nothing
                                area_val = decode_packed_state_area(packed_code, bins_for_decoding)
                        elseif haskey(state_area_map_sim, state_to)
                                area_val = state_area_map_sim[state_to]
                        end
                        area_val_discrete = if haskey(state_area_map_discrete, state_to)
                                state_area_map_discrete[state_to]
                            elseif !isnothing(discretization_bins_used) && state_to ≥ 1 && state_to ≤ length(discretization_bins_used)
                                discretization_bins_used[state_to]
                            else
                                nothing
                            end
                        if !isnothing(area_val) && 1 ≤ period_ix ≤ num_periods
                                running_area = max(running_area, area_val)
                                daily_area[period_ix] = running_area
                        end
                        if !isnothing(area_val_discrete) && 1 ≤ period_ix ≤ num_periods
                                running_area_discrete = max(running_area_discrete, area_val_discrete)
                                daily_area_discrete[period_ix] = running_area_discrete
                        end
                end

                for arc_idx in fire_arcs_used[g]
                        arc_row = fire_models[g].long_arcs[arc_idx, :]
                        time_from = arc_row[FM.TIME_FROM]
                        time_to = arc_row[FM.TIME_TO]
                        crew_count = arc_row[FM.CREWS_PRESENT]
                        state_to = arc_row[FM.STATE_TO]

                        if 1 ≤ time_from ≤ num_periods && time_from ≥ current_day
                                daily_crews[time_from] = crew_count
                        end

                        period_ix = time_to - 1
                        packed_code = get(packed_lookup, state_to, nothing)
                        area_val = nothing
                        if packed_code !== nothing
                                area_val = decode_packed_state_area(packed_code, bins_for_decoding)
                        elseif haskey(state_area_map_sim, state_to)
                                area_val = state_area_map_sim[state_to]
                        end
                        area_val_discrete = if haskey(state_area_map_discrete, state_to)
                                state_area_map_discrete[state_to]
                            elseif !isnothing(discretization_bins_used) && state_to ≥ 1 && state_to ≤ length(discretization_bins_used)
                                discretization_bins_used[state_to]
                            else
                                nothing
                            end
                        if !isnothing(area_val) && 1 ≤ period_ix ≤ num_periods
                                running_area = max(running_area, area_val)
                                daily_area[period_ix] = running_area
                        end
                        if !isnothing(area_val_discrete) && 1 ≤ period_ix ≤ num_periods
                                running_area_discrete = max(running_area_discrete, area_val_discrete)
                                daily_area_discrete[period_ix] = running_area_discrete
                        end
                end

                prev_area = nothing
                prev_area_discrete = nothing
                for period_ix in 1:num_periods
                        if daily_area[period_ix] === nothing && prev_area !== nothing
                                daily_area[period_ix] = prev_area
                        elseif daily_area[period_ix] !== nothing
                                if prev_area === nothing || daily_area[period_ix] > prev_area
                                        prev_area = daily_area[period_ix]
                                else
                                        daily_area[period_ix] = prev_area
                                end
                        end
                        if daily_area_discrete[period_ix] === nothing && prev_area_discrete !== nothing
                                daily_area_discrete[period_ix] = prev_area_discrete
                        elseif daily_area_discrete[period_ix] !== nothing
                                if prev_area_discrete === nothing || daily_area_discrete[period_ix] > prev_area_discrete
                                        prev_area_discrete = daily_area_discrete[period_ix]
                                else
                                        daily_area_discrete[period_ix] = prev_area_discrete
                                end
                        end
        end

        stats_filename = "fire_stats_$(g)_$(t).json"
        stats_payload = Dict{String,Any}(
                        "optimizer_index" => g,
                        "fire_event_id" => get(state_entry, "fire_event_id", nothing),
                        "arc_file" => get(state_entry, "arc_file", nothing),
                        "day_index" => t,
                        "daily_crews" => daily_crews,
                        "daily_area_acres" => [daily_area[i] === nothing ? nothing : daily_area[i] for i in 1:num_periods],
                        "daily_area_acres_discrete" => [daily_area_discrete[i] === nothing ? nothing : daily_area_discrete[i] for i in 1:num_periods],
                )
                open(joinpath(output_folder, stats_filename), "w") do io
                        JSON.print(io, stats_payload)
                end
                detail_ix = min(current_day, num_periods)
                crews_today = daily_crews[detail_ix]
                area_today = daily_area[detail_ix]
                area_discrete_today = daily_area_discrete[detail_ix]
                start_day_val = init_info.start_days[g]
                base_label = replace(basename(input_folder), "fire_models_" => "")
                base_date = try
                        Date(base_label, dateformat"yyyy_mm_dd")
                catch
                        nothing
                end
                start_date_str = isnothing(base_date) ? nothing : Dates.format(base_date + Day(start_day_val), dateformat"yyyy-mm-dd")
                fire_id_val = get(state_entry, "fire_event_id", init_info.fire_ids[g])
                incident_name = get(state_entry, "incident_name", nothing)
                @info "Optimizer day detail" day=current_day fire=g fire_id=fire_id_val incident=incident_name start_day=start_day_val start_date=start_date_str crews=crews_today area=area_today area_discrete=area_discrete_today
                fire_label = if incident_name === nothing || incident_name === missing || incident_name === ""
                        "fire $(g)"
                else
                        string(incident_name)
                end
                area_str = area_today === nothing ? "n/a" : string(area_today)
                summary_entry = "fire $(fire_label) (id $(fire_id_val), start day $(start_day_val)): crews=$(crews_today), area=$(area_str)"
                push!(day_summary_entries, summary_entry)

                output_files = Dict{String,Any}(
                        "arcs" => arcs_filename,
                        "arc_costs" => arc_costs_filename,
                        "stats" => stats_filename,
                )
                state_entry["output_files"] = output_files
                state_entry["day_index"] = t
                fire_manifest_entries[g] = state_entry
        end

        @info "Day $(current_day) summary" summary=day_summary_entries
        push!(optimizer_day_rollup, (day=current_day, summary=copy(day_summary_entries)))

        for j in 1:num_crews
                crew_arcs_filename = "crew_arcs_$(j)_$(t).json"
                crew_costs_filename = "crew_arc_costs_$(j)_$(t).json"
                open(joinpath(output_folder, crew_arcs_filename), "w") do io
                        JSON.print(io, crew_arcs[j])
                end
                open(joinpath(output_folder, crew_costs_filename), "w") do io
                        JSON.print(io, crew_arc_costs[j])
                end
                push!(crew_manifest_entries, Dict{String,Any}(
                        "crew_index" => j,
                        "arcs" => crew_arcs_filename,
                        "arc_costs" => crew_costs_filename,
                ))
        end

        manifest_dict = Dict{String,Any}(
                "generated_at" => Dates.format(Dates.now(), dateformat"YYYY-mm-ddTHH:MM:SS"),
                "day_index" => t,
                "crew_step" => (:crew_step in propertynames(init_info)) ? init_info.crew_step : firefighters_per_crew,
                "initial_firefighters_per_crew" => personnel_per_crew,
                "max_crews" => num_crews,
                "horizon_periods" => num_time_periods,
                "time_period_hours" => (:time_period_hours in propertynames(init_info)) ? init_info.time_period_hours : 6,
                "travel_speed_miles_per_period" => travel_speed,
                "fire_order" => fire_manifest_entries,
                "crew_outputs" => crew_manifest_entries,
        )
        if discretization_bins_used !== nothing
                manifest_dict["discretization_bins"] = discretization_bins_used
        else
                manifest_dict["discretization_bins"] = nothing
        end

        manifest_filename = joinpath(output_folder, "arc_manifest_$(t).json")
        open(manifest_filename, "w") do io
        JSON.print(io, manifest_dict)
        end
end

@info "Optimizer horizon summary" summary=[
        begin
                day_label = entry.day
                entries = entry.summary
                formatted = isempty(entries) ? "none" : join(entries, " \\ ")
                "Day $(day_label): $(formatted)"
        end
        for entry in optimizer_day_rollup
]

end # let prev_dual_warm_start

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
