include("BranchAndPrice.jl") # load the optimizer implementation so this script can initialize and run it

# EmpiricalMain.jl is the CLI driver for running end-to-end wildfire suppression
# experiments on historical data.  The script parses command-line flags, loads
# preprocessed datasets, seeds the branch-and-price solver with intuitive plans,
# iterates a rolling-horizon optimization, and exports a rich set of artifacts
# (JSON manifests, CSV rollups, per-day logs).  The goal of the extra comments
# below is to help new contributors trace how data flows from disk into the
# optimizer and back out into reports.

# Core packages that drive the command-line workflow and optimization.
# IterTools exports a groupby helper that clashes with DataFrames, so we explicitly
# import the DataFrames version below.
using JuMP, Gurobi, JSON, Profile, ArgParse, Logging, IterTools, CSV, DataFrames, Dates, DelimitedFiles, Printf
import DataFrames: groupby
import Logging: min_enabled_level, shouldlog, handle_message # grant direct access to these logging hooks

# DualLogger mirrors all log events to both console and file outputs.
struct DualLogger <: AbstractLogger
        loggers::NTuple{2,AbstractLogger}
end

# connect the wrapper to the underlying logger methods
min_enabled_level(l::DualLogger) = min(min_enabled_level(l.loggers[1]), min_enabled_level(l.loggers[2]))
shouldlog(l::DualLogger, level, _module, group, id) =
        shouldlog(l.loggers[1], level, _module, group, id) ||
        shouldlog(l.loggers[2], level, _module, group, id)
handle_message(l::DualLogger, level, message, _module, group, id, file, line; kwargs...) = begin
        handle_message(l.loggers[1], level, message, _module, group, id, file, line; kwargs...)
        handle_message(l.loggers[2], level, message, _module, group, id, file, line; kwargs...)
end

const GRB_ENV = Gurobi.Env()

# Canonical GACC name mapping so command-line abbreviations are standardized.
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
     # Keep a unique list of region names for convenience when users specify "all".

const JULIA_EXT_STATE_CODE = 1

# Utility routines used by the fire-model preprocessing logic -----------------

# Build default area discretization bins when none are supplied alongside a dataset.
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

function load_discretization_bins_file(folder::String)
        bin_path = joinpath(folder, "discretization_bins.csv")
        if isfile(bin_path)
                data = readdlm(bin_path, ',')
                return vec(Float64.(data))
        end
        return nothing
end

function decode_packed_state_area(packed_code::Int, bins::Vector{Float64})
        packed_code <= JULIA_EXT_STATE_CODE && return 0.0
        num_bins = length(bins)
        num_bins == 0 && return 0.0
        active = packed_code - (JULIA_EXT_STATE_CODE + 1) # offset accounts for Julia index base
        next_idx = active ÷ num_bins                     # decode coarse + fine bin indices
        next_idx = clamp(next_idx, 0, num_bins - 1)      # guard against malformed codes
        return bins[next_idx + 1]                        # translate back to an area value
end

# Parse command-line GACC abbreviations into canonical names.
function parse_gaccs(str::String)
        # Accept CLI values like "gb,nw" and normalize them to canonical names
        # so the rest of the code can rely on a single naming convention.
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

# Parse a mapping such as "GB:123,456;NW:789" into a Dict keyed by canonical GACC name.
function parse_fires_by_gacc(str::String)
        # Turn a compact "GB:1,2;NW:3" specification into a Dict that mirrors
        # the CSV schema (keys are canonical names, values are fire IDs).
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

# Count how many distinct fires remain after the command-line GACC / ID filters are applied.
function count_selected_fires(
        fire_gaccs::Vector{String},
        fires_by_gacc::Dict{String,Vector{Int64}},
        input_folder::String,
)
        # This helper mirrors the eventual filtering logic so we can log how
        # many fires will be in play before allocating large data structures.
        selected_fires = CSV.read(joinpath(input_folder, "selected_fires.csv"), DataFrame)
        selected_fires[!, "GACC"] = normalize_gacc.(selected_fires[!, "GACC"]) # normalize to canonical casing
        fire_gaccs = normalize_gacc.(fire_gaccs) # make the filter list consistent too
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
        start_key = "start_day_of_sim" # preferred column for the fire start within the planning window
        alt_keys = ("sim_start_day_dsfr", "day_since_first_report", "start_day") # fallbacks seen in historical exports
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
                selected_fires[!, :start_day_of_sim] = copy(selected_fires[!, start_col]) # ensure downstream uses a consistent symbol
        end
        selected_fires[!, :GACC] = normalize_gacc.(selected_fires[!, :GACC])

        if !isempty(fires_by_gacc)
                normalized = Dict{String,Set{Int64}}()
                for (gacc, fire_list) in fires_by_gacc
                        normalized[normalize_gacc(gacc)] = Set(Int64.(fire_list)) # store unique fire ids per normalized GACC
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
                return subset # no day-0 fires match the requested filters
        end

        for subdf in groupby(filtered, :GACC)
                gacc = subdf[1, :GACC]
                subset[gacc] = collect(unique(Int64.(subdf[!, :FIRE_EVENT_ID])))
        end

        return subset
end

function build_fire_arc_lookup(fire_model::TimeSpaceNetwork)
        lookup = Dict{NTuple{5,Int64},Int64}()
        raw_from = fire_model.raw_state_from
        raw_to = fire_model.raw_state_to
        for ix in 1:size(fire_model.long_arcs, 1)
                from_raw = raw_from === nothing ? fire_model.long_arcs[ix, FM.STATE_FROM] : raw_from[ix]
                to_raw = raw_to === nothing ? fire_model.long_arcs[ix, FM.STATE_TO] : raw_to[ix]
                time_from = fire_model.long_arcs[ix, FM.TIME_FROM]
                time_to = fire_model.long_arcs[ix, FM.TIME_TO]
                crews = fire_model.long_arcs[ix, FM.CREWS_PRESENT]
                lookup[(from_raw, to_raw, time_from, time_to, crews)] = ix
        end
        return lookup
end

function build_crew_arc_lookup(crew_model::TimeSpaceNetwork)
        lookup = Dict{NTuple{8,Int64},Int64}()
        for ix in 1:size(crew_model.long_arcs, 1)
                arc = crew_model.long_arcs[ix, :]
                key = (
                        arc[CM.FROM_TYPE],
                        arc[CM.LOC_FROM],
                        arc[CM.TO_TYPE],
                        arc[CM.LOC_TO],
                        arc[CM.TIME_FROM],
                        arc[CM.TIME_TO],
                        arc[CM.REST_FROM],
                        arc[CM.REST_TO],
                )
                lookup[key] = ix
        end
        return lookup
end

function final_snapshot_cost(fire_model::TimeSpaceNetwork, arcs_used::Vector{Int}, num_time_periods::Int)
        best_arc = 0
        best_t = -1
        for a_ix in arcs_used
                to_t = fire_model.long_arcs[a_ix, FM.TIME_TO]
                c = fire_model.arc_costs[a_ix]
                if (to_t <= num_time_periods + 1) && (c > 1e-12) && (to_t - 1) > best_t
                        best_arc = a_ix
                        best_t = to_t - 1
                end
        end
        if best_arc != 0
                return fire_model.arc_costs[best_arc]
        end
        return sum(fire_model.arc_costs[arcs_used])
end

function seed_fire_plans_from_arc!(
        seed_dir::String,
        init_info,
        fire_plans::FirePlanData,
        fire_models::Vector{TimeSpaceNetwork},
        num_time_periods::Int,
        firefighters_per_crew::Int,
)
        fire_csv = joinpath(seed_dir, "selected_fire_arcs.csv")
        isfile(fire_csv) || return 0
        df = CSV.read(fire_csv, DataFrame)
        if isempty(df)
                        return 0
        end
        name_lookup = Dict(lowercase(String(col)) => col for col in names(df))
        personnel_col = get(name_lookup, "personnel", nothing)
        if personnel_col === nothing
                @warn "selected_fire_arcs.csv missing personnel column" fire_csv
                return 0
        end
        idx_col = :__fire_ix
        if haskey(name_lookup, "fire_event_id")
                fire_symbol = name_lookup["fire_event_id"]
                fire_lookup = Dict(string(init_info.fire_ids[ix]) => ix for ix in 1:length(init_info.fire_ids))
                df[!, idx_col] = [get(fire_lookup, string(fid), 0) for fid in df[!, fire_symbol]]
        elseif haskey(name_lookup, "fire")
                df[!, idx_col] = [Int(round(val)) for val in df[!, name_lookup["fire"]]]
        else
                @warn "selected_fire_arcs.csv missing fire identifier columns (fire_event_id or fire)" fire_csv
                return 0
        end
        valid_mask = (df[!, idx_col] .>= 1) .& (df[!, idx_col] .<= length(init_info.fire_ids))
        if !all(valid_mask)
                @warn "Skipping fire arcs that do not map to optimizer indices" discarded = sum(.!valid_mask)
                df = df[valid_mask, :]
        end
        crew_step = (:crew_step in propertynames(init_info)) ? init_info.crew_step : firefighters_per_crew
        total_seeded = 0
        grouped = groupby(df, idx_col)
        for group in grouped
                fire_ix = first(group[!, idx_col])
                sub = sort(DataFrame(group), [:time_from, :arc_index])
                lookup = build_fire_arc_lookup(fire_models[fire_ix])
                arcs_used = Int[]
                missing_arc = false
                for row in eachrow(sub)
                        raw_from = Int(round(row.state_from_raw))
                        raw_to = Int(round(row.state_to_raw))
                        time_from = Int(round(row.time_from))
                        time_to = Int(round(row.time_to))
                        crew_count = Int(round(row[personnel_col] / crew_step))
                        key = (raw_from, raw_to, time_from, time_to, crew_count)
                        arc_ix = get(lookup, key, 0)
                        if arc_ix == 0
                                missing_arc = true
                                fire_label = (:fire_event_id in keys(name_lookup)) ? group[1, name_lookup["fire_event_id"]] : init_info.fire_ids[fire_ix]
                                @warn "Could not map network-flow fire arc to BPC arc" fire_id=fire_label key=key
                                break
                        end
                        push!(arcs_used, arc_ix)
                end
                if missing_arc || isempty(arcs_used)
                        continue
                end
                crew_demands = zeros(Int, num_time_periods)
                for arc_ix in arcs_used
                        tf = fire_models[fire_ix].long_arcs[arc_ix, FM.TIME_FROM]
                        if 1 <= tf <= num_time_periods
                                crew_demands[tf] = fire_models[fire_ix].long_arcs[arc_ix, FM.CREWS_PRESENT]
                        end
                end
                cost = final_snapshot_cost(fire_models[fire_ix], arcs_used, num_time_periods)
                add_column_to_plan_data!(fire_plans, fire_ix, cost, crew_demands, arcs_used)
                total_seeded += 1
                @debug "Seeded fire plan from network_flow_direct" fire=fire_ix crew_demands=crew_demands cost=cost arcs=length(arcs_used)
        end
        return total_seeded
end

function seed_crew_routes_from_arc!(
        seed_dir::String,
        crew_routes::CrewRouteData,
        crew_models::Vector{TimeSpaceNetwork},
        num_fires::Int,
        num_time_periods::Int,
)
        num_crews = length(crew_models)
        seeded = 0
        for crew_ix in 1:num_crews
                base = @sprintf("crew_%03d_selected_arcs_post_processed.csv", crew_ix)
                path = joinpath(seed_dir, base)
                if !isfile(path)
                        alt = @sprintf("crew_%03d_selected_arcs.csv", crew_ix)
                        path = joinpath(seed_dir, alt)
                        isfile(path) || continue
                end
                df = CSV.read(path, DataFrame)
                has_value_col = :value in names(df)
                rows = has_value_col ?
                        [row for row in eachrow(df) if row[:value] > 1e-6] :
                        collect(eachrow(df))
                if isempty(rows)
                        continue
                end
                lookup = build_crew_arc_lookup(crew_models[crew_ix])
                arcs_used = Int[]
                missing_arc = false
                for row in rows
                        key = (
                                Int(round(row[:from_type])),
                                Int(round(row[:loc_from])),
                                Int(round(row[:to_type])),
                                Int(round(row[:loc_to])),
                                Int(round(row[:time_from])),
                                Int(round(row[:time_to])),
                                Int(round(row[:rest_from])),
                                Int(round(row[:rest_to])),
                        )
                        arc_ix = get(lookup, key, 0)
                        if arc_ix == 0
                                missing_arc = true
                                @warn "Could not map network-flow crew arc to BPC arc" crew_ix key
                                break
                        end
                        push!(arcs_used, arc_ix)
                end
                if missing_arc || isempty(arcs_used)
                        continue
                end
                fires_fought = get_fires_fought(
                        crew_models[crew_ix].wide_arcs,
                        arcs_used,
                        (num_fires, num_time_periods),
                )
                cost = sum(crew_models[crew_ix].arc_costs[arcs_used])
                add_column_to_route_data!(crew_routes, crew_ix, cost, fires_fought, arcs_used)
                seeded += 1
                @debug "Seeded crew route from network_flow_direct" crew=crew_ix fires_fought=fires_fought cost=cost arcs=length(arcs_used)
        end
        return seeded
end

function seed_from_arc_solution!(
        seed_dir::String,
        init_info,
        crew_routes::CrewRouteData,
        fire_plans::FirePlanData,
        crew_models::Vector{TimeSpaceNetwork},
        fire_models::Vector{TimeSpaceNetwork},
        num_time_periods::Int,
        firefighters_per_crew::Int,
)
        path = strip(seed_dir)
        if isempty(path)
                return
        end
        arc_dir = abspath(path)
        if !isdir(arc_dir)
                @warn "Seed arc directory does not exist" arc_dir
                return
        end
        fire_seeded = seed_fire_plans_from_arc!(
                arc_dir,
                init_info,
                fire_plans,
                fire_models,
                num_time_periods,
                firefighters_per_crew,
        )
        crew_seeded = seed_crew_routes_from_arc!(
                arc_dir,
                crew_routes,
                crew_models,
                length(fire_models),
                num_time_periods,
        )
        @info "Seeded columns from network_flow_direct output" arc_dir fire_plans=fire_seeded crew_routes=crew_seeded
end

function resolve_input_folder(folder::String)
        # try the user-specified path directly
        # if the provided path already points to a folder with selected_fires.csv, use it
        if isfile(joinpath(folder, "selected_fires.csv"))
                return folder
        end
        # check for an arc_arrays subdirectory (common layout for prepared datasets)
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
        # ArgParse describes all supported CLI flags for the empirical workflow.
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
                "--perfect-info"
                help = "solve the entire horizon once with full knowledge of future fires (no rolling re-solve)"
                action = :store_true
                "--crew-costs"
                help = "crew cost mode: 'on' (default), 'off', or 'rest' to keep only rest penalties"
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
                "--baseline-label"
                help = "Label to store in optimizer arc CSV outputs for downstream comparisons"
                default = "Optimizer"
                "--time-limit"
                help = "Time limit in seconds for the branch-and-price algorithm"
                arg_type = Float64
                default = 1800.0
                "--rest-periods"
                help = "Number of consecutive periods a crew must rest once they return to base"
                arg_type = Int
                default = 3
                "--seed-arc-output-dir"
                help = "Optional path to a network_flow_direct output folder whose solution will be seeded into the master problem"
                default = ""
                "--input-folder"
                help = "Directory containing input files (full path or dataset under data/empirical_fire_models)"
                default = "raw"
                "--output-folder"
                help = "Directory to store output files"
                default = "data/output"
                "--scalability-benchmark"
                help = "also write scalability_output.json in the format expected by process_outputs.py scalability table"
                action = :store_true
        end
        return parse_args(arg_parse_settings)
end


args = get_command_line_args() # parse and validate all CLI input
crew_gaccs = parse_gaccs(args["crew-gaccs"]) # map crew GACC abbreviations to canonical names
fire_gaccs = isempty(args["fire-gaccs"]) ? crew_gaccs : parse_gaccs(args["fire-gaccs"]) # default fires to same set
seed_fastest_extinguish = args["seed-fastest-extinguish"]
seed_frontloaded = args["seed-frontloaded"]
seed_max_daily = args["seed-max-daily"]
seed_crew_routes = args["seed-crew-routes"]
final_snapshot_only = args["final-snapshot-only"]
day_one_only = args["day-1-only"] # new flag that activates the day-0 fire filter
perfect_info_mode = args["perfect-info"]
crew_costs_mode = lowercase(String(args["crew-costs"]))
rest_penalties_only = crew_costs_mode in ("rest", "rest-only", "rest_only")
zero_crew_costs = rest_penalties_only || (crew_costs_mode in ("off","0","false","no")) # treat rest-only like "off" for travel/fight costs
firefighters_per_crew = args["firefighters-per-crew"]
personnel_per_crew = args["personnel-per-crew"]
rest_periods = Int(args["rest-periods"])
fires_by_gacc = parse_fires_by_gacc(args["fires"]) # optional explicit fire list
time_limit = args["time-limit"]
input_folder = resolve_input_folder(args["input-folder"]) # locate arc_arrays directory automatically if needed
output_folder = args["output-folder"]
baseline_label = String(args["baseline-label"])
seed_arc_output_dir = String(strip(String(args["seed-arc-output-dir"])))

if day_one_only
        # Build a reduced GACC -> fire list containing only incidents active on day 0.
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

mkpath(output_folder) # ensure the output directory exists before writing logs/artifacts

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
@info "Perfect info mode" perfect_info_mode
@info "Crew costs mode" crew_costs_mode rest_penalties_only=rest_penalties_only rest_periods=rest_periods
@info "Arc CSV baseline label" baseline_label

num_fires = count_selected_fires(fire_gaccs, fires_by_gacc, input_folder)
num_crews = 0

num_time_periods = 14
travel_speed = 40.0 * 6.0
GC.gc()

# Load all arc arrays, resource pools, and helper lookups from disk.  The
# initialize_* routine returns both the JuMP column data structures and the raw
# time-space networks that feed subproblems.
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
        enforce_rest_penalties = rest_penalties_only,
        rest_periods = rest_periods,
)

num_crews = length(crew_models)
num_fires = length(fire_models) # dataset-driven counts might differ from initial guesses

if !isempty(seed_arc_output_dir)
        try
                seed_from_arc_solution!(
                        seed_arc_output_dir,
                        init_info,
                        crew_routes,
                        fire_plans,
                        crew_models,
                        fire_models,
                        num_time_periods,
                        firefighters_per_crew,
                )
        catch err
                @warn "Failed to seed columns from arc solution" seed_arc_output_dir err
        end
end

crew_names = hasproperty(init_info, :crew_names) ? init_info.crew_names : nothing
if crew_names !== nothing && !isempty(output_folder)
        try
                mapping = DataFrame(
                        crew_index = collect(1:length(crew_names)),
                        crew_name = crew_names,
                )
                CSV.write(joinpath(output_folder, "crew_index_mapping.csv"), mapping)
        catch
        end
end

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
# for j in 1:num_crews
# 	no_fire_anticipation!(crew_models[j], [fsp.start_time_period for fsp in fire_models]) # ensure crew subproblems respect fire activation timing
# end

# Track committed arcs (history) so past-day decisions are preserved across re-solves
committed_fire_arcs = [Set{Int64}() for _ in 1:num_fires]
committed_crew_arcs = [Set{Int64}() for _ in 1:num_crews]

let prev_dual_warm_start = nothing
optimizer_day_rollup = Any[]
optimizer_day_arc_records = Vector{Dict{Symbol,Any}}()
final_fire_arcs = nothing
final_fire_arc_costs = nothing
final_crew_arcs = nothing
final_crew_arc_costs = nothing

# Main rolling-horizon loop: each iteration finalizes decisions for the current day,
# seeds additional columns, and re-optimizes the branch-and-price master problem.
# Perfect-info mode only executes the first iteration (t = 0) with omniscient foresight.
rolling_loop_end = perfect_info_mode ? 0 : num_time_periods
for t in 0:rolling_loop_end

    global crew_routes, fire_plans, crew_models, fire_models, cut_data

    # Start each day with fresh column containers (the problem structure changes
    # as fires start/finish), but keep the underlying TSN models so we can prune
    # arcs based on committed history.
    crew_routes = CrewRouteData(Int(floor(6 * 1e6 / num_crews)), num_fires, num_crews, num_time_periods)
    fire_plans = FirePlanData(Int(floor(6 * 1e6  / num_crews)), num_fires, num_time_periods)
	cut_data = CutData(num_crews, num_fires, num_time_periods)

    # Seed dummy plan/route columns so the restricted master is always feasible.
    # These enormous-cost columns act as big-M slack so Gurobi has something to
    # work with before legitimate routes/plans are priced.
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
    fires_active = perfect_info_mode ?
        collect(1:num_fires) :
        [g for g in 1:num_fires if isnothing(fire_start_periods[g]) || fire_start_periods[g] <= current_day]
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
        # (Some seeding strategies use a single representative arc from the terminal day.)
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
        # Build a dynamic program that finds the earliest path to the extinguish
        # state (state 1).  This encourages the master problem to explore plans
        # that act aggressively as soon as the fire is active.
        if seed_fastest_extinguish
            for g in fires_active
                fm = fire_models[g]
                states, times = size(fm.state_in_arcs)
                # Earliest reachability to extinguish state (assume state 1)
                prev_arc = fill(0, states, times)
                reachable = falses(states, times)
                ext_state_id = 1
                found_t = nothing
                # Forward pass: mark which states can be reached by time tt given
                # the currently known arcs (ignoring reduced costs entirely).
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
                # (Reverse pointers reconstruct the minimum-time plan.)
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
                # Require at least one post-start arc so the plan extends into the active window.
                start_day = fire_models[g].start_time_period
                if !isnothing(start_day)
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
                # Require at least one post-start arc so the plan extends into the active window.
                start_day = fire_models[g].start_time_period
                if !isnothing(start_day)
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
                # Require at least one post-start arc so the plan extends into the active window.
                start_day = fire_models[g].start_time_period
                if !isnothing(start_day)
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
        # by taking a single high-crew arc at that day and then greedily continuing.
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

        # Reuse duals from the previous day whenever the matrix dimensions match;
        # this gives the next branch-and-price solve a head start.
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

        # Solve the full branch-and-price model for the current rolling-horizon
        # day.  Most arguments mirror CLI options so experiments can be reproduced.
        result = branch_and_price(num_fires,
                num_crews,
                num_time_periods,
                current_time = t,
                from_empirical = true,
                gaccs = crew_gaccs,
                fire_gaccs = fire_gaccs,
                travel_speed = travel_speed,
                rest_periods = rest_periods,
                firefighters_per_crew = firefighters_per_crew,
                initial_firefighters_per_crew = personnel_per_crew,
                fires_by_gacc = fires_by_gacc,
                input_folder = input_folder,
                crew_routes = crew_routes,
                fire_plans = fire_plans,
                crew_models = crew_models,
                fire_models = fire_models,
                cut_data = cut_data,
                algo_tracking = true,
                total_time_limit = time_limit,
                output_folder = output_folder,
                dual_warm_start = warm_start_to_use,
                final_snapshot_only = final_snapshot_only,
                include_future_fires = perfect_info_mode,
                )
        # unpack the full return tuple (node stats, incumbent bounds, selected arcs, warm starts, …)
        explored_nodes, ubs, lbs, columns, heuristic_times, times, time_1, root_node_ip_sol, root_node_ip_sol_time, fire_arcs_used, crew_arcs_used, root_dual_warm_start, final_ub, final_lb, total_solve_time, algo_summary = result
        @info "BPC explored_nodes" explored_nodes = explored_nodes
        @info "BPC upper bounds" ubs = ubs
        @info "BPC lower bounds" lbs = lbs
        @info "BPC column counts" columns = columns
        @info "BPC heuristic times" heuristic_times = heuristic_times
        @info "BPC cumulative times" times = times
        @info "BPC time_1 snapshot" time_1 = time_1
        @info "BPC root node IP objective" root_node_ip_sol = root_node_ip_sol root_node_ip_sol_time = root_node_ip_sol_time
        @info "BPC final bounds" final_ub = final_ub final_lb = final_lb
        @info "BPC total solve time" total_solve_time = total_solve_time
        @debug "final arcs used" fire_arcs_used, crew_arcs_used

        if fire_arcs_used === nothing || crew_arcs_used === nothing
                @warn "No arc information returned from branch_and_price; stopping early"
                break
        end

        prev_dual_warm_start = root_dual_warm_start

        if !perfect_info_mode
                # Commit decisions up to current day (t+1): preserve those arcs in future iterations
                # so the rolling horizon respects previously executed suppression actions.
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
                        # by pruning away trajectories we never plan to use again.  This keeps
                        # subsequent DP solves small even as the horizon marches forward.
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
        end

	# After the day is solved we snapshot the exact arcs we used so downstream
	# visualization/logging code can refer to them without re-solving DP problems.
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
        if t == rolling_loop_end
                final_fire_arcs = deepcopy(fire_arcs)
                final_fire_arc_costs = deepcopy(fire_arc_costs)
                final_crew_arcs = deepcopy(crew_arcs)
                final_crew_arc_costs = deepcopy(crew_arc_costs)
        end

        discretization_bins_used = (
                (init_info !== nothing) && (:discretization_bins in propertynames(init_info))
        ) ? init_info.discretization_bins : nothing
        if discretization_bins_used === nothing
                discretization_bins_used = load_discretization_bins_file(input_folder)
        end
        bins_for_decoding = discretization_bins_used === nothing ? default_discretization_bins() : collect(Float64.(discretization_bins_used))
        # Fire manifests capture per-fire metadata (IDs, state info, etc.) and
        # get enriched each day with references to the newly generated files.
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
                if !(state_meta isa AbstractVector)
                        state_meta = Any[]
                        state_entry["state_metadata"] = state_meta
                end
                if state_meta isa AbstractVector
                        # Normalize metadata into quick lookup tables so we can
                        # map optimizer state indices back to geospatial/area info.
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
                        # Fill missing discrete areas by decoding packed codes so
                        # downstream logs never see a zero/unknown value.
                        for (sid, sm) in state_lookup
                                if !haskey(state_area_map_discrete, sid)
                                        packed_code = get(packed_lookup, sid, 0)
                                        decoded_area = decode_packed_state_area(packed_code, bins_for_decoding)
                                        if decoded_area > 0
                                                state_area_map_discrete[sid] = decoded_area
                                        end
                                end
                                if !haskey(state_area_map_sim, sid)
                                        packed_code = get(packed_lookup, sid, 0)
                                        decoded_area = decode_packed_state_area(packed_code, bins_for_decoding)
                                        if decoded_area > 0
                                                state_area_map_sim[sid] = decoded_area
                                        end
                                end
                        end
                end
                begin
                        existing_packed_codes = Set(values(packed_lookup))
                        raw_state_codes = isnothing(fire_models[g].raw_state_to) ? Int[] : unique(Int.(fire_models[g].raw_state_to))
                        next_state_id = isempty(state_lookup) ? 1 : (maximum(keys(state_lookup)) + 1)
                        for code in raw_state_codes
                                code <= 0 && continue
                                if code in existing_packed_codes
                                        continue
                                end
                                area_val = decode_packed_state_area(code, bins_for_decoding)
                                area_val <= 0 && continue
                                entry = Dict{String,Any}(
                                        "state_id" => next_state_id,
                                        "packed_code" => code,
                                        "area_acres_sim" => area_val,
                                        "area_acres_discrete" => area_val,
                                )
                                push!(state_meta, entry)
                                state_lookup[next_state_id] = entry
                                state_area_map_sim[next_state_id] = area_val
                                state_area_map_discrete[next_state_id] = area_val
                                packed_lookup[next_state_id] = code
                                push!(existing_packed_codes, code)
                                next_state_id += 1
                        end
                        state_entry["state_metadata"] = state_meta
                end

                fire_arcs_export = fire_arcs[g]
                if !isempty(state_lookup)
                        # When metadata contains "packed state codes" we translate the
                        # raw optimizer node IDs into those codes so GIS tools can
                        # align areas-of-origin with the original dataset.
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

                # Convert optimizer states into interpretable acreage values.
                state_area_map = state_area_map_sim
                state_area_map_discrete = isempty(state_area_map_discrete) ? Dict{Int,Float64}() : state_area_map_discrete

                num_periods = num_time_periods
                # Preallocate daily time-series arrays for crews and areas; they
                # will be filled using both historical commitments and new plans.
                daily_crews = fill(0, num_periods)
                daily_area = Vector{Union{Nothing, Float64}}(undef, num_periods)
                daily_area_discrete = Vector{Union{Nothing, Float64}}(undef, num_periods)
                running_area = 0.0
                running_area_discrete = 0.0
                period_area_exact = Dict{Int, Tuple{Int, Float64}}()
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
                        area_for_summary = get(state_area_map_sim, state_to, nothing)
                        if area_for_summary === nothing
                                area_for_summary = area_val
                        end
                        if !isnothing(area_val) && 1 ≤ period_ix ≤ num_periods
                                running_area = max(running_area, area_val)
                                daily_area[period_ix] = running_area
                        end
                        if !isnothing(area_val_discrete) && 1 ≤ period_ix ≤ num_periods
                                running_area_discrete = max(running_area_discrete, area_val_discrete)
                                daily_area_discrete[period_ix] = running_area_discrete
                        end
                        if !isnothing(area_for_summary) && 0 ≤ period_ix ≤ num_periods
                                existing = get(period_area_exact, period_ix, nothing)
                                if existing === nothing || time_from ≥ existing[1]
                                        period_area_exact[period_ix] = (time_from, area_for_summary)
                                end
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
                        area_for_summary = get(state_area_map_sim, state_to, nothing)
                        if area_for_summary === nothing
                                area_for_summary = area_val
                        end
                        if !isnothing(area_val) && 1 ≤ period_ix ≤ num_periods
                                running_area = max(running_area, area_val)
                                daily_area[period_ix] = running_area
                        end
                        if !isnothing(area_val_discrete) && 1 ≤ period_ix ≤ num_periods
                                running_area_discrete = max(running_area_discrete, area_val_discrete)
                                daily_area_discrete[period_ix] = running_area_discrete
                        end
                        if !isnothing(area_for_summary) && 0 ≤ period_ix ≤ num_periods
                                existing = get(period_area_exact, period_ix, nothing)
                                if existing === nothing || time_from ≥ existing[1]
                                        period_area_exact[period_ix] = (time_from, area_for_summary)
                                end
                        end
                end

                # Forward-fill any missing area data to keep plots monotone.
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
                        if daily_area_discrete[period_ix] === nothing && daily_area[period_ix] !== nothing
                                daily_area_discrete[period_ix] = daily_area[period_ix]
                        end
                end

        # Persist a per-fire stats JSON summarizing the day's crew assignments
        # and resulting area trajectories.  Downstream dashboards consume this.
        stats_filename = "fire_stats_$(g)_$(t).json"
        stats_payload = Dict{String,Any}(
                        "optimizer_index" => g,
                        "fire_event_id" => get(state_entry, "fire_event_id", nothing),
                        "arc_file" => get(state_entry, "arc_file", nothing),
                        "day_index" => t,
                        "daily_crews" => daily_crews,
                        "daily_area_acres" => [daily_area_discrete[i] === nothing ? nothing : daily_area_discrete[i] for i in 1:num_periods],
                        "daily_area_acres_discrete" => [daily_area_discrete[i] === nothing ? nothing : daily_area_discrete[i] for i in 1:num_periods],
                )
                open(joinpath(output_folder, stats_filename), "w") do io
                        JSON.print(io, stats_payload)
                end
                # Compose a short narrative for the day's decisions to aid debugging
                if perfect_info_mode
                        detail_ix = num_periods
                else
                        detail_ix = min(current_day, num_periods)
                end
                crews_today = daily_crews[detail_ix]
                area_today = daily_area[detail_ix]
                area_discrete_today = daily_area_discrete[detail_ix]
                if area_discrete_today === nothing && area_today !== nothing
                        area_discrete_today = Float64(area_today)
                end
                area_today_exact = begin
                        candidate_ix = detail_ix
                        found = nothing
                        while candidate_ix ≥ 0
                                entry = get(period_area_exact, candidate_ix, nothing)
                                if entry !== nothing
                                        found = entry[2]
                                        break
                                end
                                candidate_ix -= 1
                        end
                        found
                end
                if area_today_exact !== nothing
                        area_today = area_today_exact
                        if area_discrete_today === nothing
                                area_discrete_today = Float64(area_today_exact)
                        end
                end
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

                fire_id_field = if fire_id_val === nothing || fire_id_val === missing
                        missing
                else
                        Int(fire_id_val)
                end
                start_day_field = if start_day_val === nothing || start_day_val === missing
                        missing
                else
                        Int(start_day_val)
                end
                crews_field = if crews_today === nothing || crews_today === missing
                        missing
                else
                        Int(crews_today)
                end
                area_field = area_today === nothing ? missing : area_today

                area_str = area_today === nothing ? "n/a" : string(area_today)
                summary_entry = "fire $(fire_label) (id $(fire_id_val), start day $(start_day_val)): crews=$(crews_today), area=$(area_str)"
                push!(day_summary_entries, summary_entry)
                if current_day <= num_time_periods
                        push!(optimizer_day_arc_records, Dict(
                                :baseline => baseline_label,
                                :day_index => t,
                                :day_number => current_day,
                                :fire_id => fire_id_field,
                                :start_day => start_day_field,
                                :current_area => area_field,
                                :assigned_crews => crews_field,
                        ))
                end

                output_files = Dict{String,Any}(
                        "arcs" => arcs_filename,
                        "arc_costs" => arc_costs_filename,
                        "stats" => stats_filename,
                )
                state_entry["output_files"] = output_files
                state_entry["day_index"] = t
                fire_manifest_entries[g] = state_entry
        end

        # Keep a rolling log of per-day textual summaries for CLI output.
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

        # The manifest ties together every JSON/CSV asset generated for this day
        # so downstream tools can locate the correct files without recomputing.
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

        objective_val = isfinite(final_ub) ? final_ub : nothing
        bound_val = isfinite(final_lb) ? final_lb : nothing
        gap_val = nothing
        if objective_val !== nothing && bound_val !== nothing
                denom = max(1.0, abs(objective_val))
                gap_val = abs(objective_val - bound_val) / denom
        end
        heuristic_best = get(algo_summary, :heuristic_best_objective, nothing)
        if heuristic_best isa Real && !isfinite(heuristic_best)
                heuristic_best = nothing
        end
        summary_payload = Dict{String,Any}(
                "objective" => objective_val,
                "bound" => bound_val,
                "gap" => gap_val,
                "solve_seconds" => total_solve_time,
                "day_index" => t,
                "explored_nodes_series" => explored_nodes,
                "upper_bounds_series" => ubs,
                "lower_bounds_series" => lbs,
                "column_counts_series" => columns,
                "heuristic_times_series" => heuristic_times,
                "cumulative_times_series" => times,
                "time_1_snapshot" => time_1,
                "root_node_ip_objective" => root_node_ip_sol,
                "root_node_ip_time" => root_node_ip_sol_time,
                "final_upper_bound" => final_ub,
                "final_lower_bound" => final_lb,
                "algo_crew_subproblem_time" => get(algo_summary, :crew_subproblem_time, nothing),
                "algo_fire_subproblem_time" => get(algo_summary, :fire_subproblem_time, nothing),
                "algo_master_problem_time" => get(algo_summary, :master_problem_time, nothing),
                "algo_crew_subproblem_solves" => get(algo_summary, :crew_subproblem_solves, nothing),
                "algo_fire_subproblem_solves" => get(algo_summary, :fire_subproblem_solves, nothing),
                "algo_crew_columns_added" => get(algo_summary, :crew_columns_added, nothing),
                "algo_fire_columns_added" => get(algo_summary, :fire_columns_added, nothing),
                "algo_cuts_added_total" => get(algo_summary, :cuts_added, nothing),
                "algo_cut_separation_time" => get(algo_summary, :cut_separation_time, nothing),
                "algo_dcg_iterations" => get(algo_summary, :dcg_iterations, nothing),
                "algo_crew_route_pool_size" => get(algo_summary, :crew_route_pool_size, nothing),
                "algo_fire_plan_pool_size" => get(algo_summary, :fire_plan_pool_size, nothing),
                "algo_total_columns_active" => get(algo_summary, :total_columns_active, nothing),
                "algo_total_cuts_generated" => get(algo_summary, :total_cuts_generated, nothing),
                "algo_final_binding_cuts" => get(algo_summary, :final_binding_cuts, nothing),
                "algo_heuristic_runs" => get(algo_summary, :heuristic_runs, nothing),
                "algo_heuristic_total_time" => get(algo_summary, :heuristic_total_time, nothing),
                "algo_heuristic_best_objective" => heuristic_best,
                "algo_warm_start_used" => get(algo_summary, :warm_start_used, nothing),
                "algo_warm_start_rows" => get(algo_summary, :warm_start_rows, nothing),
                "algo_warm_start_cols" => get(algo_summary, :warm_start_cols, nothing),
        )
        if !isempty(explored_nodes)
                summary_payload["explored_nodes"] = explored_nodes[end]
        end
        open(joinpath(output_folder, "run_summary.json"), "w") do io
                JSON.print(io, summary_payload)
        end

        if args["scalability-benchmark"]
                scalability_payload = Dict{String,Any}(
                        "explored_nodes" => explored_nodes,
                        "upper_bounds"   => ubs,
                        "lower_bounds"   => lbs,
                        "times"          => times,
                        "heuristic_times" => heuristic_times,
                        "num_columns"    => columns,
                        "init_time"      => time_1,
                )
                open(joinpath(output_folder, "scalability_output.json"), "w") do io
                        JSON.print(io, scalability_payload, 4)
                end
        end
end

if !isempty(optimizer_day_arc_records)
        # Aggregate day-level metrics across the whole horizon for quick plotting.
        arc_df = DataFrame(optimizer_day_arc_records)
        select!(arc_df, [:baseline, :day_index, :day_number, :fire_id, :start_day, :current_area, :assigned_crews])
        CSV.write(joinpath(output_folder, "optimizer_arc_summary.csv"), arc_df)
end

if (final_fire_arcs !== nothing) && (final_fire_arc_costs !== nothing)
        fire_rows = Vector{Dict{Symbol,Any}}()
        selected_fire_arc_records = DataFrame(
                fire = Int[],
                arc_index = Int[],
                value = Float64[],
                cost = Float64[],
                state_from_raw = Int[],
                time_from = Int[],
                time_to = Int[],
                state_to_raw = Int[],
                personnel = Float64[],
        )
        final_arcs = final_fire_arcs
        final_costs = final_fire_arc_costs
        num_fires_final = length(final_arcs)
        crew_step_for_personnel = (init_info !== nothing && (:crew_step in propertynames(init_info))) ?
                init_info.crew_step : firefighters_per_crew
        fire_order_entries = (init_info !== nothing && (:fire_order in propertynames(init_info))) ? init_info.fire_order : nothing
        for g in 1:num_fires_final
                arc_matrix = final_arcs[g]
                arc_costs = final_costs[g]
                num_arcs = size(arc_matrix, 2)
                fire_event_id = (!isnothing(init_info) && (:fire_ids in propertynames(init_info))) ? init_info.fire_ids[g] : missing
                state_entry = (fire_order_entries !== nothing && g ≤ length(fire_order_entries)) ? fire_order_entries[g] : Dict{String,Any}()
                state_meta = get(state_entry, "state_metadata", nothing)
                packed_lookup = Dict{Int,Int}()
                if state_meta isa AbstractVector
                        for sm in state_meta
                                if sm isa Dict && haskey(sm, "state_id")
                                        sid = Int(sm["state_id"])
                                        raw = get(sm, "packed_code", nothing)
                                        if !(raw === nothing || raw === missing)
                                                packed_lookup[sid] = Int(raw)
                                        end
                                end
                        end
                end
                map_state = s -> get(packed_lookup, s, s)
                for k in 1:num_arcs
                        acres_val = arc_costs[k] * 1.0e4
                        crews_val = arc_matrix[FM.CREWS_PRESENT, k]
                        time_from = Int(arc_matrix[FM.TIME_FROM, k])
                        time_to = Int(arc_matrix[FM.TIME_TO, k])
                        day_idx = time_to - 1
                        state_from_raw = map_state(Int(arc_matrix[FM.STATE_FROM, k]))
                        state_to_raw = map_state(Int(arc_matrix[FM.STATE_TO, k]))
                        push!(fire_rows, Dict(
                                :fire => g,
                                :fire_event_id => fire_event_id,
                                :day => day_idx,
                                :acres => acres_val,
                                :crews => crews_val * (crew_step_for_personnel / 50.0),
                        ))
                        push!(selected_fire_arc_records, (
                                fire = g,
                                arc_index = k,
                                value = 1.0,
                                cost = arc_costs[k],
                                state_from_raw = state_from_raw,
                                time_from = time_from,
                                time_to = time_to,
                                state_to_raw = state_to_raw,
                                personnel = Float64(crews_val * crew_step_for_personnel),
                        ))
                end
        end
        fire_progression_df = DataFrame(fire_rows)
        CSV.write(joinpath(output_folder, "fire_progression_summary.csv"), fire_progression_df)
        CSV.write(joinpath(output_folder, "selected_fire_arcs.csv"), selected_fire_arc_records)

        if final_crew_arcs !== nothing
                fire_daily_counts = [Int[] for _ in 1:num_fires_final]
                crews_on_fire_counts = Int[]
                crews_in_transit_counts = Int[]
                crews_resting_counts = Int[]
                crews_at_base_counts = Int[]
                crew_status_records = [Vector{Symbol}() for _ in 1:length(final_crew_arcs)]
                crew_cost_vectors = final_crew_arc_costs === nothing ? [Float64[] for _ in final_crew_arcs] : final_crew_arc_costs

                function ensure_length!(vec::Vector{Int}, len::Int)
                        while length(vec) < len
                                push!(vec, 0)
                        end
                end

                function accumulate_range!(vec::Vector{Int}, start_day::Int, end_day::Int, amount::Int)
                        if end_day <= start_day
                                return
                        end
                        for day in max(start_day, 0):(end_day - 1)
                                idx = day + 1
                                ensure_length!(vec, idx)
                                vec[idx] += amount
                        end
                end

                function classify_arc_status(arc)::Symbol
                        ft = arc[CM.FROM_TYPE]
                        tt = arc[CM.TO_TYPE]
                        lf = arc[CM.LOC_FROM]
                        lt = arc[CM.LOC_TO]
                        rt = arc[CM.REST_TO]
                        if tt == CM.FIRE_CODE && ft == CM.FIRE_CODE && lf == lt
                                return :fire
                        elseif (ft == CM.BASE_CODE && tt == CM.FIRE_CODE) || (ft == CM.FIRE_CODE && tt == CM.BASE_CODE)
                                return :travel
                        elseif tt == CM.FIRE_CODE && ft == CM.FIRE_CODE && lf != lt
                                return :travel
                        elseif tt == CM.BASE_CODE && rt > 0
                                return :rest
                        elseif ft == CM.BASE_CODE && tt == CM.BASE_CODE && rt > 0
                                return :rest
                        elseif tt == CM.BASE_CODE
                                return :base
                        else
                                return :base
                        end
                end

                for (crew_idx, crew_matrix) in enumerate(final_crew_arcs)
                        num_arcs = size(crew_matrix, 2)
                        statuses = Vector{Symbol}(undef, max(num_arcs, 0))
                        for k in 1:num_arcs
                                arc = crew_matrix[:, k]
                                status = classify_arc_status(arc)
                                statuses[k] = status
                                start_day = Int(max(0, arc[CM.TIME_FROM]))
                                end_day = Int(max(0, arc[CM.TIME_TO]))
                                if status == :travel
                                        accumulate_range!(crews_in_transit_counts, start_day, end_day, 1)
                                elseif status == :rest
                                        accumulate_range!(crews_resting_counts, start_day, end_day, 1)
                                elseif status == :base
                                        accumulate_range!(crews_at_base_counts, start_day, end_day, 1)
                                else
                                        accumulate_range!(crews_on_fire_counts, start_day, end_day, 1)
                                        loc = Int(arc[CM.LOC_TO])
                                        if 1 <= loc <= length(fire_daily_counts)
                                                vec = fire_daily_counts[loc]
                                                for day in max(0, start_day):(end_day - 1)
                                                        idx = day + 1
                                                        ensure_length!(vec, idx)
                                                        vec[idx] += 1
                                                end
                                        end
                                end
                        end
                        crew_status_records[crew_idx] = statuses
                end

                fire_daily_records = DataFrame(
                        fire = Int[],
                        fire_event_id = String[],
                        day = Int[],
                        crews = Int[],
                )
                for (fire_idx, counts) in enumerate(fire_daily_counts)
                        fire_event_id_val = (!isnothing(init_info) && (:fire_ids in propertynames(init_info))) ? init_info.fire_ids[fire_idx] : fire_idx
                        for (day_idx, crew_val) in enumerate(counts)
                                crew_val <= 0 && continue
                                push!(fire_daily_records, (
                                        fire = fire_idx,
                                        fire_event_id = string(fire_event_id_val),
                                        day = day_idx - 1,
                                        crews = crew_val,
                                ))
                        end
                end
                CSV.write(joinpath(output_folder, "fire_daily_crews.csv"), fire_daily_records)

                max_days = max(1, num_time_periods)
                for vec in (crews_on_fire_counts, crews_in_transit_counts, crews_resting_counts, crews_at_base_counts)
                        ensure_length!(vec, max_days)
                        if length(vec) > max_days
                                resize!(vec, max_days)
                        end
                end
                status_df = DataFrame(
                        baseline = fill(baseline_label, max_days),
                        day_index = collect(0:(max_days - 1)),
                        day_number = collect(1:max_days),
                        crews_on_fire = crews_on_fire_counts,
                        crews_in_transit = crews_in_transit_counts,
                        crews_resting = crews_resting_counts,
                        crews_at_base = crews_at_base_counts,
                )
                CSV.write(joinpath(output_folder, "crew_status.csv"), status_df)

                for (crew_idx, crew_matrix) in enumerate(final_crew_arcs)
                        num_arcs = size(crew_matrix, 2)
                        num_arcs == 0 && continue
                        cost_vec = crew_cost_vectors[crew_idx]
                        if length(cost_vec) != num_arcs
                                cost_vec = collect(cost_vec)
                                while length(cost_vec) < num_arcs
                                        push!(cost_vec, 0.0)
                                end
                                if length(cost_vec) > num_arcs
                                        resize!(cost_vec, num_arcs)
                                end
                        end
                        df = DataFrame(
                                arc_index = collect(1:num_arcs),
                                value = ones(Float64, num_arcs),
                                cost = cost_vec,
                                crew_number = fill(crew_idx, num_arcs),
                                from_type = Int.(vec(crew_matrix[CM.FROM_TYPE, :])),
                                loc_from = Int.(vec(crew_matrix[CM.LOC_FROM, :])),
                                to_type = Int.(vec(crew_matrix[CM.TO_TYPE, :])),
                                loc_to = Int.(vec(crew_matrix[CM.LOC_TO, :])),
                                time_from = Int.(vec(crew_matrix[CM.TIME_FROM, :])),
                                time_to = Int.(vec(crew_matrix[CM.TIME_TO, :])),
                                rest_from = Int.(vec(crew_matrix[CM.REST_FROM, :])),
                                rest_to = Int.(vec(crew_matrix[CM.REST_TO, :])),
                        )
                        statuses = crew_status_records[crew_idx]
                        df[!, :status] = [String(statuses[k]) for k in 1:num_arcs]
                        filename = string("crew_", lpad(string(crew_idx), 3, '0'), "_selected_arcs_post_processed.csv")
                        CSV.write(joinpath(output_folder, filename), df)
                end
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
