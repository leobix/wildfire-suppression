include("Subproblems.jl")  # Load common TSN builders/constants.

# This variant mirrors the single-shot solver but wraps it inside a rolling
# horizon where future fires remain invisible until their modeled start day.
# Every day we freeze arcs that have already been executed, re-solve the full
# network-flow MIP for the remaining horizon, export per-day artifacts, advance
# the horizon one day forward, and only then allow newly ignited fires to enter
# the optimizer so there is no future knowledge.

using JuMP, Gurobi, JSON, ArgParse
using DataFrames, CSV

const GRB_ENV = Gurobi.Env()

raw_state_lookup = Vector{Dict{Int,Int}}()
crew_step = 50

# Count fires selected for a run (copied from network_flow_direct.jl).
"""
    count_selected_fires(fire_gaccs, fires_by_gacc, input_folder)

Reads `selected_fires.csv` and counts how many unique fires should appear in the
run given the GACC filters. The logic matches `network_flow_direct.jl` so both
drivers agree on dataset selection.
"""
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
            mask .|= (selected_fires[:, "GACC"] .== gacc) .&
                in.(selected_fires[:, "FIRE_EVENT_ID"], Ref(fires))
        end
        selected_fires = selected_fires[mask, :]
    else
        selected_fires = selected_fires[in.(selected_fires[:, "GACC"], Ref(fire_gaccs)), :]
    end
    return length(unique(selected_fires[:, "FIRE_EVENT_ID"]))
end

slugify(str) = replace(lowercase(strip(str)), r"[^0-9a-z]+" => "_")

function get_command_line_args()
    """
    Build/parse the CLI flags. All flags match network_flow_direct.jl plus
    `--time_limit`, which controls the per-day solver wall clock limit.
    """
    settings = ArgParseSettings()
    @add_arg_table settings begin
        "--debug"
        help = "enable verbose logging"
        action = :store_true
        "--directory_output"
        "-d"
        help = "directory for run outputs"
        arg_type = String
        default = "data/experiment_outputs/network_flow_direct/"
        "--date"
        help = "fire-model bundle date, e.g. 2018_08_01"
        arg_type = String
        default = "2018_08_01"
        "--gaccs"
        help = "comma-separated list of GACCs"
        arg_type = String
        default = "SW"
        "--run_label"
        help = "label appended to output folders"
        arg_type = String
        default = "baseline"
        "--time_limit"
        help = "solver time limit per rolling-horizon day (seconds)"
        arg_type = Float64
        default = 300.0
        "--rest_periods"
        help = "number of consecutive rest periods required (set 0 to disable rest arcs)"
        arg_type = Int
        default = 3
    end
    return parse_args(settings)
end

function build_arc_time_lookup(
    time_space_networks::Vector{TimeSpaceNetwork},
    time_index::Int,
)
    """
    Precomputes `Dict(time_from => arc_indices)` for each TSN.  We use this
    lookup when converting committed arc sets into `fix()` constraints for past
    days so the master problem cannot change history.
    """
    lookups = Vector{Dict{Int64, Vector{Int64}}}(undef, length(time_space_networks))
    for (idx, tsn) in enumerate(time_space_networks)
        dict = Dict{Int64, Vector{Int64}}()
        arc_array = tsn.long_arcs
        for arc_ix in 1:size(arc_array, 1)
            t = arc_array[arc_ix, time_index]
            vec = get!(() -> Int64[], dict, t)
            push!(vec, arc_ix)
        end
        lookups[idx] = dict
    end
    return lookups
end

function build_fixed_arc_values(
    arc_lookup::Vector{Dict{Int64, Vector{Int64}}}, 
    committed_sets::Vector{Set{Int64}},
    cutoff_time::Int;
    active_indices::Union{Nothing,Vector{Int}} = nothing,
)
    """
    Given the arc/time lookup and the committed history sets, return a vector of
    dictionaries mapping `arc_ix => 0/1` for every arc whose `TIME_FROM` is
    strictly before `cutoff_time`.  Returning `nothing` signals that no arcs
    need to be fixed (e.g., on day 1).
    """
    if cutoff_time <= 0
        return nothing
    end
    fixed = [Dict{Int64, Float64}() for _ in arc_lookup]
    any_fixed = false
    allowed = active_indices === nothing ? nothing : Set(active_indices)
    for idx in eachindex(arc_lookup)
        if allowed !== nothing && !(idx in allowed)
            continue
        end
        local_lookup = arc_lookup[idx]
        local_fixed = fixed[idx]
        committed = committed_sets[idx]
        for (time_from, arcs) in local_lookup
            if time_from >= cutoff_time
                continue
            end
            any_fixed = true
            for arc_ix in arcs
                local_fixed[arc_ix] = arc_ix in committed ? 1.0 : 0.0
            end
        end
    end
    return any_fixed ? fixed : nothing
end

function active_fire_indices(
    fire_models::Vector{TimeSpaceNetwork},
    current_day::Int,
)
    """
    Returns a vector of fire indices that have started strictly before the end
    of the 1-indexed `current_day`. Empirical start-day metadata is zero-based
    (0 ⇒ the first simulation day), so a start period of 1 means the fire
    should first appear on day 2, and so on.
    """
    active = Int[]
    for (idx, fm) in enumerate(fire_models)
        start_period = fm.start_time_period
        if isnothing(start_period) || start_period < current_day
            push!(active, idx)
        end
    end
    return active
end

function full_network_flow(
    crew_models::Vector{TimeSpaceNetwork},
    fire_models::Vector{TimeSpaceNetwork};
    integer = true,
    verbose = false,
    time_limit = 180.0,
    output_dir::Union{Nothing,String} = nothing,
    fire_meta = nothing,
    crew_rest_deadlines::Union{Nothing,Vector{Int}} = nothing,
    fixed_fire_values::Union{Nothing,Vector{Dict{Int64, Float64}}} = nothing,
    fixed_crew_values::Union{Nothing,Vector{Dict{Int64, Float64}}} = nothing,
    active_fires::Vector{Int} = collect(1:length(fire_models)),
)
    ub = Inf
    lb = 0.0
    num_crews = length(crew_models)
    num_fires = length(fire_models)
    active_fire_list = unique(active_fires)
    # Track fires that are currently hidden so we can mask their crew arcs.
    inactive_fire_set = Set{Int}()
    if length(active_fire_list) < num_fires
        all_fires = collect(1:num_fires)
        inactive_fire_set = Set(setdiff(all_fires, active_fire_list))
    end
    _, num_times, _ = size(crew_models[1].state_in_arcs)

    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    if !verbose
        set_optimizer_attribute(m, "OutputFlag", 0)
    end
    set_optimizer_attribute(m, "TimeLimit", time_limit)

    fire_vars = Vector{Union{Nothing, Vector{VariableRef}}}(undef, num_fires)
    for fire in 1:num_fires
        fire_vars[fire] = nothing
    end
    # Fire arcs: one binary per long arc. Past-day arcs can be fixed via
    # `fixed_fire_values` so the rolling horizon keeps history immutable.
    for fire in active_fire_list
        fire_model = fire_models[fire]
        y = @variable(
            m,
            [1:size(fire_model.long_arcs, 1)],
            lower_bound = 0,
            upper_bound = 1,
        )
        if integer
            set_binary.(y)
        end
        if fixed_fire_values !== nothing
            for (arc_ix, val) in fixed_fire_values[fire]
                fix(y[arc_ix], val; force = true)
            end
        end
        fire_vars[fire] = y
    end

    crew_vars = Vector{Vector{VariableRef}}(undef, num_crews)
    # Crew arcs: identical treatment, but domains are filtered so each vector
    # only contains arcs belonging to that specific crew.
    for crew in 1:num_crews
        crew_model = crew_models[crew]
        ixs = findall(crew_model.long_arcs[:, CM.CREW_NUMBER] .== crew)
        z = @variable(
            m,
            [ixs],
            lower_bound = 0,
            upper_bound = 1,
        )
        if integer
            set_binary.(z)
        end
        if fixed_crew_values !== nothing
            for (arc_ix, val) in fixed_crew_values[crew]
                if arc_ix ∈ ixs
                    fix(z[arc_ix], val; force = true)
                end
            end
        end
        # Crew networks still cover the full geography, so explicitly zero out
        # arcs that travel to fires which have not started yet.
        if !isempty(inactive_fire_set)
            for arc_ix in ixs
                arc = crew_model.long_arcs[arc_ix, :]
                if arc[CM.TO_TYPE] == CM.FIRE_CODE
                    fire_ix = arc[CM.LOC_TO]
                    if fire_ix in inactive_fire_set
                        fix(z[arc_ix], 0.0; force = true)
                    end
                end
            end
        end
        crew_vars[crew] = z
    end

    @objective(
        m,
        Min,
        sum(
            fire_models[fire].arc_costs[ix] * fire_vars[fire][ix]
            for fire in active_fire_list,
                ix in 1:size(fire_models[fire].long_arcs, 1)
            if fire_models[fire].long_arcs[ix, FM.TIME_TO] == num_times + 1
        ),
    )

    if !isempty(active_fire_list)
        @constraint(
            m,
            fire_flow[fire = active_fire_list, t = 1:num_times, s = 1:size(fire_models[fire].state_out_arcs, 1)],
            sum(fire_vars[fire][ix] for ix in fire_models[fire].state_out_arcs[s, t]) ==
            sum(fire_vars[fire][ix] for ix in fire_models[fire].state_in_arcs[s, t]),
        )

        @constraint(
            m,
            fire_start[fire = active_fire_list],
            fire_vars[fire][1] == 1,
        )
    end

    locs, times, rests = size(crew_models[1].state_in_arcs)
    @constraint(
        m,
        crew_flow[crew = 1:num_crews, l = 1:locs, t = 1:times, r = 1:rests],
        sum(crew_vars[crew][ix] for ix in crew_models[crew].state_in_arcs[l, t, r]) ==
        sum(crew_vars[crew][ix] for ix in crew_models[crew].state_out_arcs[l, t, r]),
    )

    @constraint(
        m,
        crew_start[crew = 1:num_crews],
        sum(
            crew_vars[crew][ix]
            for ix in findall(crew_models[crew].long_arcs[:, CM.TIME_FROM] .== 0)
            if crew_models[crew].long_arcs[ix, CM.CREW_NUMBER] == crew
        ) == 1,
    )

    if !isempty(active_fire_list)
        @constraint(
            m,
            linking[fire = active_fire_list, t = 1:num_times],
            sum(
                crew_vars[crew][ix]
                for crew in 1:num_crews,
                    ix in vcat(
                        crew_models[crew].state_in_arcs[fire, t, 1],
                        crew_models[crew].state_in_arcs[fire, t, 2],
                    )
            ) >=
            sum(
                fire_models[fire].long_arcs[ix, FM.CREWS_PRESENT] * fire_vars[fire][ix]
                for ix in eachindex(fire_vars[fire])
                if fire_models[fire].long_arcs[ix, FM.TIME_FROM] == t
            ),
        )
    end

    # Solve a standard JuMP model. Rolling-horizon behavior comes entirely from
    # the fixed arc values passed in above.
    optimize!(m)
    termination = termination_status(m)

    fire_selected = nothing
    fire_selected_values = nothing
    crew_selected = nothing
    crew_selected_values = nothing

    if has_values(m)
        ub = objective_value(m)
        lb = objective_bound(m)
        solve_seconds = JuMP.solve_time(m)
        gap = (ub - lb) / max(1, abs(ub))
        @info "Solve complete" objective = ub bound = lb gap = gap time = solve_seconds

        fire_selected = [Int[] for _ in 1:num_fires]
        fire_selected_values = [Float64[] for _ in 1:num_fires]
        for fire in active_fire_list
            fire_var = fire_vars[fire]
            fire_var === nothing && continue
            vals = value.(fire_var)
            selected = [ix for ix in axes(vals, 1) if vals[ix] > 1e-6]
            fire_selected[fire] = selected
            fire_selected_values[fire] = [vals[ix] for ix in selected]
            for ix in selected
                arc = fire_models[fire].long_arcs[ix, :]
                cost = fire_models[fire].arc_costs[ix]
                raw_from = fire_models[fire].raw_state_from === nothing ?
                    get(raw_state_lookup[fire], arc[FM.STATE_FROM], arc[FM.STATE_FROM]) :
                    fire_models[fire].raw_state_from[ix]
                raw_to = fire_models[fire].raw_state_to === nothing ?
                    get(raw_state_lookup[fire], arc[FM.STATE_TO], arc[FM.STATE_TO]) :
                    fire_models[fire].raw_state_to[ix]
                personnel = arc[FM.CREWS_PRESENT] * crew_step
                @debug "Fire arc" fire = fire index = ix cost = cost time_from = arc[FM.TIME_FROM] time_to = arc[FM.TIME_TO] raw_from = raw_from raw_to = raw_to personnel = personnel
            end
        end

        crew_selected = Vector{Vector{Int64}}(undef, num_crews)
        crew_selected_values = Vector{Vector{Float64}}(undef, num_crews)
        for crew in 1:num_crews
            vals = value.(crew_vars[crew])
            selected = [ix for ix in axes(vals, 1) if vals[ix] > 1e-6]
            crew_selected[crew] = selected
            crew_selected_values[crew] = [vals[ix] for ix in selected]
        end

        # Serialize the selected fire arcs regardless of whether we export
        # per-day folders.  These are also used downstream when generating
        # per-fire CSVs.
        records = DataFrame(
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
        for fire in 1:num_fires
            vals = fire_selected_values[fire]
            selected = fire_selected[fire]
            for (pos, ix) in enumerate(selected)
                arc = fire_models[fire].long_arcs[ix, :]
                cost = fire_models[fire].arc_costs[ix]
                raw_from = fire_models[fire].raw_state_from === nothing ?
                    get(raw_state_lookup[fire], arc[FM.STATE_FROM], arc[FM.STATE_FROM]) :
                    fire_models[fire].raw_state_from[ix]
                raw_to = fire_models[fire].raw_state_to === nothing ?
                    get(raw_state_lookup[fire], arc[FM.STATE_TO], arc[FM.STATE_TO]) :
                    fire_models[fire].raw_state_to[ix]
                personnel = arc[FM.CREWS_PRESENT] * crew_step
                push!(records, (
                    fire = fire,
                    arc_index = ix,
                    value = vals[pos],
                    cost = cost,
                    state_from_raw = raw_from,
                    time_from = arc[FM.TIME_FROM],
                    time_to = arc[FM.TIME_TO],
                    state_to_raw = raw_to,
                    personnel = personnel,
                ))
            end
        end

        if isnothing(output_dir)
            CSV.write("selected_fire_arcs.csv", records)
        else
            mkpath(output_dir)
            CSV.write(joinpath(output_dir, "selected_fire_arcs.csv"), records)
            fire_count = length(fire_models)
        progression_summary = DataFrame(
            fire = Int[],
            fire_event_id = String[],
            day = Int[],
            acres = Float64[],
            crews = Float64[],
        )

        crews_on_fire_counts = Int[]
        crews_in_transit_counts = Int[]
        crews_resting_counts = Int[]
        crews_at_base_counts = Int[]
        fire_daily_counts = [Int[] for _ in 1:fire_count]

        function ensure_length!(vec::Vector{T}, len::Int) where {T}
            current = length(vec)
            if current < len
                resize!(vec, len)
                for idx in (current + 1):len
                    vec[idx] = zero(T)
                end
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
            for fire in 1:fire_count
                fire_rows = records[records.fire .== fire, :]
                nrow(fire_rows) == 0 && continue
                arc_info = fire_meta !== nothing ? fire_meta.fire_order[fire] : Dict{String,Any}()
                fire_id = get(arc_info, "fire_event_id", fire)
                arc_file = get(arc_info, "arc_file", "arc_file")
                filename = string(
                    "fire_",
                    slugify(string(fire_id)),
                    "_",
                    slugify(splitext(basename(string(arc_file)))[1]),
                    ".csv",
                )
                selected = fire_rows.arc_index
                subset = DataFrame(fire_models[fire].long_arcs[selected, :], :auto)
                rename!(
                    subset,
                    [:col1, :state_from, :time_from_model, :time_to_model, :state_to, :crews_present],
                )
                fire_export = hcat(fire_rows, subset[:, Not(:col1)])
                CSV.write(joinpath(output_dir, filename), fire_export)
                for ix in selected
                    arc = fire_models[fire].long_arcs[ix, :]
                    day = arc[FM.TIME_TO] - 1
                    acres = round(fire_models[fire].arc_costs[ix] * 1e4)
                    crews = arc[FM.CREWS_PRESENT] * crew_step / 50.0
                    push!(progression_summary, (
                        fire = fire,
                        fire_event_id = string(fire_id),
                        day = day,
                        acres = acres,
                        crews = crews,
                    ))
                end
            end
            CSV.write(joinpath(output_dir, "fire_progression_summary.csv"), progression_summary)

            for crew in 1:num_crews
                selected = crew_selected[crew]
                if isempty(selected)
                    accumulate_range!(crews_at_base_counts, 0, num_times, 1)
                    continue
                end
                vals = crew_selected_values[crew]
                crew_arc_data = DataFrame(
                    crew_models[crew].long_arcs[selected, :],
                    [:crew_number, :from_type, :loc_from, :to_type, :loc_to, :time_from, :time_to, :rest_from, :rest_to],
                )
                crew_export = DataFrame(
                    arc_index = selected,
                    value = vals,
                    cost = crew_models[crew].arc_costs[selected],
                )
                crew_export = hcat(crew_export, crew_arc_data)
                crew_filename = string("crew_", lpad(string(crew), 3, '0'), "_selected_arcs.csv")
                CSV.write(joinpath(output_dir, crew_filename), crew_export)
                post_status = Vector{String}(undef, length(selected))
                has_positive_duration = false
                rest_requirement_pending = crew_rest_deadlines === nothing ? false : (crew_rest_deadlines[crew] <= num_times)
                for (pos, ix) in enumerate(selected)
                    arc = crew_models[crew].long_arcs[ix, :]
                    status = classify_arc_status(arc)
                    if status == :rest
                        if rest_requirement_pending
                            rest_requirement_pending = false
                        else
                            status = :base
                        end
                    end
                    post_status[pos] = string(status)
                    start_day = Int(max(0, arc[CM.TIME_FROM]))
                    end_day = Int(max(0, arc[CM.TIME_TO]))
                    if end_day > start_day
                        has_positive_duration = true
                    end
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
                crew_post_export = copy(crew_export)
                crew_post_export[!, :status] = post_status
                post_filename = string("crew_", lpad(string(crew), 3, '0'), "_selected_arcs_post_processed.csv")
                CSV.write(joinpath(output_dir, post_filename), crew_post_export)
                if !has_positive_duration
                    accumulate_range!(crews_at_base_counts, 0, num_times, 1)
                end
            end

            fire_daily_records = DataFrame(
                fire = Int[],
                fire_event_id = String[],
                day = Int[],
                crews = Int[],
            )
            for (fire_idx, counts) in enumerate(fire_daily_counts)
                fire_info = fire_meta !== nothing ? fire_meta.fire_order[fire_idx] : Dict{String,Any}()
                fire_id = haskey(fire_info, "fire_event_id") ? fire_info["fire_event_id"] : fire_idx
                for (day_idx, crews) in enumerate(counts)
                    if crews <= 0
                        continue
                    end
                    push!(fire_daily_records, (
                        fire = fire_idx,
                        fire_event_id = string(fire_id),
                        day = day_idx - 1,
                        crews = crews,
                    ))
                end
            end
            CSV.write(joinpath(output_dir, "fire_daily_crews.csv"), fire_daily_records)

            max_days = max(1, num_times)
            for vec in (crews_on_fire_counts, crews_in_transit_counts, crews_resting_counts, crews_at_base_counts)
                if length(vec) > max_days
                    resize!(vec, max_days)
                else
                    ensure_length!(vec, max_days)
                end
            end
            status_df = DataFrame(
                baseline = fill("Optimizer RH", max_days),
                day_index = collect(0:(max_days - 1)),
                day_number = collect(1:max_days),
                crews_on_fire = crews_on_fire_counts,
                crews_in_transit = crews_in_transit_counts,
                crews_resting = crews_resting_counts,
                crews_at_base = crews_at_base_counts,
            )
            CSV.write(joinpath(output_dir, "crew_status.csv"), status_df)

            summary_payload = Dict(
                "objective" => round(ub * 1e4),
                "bound" => round(lb * 1e4),
                "gap" => gap,
                "solve_seconds" => solve_seconds,
            )
            summary_path = joinpath(output_dir, "run_summary.json")
            open(summary_path, "w") do io
                write(io, JSON.json(summary_payload))
            end
        end
    end

       return (
        lb = lb,
        ub = ub,
        status = termination,
        fire_selected = fire_selected,
        fire_selected_values = fire_selected_values,
        crew_selected = crew_selected,
        crew_selected_values = crew_selected_values,
    )
end

function rolling_horizon_network_flow()
    """
    Entry point that runs the rolling horizon. High-level steps:

      1. Load the crew/fire TSNs exactly like the single-shot script.
      2. Precompute arc-by-time lookups and empty committed arc sets.
      3. For each day in the horizon:
          a. Build dictionaries that fix all past arcs to 0/1.
          b. Call `full_network_flow` with those fixed values.
          c. Commit any arcs whose `TIME_FROM` < current_day.
    """
    args = get_command_line_args()
    dataset = joinpath(@__DIR__, "..", "ai_wildfire", "fire_models_" * args["date"])
    raw_gaccs = split(strip(args["gaccs"]), ',')
    target_gaccs = [String(normalize_gacc(strip(g))) for g in raw_gaccs if !isempty(strip(g))]

    gacc_slug = join([uppercase(strip(g)) for g in raw_gaccs if !isempty(strip(g))], "-")
    run_label_slug = slugify(args["run_label"])
    run_folder = string(args["date"], "_", gacc_slug, "_", run_label_slug)
    run_output_dir = joinpath(args["directory_output"], run_folder)
    mkpath(run_output_dir)

    num_time_periods = 14
    crew_speed = 40.0 * 6.0

    num_fires = count_selected_fires(target_gaccs, Dict{String,Vector{Int64}}(), dataset)
    crew_models, crew_info = build_crew_models_from_empirical(
        num_fires,
        num_time_periods,
        crew_speed;
        crew_gaccs = target_gaccs,
        fire_gaccs = target_gaccs,
        fire_folder = dataset,
        rest_periods = args["rest_periods"],
    )
    num_crews = length(crew_models)
    start_base_arc_indices = Vector{Union{Nothing,Int64}}(undef, num_crews)
    for c in 1:num_crews
        arcs = crew_models[c].long_arcs
        substitute = nothing
        for arc_ix in 1:size(arcs, 1)
            if arcs[arc_ix, CM.FROM_TYPE] == CM.BASE_CODE &&
               arcs[arc_ix, CM.TO_TYPE] == CM.BASE_CODE &&
               arcs[arc_ix, CM.TIME_FROM] == 0 &&
               arcs[arc_ix, CM.TIME_TO] == 1 &&
               arcs[arc_ix, CM.REST_FROM] == 0
                substitute = arc_ix
                break
            end
        end
        start_base_arc_indices[c] = substitute
    end
    crew_names = hasproperty(crew_info, :crew_names) ? crew_info.crew_names : nothing
    if crew_names !== nothing && !isempty(run_output_dir)
        mapping = DataFrame(
            crew_index = collect(1:length(crew_names)),
            crew_name = crew_names,
        )
        CSV.write(joinpath(run_output_dir, "crew_index_mapping.csv"), mapping)
    end

    fire_models, fire_meta = build_fire_models_from_empirical(
        num_fires,
        num_crews,
        num_time_periods;
        fire_gaccs = target_gaccs,
        fires_by_gacc = Dict{String,Vector{Int64}}(),
        fire_folder = dataset,
    )

    @info "Dataset" dataset
    @info "Target GACCs" target_gaccs
    @info "Output directory" run_output_dir
    @info "Num fires" num_fires
    @info "Num crews" num_crews

    global raw_state_lookup
    raw_state_lookup = Vector{Dict{Int, Int}}(undef, length(fire_models))
    for (f, entry) in enumerate(fire_meta.fire_order)
        state_meta = get(entry, "state_metadata", nothing)
        lookup = Dict{Int, Int}()
        if state_meta isa AbstractVector
            for sm in state_meta
                sm isa Dict || continue
                sid = Int(sm["state_id"])
                raw = get(sm, "packed_code", nothing)
                if !(raw === nothing || raw === missing)
                    lookup[sid] = Int(raw)
                else
                    lookup[sid] = sid
                end
            end
        end
        raw_state_lookup[f] = lookup
    end
    global crew_step
    crew_step = hasproperty(fire_meta, :crew_step) ? fire_meta.crew_step : 50

    fire_start_periods = [fsp.start_time_period for fsp in fire_models]
    # Prevent crews from departing toward fires before their modeled start day.
    for j in 1:num_crews
        no_fire_anticipation!(crew_models[j], fire_start_periods)
    end

    # Precompute which arc indices originate at each time stage.  This lets us
    # translate “everything before current_day must stay fixed” into `fix()` calls
    # without scanning the entire arc array every iteration.
    fire_arcs_by_time = build_arc_time_lookup(fire_models, FM.TIME_FROM)
    crew_arcs_by_time = build_arc_time_lookup(crew_models, CM.TIME_FROM)
    committed_fire_arcs = [Set{Int64}() for _ in 1:num_fires]
    committed_crew_arcs = [Set{Int64}() for _ in 1:num_crews]
    fire_has_history = falses(num_fires)

    for t in 0:num_time_periods
        current_day = t + 1
        total_days = num_time_periods + 1
        @info "##### Rolling-Horizon Day $(current_day) of $(total_days) #####"

        # UI-style logging so the user can see how the horizon is sliding.
        fires_active = active_fire_indices(fire_models, current_day)
        fires_starting_today = [
            g for g in 1:num_fires
            if fire_start_periods[g] == current_day - 1
        ]
        @info "Active fires" total = length(fires_active) fires = fires_active
        if isempty(fires_active)
            @info "No active fires yet; solver will ignore future ignitions until they start"
        end
        if !isempty(fires_starting_today)
            @info "Fires starting today" fires_starting_today
        end

        # Convert the committed sets into dictionaries passed to `fix()`.  Arcs
        # with `TIME_FROM < current_day` are forced to 1 if previously used or 0
        # if we know they were skipped, which mimics the branch-price rolling
        # horizon exactly.
        # Allow a newly started fire to plan freely on its first day by only
        # fixing history for fires that have already appeared in prior solves.
        fires_with_history = [g for g in fires_active if fire_has_history[g]]
        fixed_fire = build_fixed_arc_values(
            fire_arcs_by_time,
            committed_fire_arcs,
            t;
            active_indices = fires_with_history,
        )
        fixed_crew = build_fixed_arc_values(crew_arcs_by_time, committed_crew_arcs, t)

        day_output_dir = joinpath(run_output_dir, string("day_", lpad(string(current_day), 3, '0')))
        result = full_network_flow(
            crew_models,
            fire_models;
            integer = true,
            verbose = args["debug"],
            time_limit = args["time_limit"],
            output_dir = day_output_dir,
            fire_meta = fire_meta,
            crew_rest_deadlines = hasproperty(crew_info, :rest_by) ? crew_info.rest_by : nothing,
            fixed_fire_values = fixed_fire,
            fixed_crew_values = fixed_crew,
            active_fires = fires_active,
        )

        if result.fire_selected === nothing || result.crew_selected === nothing
            @warn "Solver did not return a feasible solution on day $(current_day); terminating rolling horizon"
            break
        end
        for g in fires_active
            fire_has_history[g] = true
        end

        # Commit any arcs whose `TIME_FROM` lies before the upcoming day so the
        # next iteration knows which decisions are immutable.
        commit_cutoff = t + 1
        for g in 1:num_fires
            for arc_ix in result.fire_selected[g]
                tf = fire_models[g].long_arcs[arc_ix, FM.TIME_FROM]
                if tf < commit_cutoff
                    push!(committed_fire_arcs[g], arc_ix)
                end
            end
        end
        for c in 1:num_crews
            for arc_ix in result.crew_selected[c]
                tf = crew_models[c].long_arcs[arc_ix, CM.TIME_FROM]
                tt = crew_models[c].long_arcs[arc_ix, CM.TIME_TO]
                if tf >= commit_cutoff
                    continue
                end
                if tf == 0 && tt == 0
                    substitute = start_base_arc_indices[c]
                    if substitute !== nothing
                        push!(committed_crew_arcs[c], substitute)
                    else
                        push!(committed_crew_arcs[c], arc_ix)
                    end
                else
                    push!(committed_crew_arcs[c], arc_ix)
                end
            end
        end
    end
end

rolling_horizon_network_flow()
