include("Subproblems.jl")  # Load time-space network builders and shared constants.

using JuMP, Gurobi, JSON, ArgParse  # Modeling, solver, serialization, and CLI parsing packages.
using DataFrames, CSV #Ryne added to output arcs

const GRB_ENV = Gurobi.Env()  # Set up a single reusable Gurobi environment.

# Optimization overview for full_network_flow:
#   Variables:
#     - y_fire[i]: continuous/binary selection variables for each long arc in every fire network.
#     - z_crew[i]: continuous/binary selection variables for each long arc tied to a particular crew.
#   Objective:
#     - Minimize total modeled suppression cost by summing fire arc costs plus crew arc costs.
#   Constraints:
#     - Fire network flow conservation keeps ingress and egress balanced per fire state and time.
#     - Fire start constraints force ignition arcs to be active in every fire network.
#     - Crew network flow conservation balances crews at each (location, time, rest state).
#     - Crew start constraints ensure each crew launches exactly one unit of flow at time zero.
#     - Linking constraints guarantee enough crew flow is assigned whenever a fire arc requires crews.

#Ryne added: copied from EmpiricalMain.jl
function count_selected_fires(
        fire_gaccs::Vector{String},
        fires_by_gacc::Dict{String,Vector{Int64}},
        input_folder::String,
)
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

function get_command_line_args()
    arg_parse_settings = ArgParseSettings()  # Initialize the argument parser configuration.
    @add_arg_table arg_parse_settings begin  # Declare the supported CLI switches.
        "--debug"  # Flag name.
        help = "run in debug mode, exposing all logging that uses @debug macro"  # Help text for --debug.
        action = :store_true  # Interpret --debug as a boolean toggle.
        "--directory_output", "-d"  # Name and short alias for the output directory argument.
        help = "directory to write outputs, must exist"  # Describe expected usage of -d.
        arg_type = String  # Parse -d as a string.
        default = "data/experiment_outputs/network_flow_direct/"  # Provide the default output path.

		#Ryne added
		"--date"
			help = "fire-model bundle date, e.g. 2018_08_01"
			arg_type = String
			default = "2018_08_01"
		"--gaccs"
			help = "comma-separated list of GACCs (abbreviations like SW,GB,SA are fine)"
			arg_type = String
			default = "SW"
		"--run_label"
			help = "label appended to output folders for bookkeeping"
			arg_type = String
			default = "baseline"
    end
    return parse_args(arg_parse_settings)  # Execute parsing and return a dictionary of arguments.
end

function full_network_flow(
	crew_models::Vector{TimeSpaceNetwork},  # Crew-side time-space networks.
	fire_models::Vector{TimeSpaceNetwork};  # Fire-side time-space networks.
	integer = true,  # Toggle between LP relaxation and MIP.
	verbose = false,  # Control solver logging intensity.
	time_limit = 180.0)  # Wall-clock limit for the solve.

	ub = Inf  # Initialize the best known feasible objective.
	lb = 0  # Initialize the best known lower bound.
	num_crews = length(crew_models)  # Count the crew networks provided.
	num_fires = length(fire_models)  # Count the fire networks provided.
	_, num_times, _ = size(crew_models[1].state_in_arcs)  # Infer horizon length from state tensors.


	# intialize model
	m = Model(() -> Gurobi.Optimizer(GRB_ENV))  # Create a JuMP model backed by Gurobi.
	if ~verbose  # Suppress solver logs unless explicitly requested.
		set_optimizer_attribute(m, "OutputFlag", 0)  # Silence Gurobi console output.
	end
	set_optimizer_attribute(m, "TimeLimit", time_limit)  # Apply the requested time limit.

	fire_vars = []  # Store handles to each fire’s variables.
	for fire ∈ 1:num_fires  # Iterate over fire networks.
		fire_model = fire_models[fire]  # Shortcut for readability.
		y = @variable(
			m,
			[1:size(fire_model.long_arcs)[1]],
			lower_bound = 0,
			upper_bound = 1,
		)  # Introduce arc-selection variables for every fire arc.
		# println(size(fire_model.long_arcs)[1])
		# for i ∈ 1:size(fire_model.long_arcs)[1]
		#     set_objective_coefficient(m, y[i], fire_model.arc_costs[i])
		# end
		if integer  # Switch to binaries if the caller requested a MIP.
			set_binary.(y)  # Enforce {0,1} domain for the fire arc decisions.
		end
		push!(fire_vars, y)  # Keep track of the variable vector so constraints can reference it.
	end

	crew_vars = []  # Store handles to each crew’s variables.
	for crew ∈ 1:num_crews  # Iterate over crew networks.
		crew_model = crew_models[crew]  # Local alias for readability.
		ixs = findall(crew_model.long_arcs[:, CM.CREW_NUMBER] .== crew)  # Identify indices belonging to this crew.
		z = @variable(
			m,
			[ixs],
			lower_bound = 0,
			upper_bound = 1,
		)  # Introduce arc-selection variables limited to this crew’s arcs.
		# println(length(crew_model.arc_costs[ixs]))
		# set_objective_coefficient.(m, z, crew_model.arc_costs[ixs])
		if integer  # Switch domain when solving the integer program.
			set_binary.(z)  # Enforce binary decisions for crew arcs.
		end
		push!(crew_vars, z)  # Store the handle for constraint construction.
	end

	@objective(
		#Ryne added
         m, 
		 Min,
         sum(
             fire_models[fire].arc_costs[ix] * fire_vars[fire][ix]
             for fire ∈ 1:num_fires,
                 ix ∈ 1:size(fire_models[fire].long_arcs)[1]
             if fire_models[fire].long_arcs[ix, FM.TIME_TO] == num_times + 1
         )
     

		# Ryne added
		# I think this sums each fire day so it is like final-snapshot-only being off and alpha = 1
		# m,
		# Min,
		# sum(
		# 	fire_models[fire].arc_costs[ix] * fire_vars[fire][ix] for
		# 	fire ∈ 1:num_fires,
		# 	ix ∈ 1:size(fire_models[fire].long_arcs)[1]
		# ) 


		#Ryne removed so that there are zero crew costs
		# + sum(
		# 	crew_models[crew].arc_costs[ix] * crew_vars[crew][ix] for
		# 	crew ∈ 1:num_crews,
		# 	ix ∈ findall(crew_models[crew].long_arcs[:, CM.CREW_NUMBER] .== crew)
		# )
	)  # Minimize combined fire and crew arc costs.

	# fire network flow
	@constraint(
		m,
		fire_flow[
			fire = 1:num_fires,
			t = 1:num_times,
			s = 1:size(fire_models[fire].state_out_arcs)[1],
		],
		sum(
			fire_vars[fire][ix] for ix ∈ fire_models[fire].state_out_arcs[s, t]
		) ==
		sum(fire_vars[fire][ix] for ix ∈ fire_models[fire].state_in_arcs[s, t])
	)  # Enforce flow balance at every fire state and time.

	# fire start
	@constraint(
		m,
		fire_start[fire = 1:num_fires],
		fire_vars[fire][1] == 1
	)  # Ensure each fire network ignites (arc index 1 represents ignition).

	# crew network flow
	locs, times, rests = size(crew_models[1].state_in_arcs)  # Extract dimensions for crew state tensors.
	@constraint(
		m,
		crew_flow[crew = 1:num_crews, l = 1:locs, t = 1:times, r = 1:rests],
		sum(crew_vars[crew][ix] for ix ∈ crew_models[crew].state_in_arcs[l, t, r]) ==
		sum(crew_vars[crew][ix] for ix ∈ crew_models[crew].state_out_arcs[l, t, r])
	)  # Preserve crew flow at every (location, time, rest) state.
	# crew start
	@constraint(
		m,
		crew_start[crew = 1:num_crews],
		sum(
			crew_vars[crew][ix] for
			ix ∈ findall(crew_models[crew].long_arcs[:, CM.TIME_FROM] .== 0) if
			crew_models[crew].long_arcs[ix, CM.CREW_NUMBER] == crew)
		==
		1
	)  # Force each crew to dispatch exactly once at time 0.

	# linking
	@constraint(
		m,
		linking[fire = 1:num_fires, t = 1:num_times],
		sum(
			crew_vars[crew][ix] for crew ∈ 1:num_crews,
			ix ∈ vcat(
				crew_models[crew].state_in_arcs[fire, t, 1],
				crew_models[crew].state_in_arcs[fire, t, 2],
			)
		) >=
		sum(
			fire_models[fire].long_arcs[ix, FM.CREWS_PRESENT] .* fire_vars[fire][ix]
			for ix ∈ eachindex(fire_vars[fire]) if
			fire_models[fire].long_arcs[ix, FM.TIME_FROM] == t
		)
	)  # Guarantee enough crew flow exists whenever a fire arc needs crews at time t.

	optimize!(m)  # Solve the constructed JuMP model.

	if has_values(m)  # Only read solution values when the solver finished cleanly.
		ub = objective_value(m)  # Pull the incumbent objective value.
		lb = objective_bound(m)  # Pull the reported dual bound.
		solve_seconds = JuMP.solve_time(m)

		#Ryne added
		@info "Solve complete cleanly" objective=ub bound=lb gap=(ub - lb)/max(1, abs(ub)) time=solve_seconds

		for fire in 1:num_fires
			vals = value.(fire_vars[fire])
			selected = findall(>(1e-6), vals)
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
				@info "Fire arc" fire=fire index=ix value=vals[ix] cost=cost arc = repr(arc) data=(
					state_from_raw = raw_from,
					time_from = arc[FM.TIME_FROM],
					time_to = arc[FM.TIME_TO],
					state_to_raw = raw_to,
					personnel = personnel,
				)
			end
		end
		# for crew in 1:num_crews
		# 	vals = value.(crew_vars[crew])
		# 	selected = findall(>(1e-6), vals)
		# 	for ix in selected
		# 		arc = crew_models[crew].long_arcs[ix, :]
		# 		cost = crew_models[crew].arc_costs[ix]
		# 		@info "Crew arc" crew=crew index=ix value=vals[ix] cost=cost data=arc
		# 	end
		# end
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
			vals = value.(fire_vars[fire])
			selected = findall(>(1e-6), vals)
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
				push!(records, (
					fire = fire,
					arc_index = ix,
					value = vals[ix],
					cost = cost,
					state_from_raw = raw_from,
					time_from = arc[FM.TIME_FROM],
					time_to = arc[FM.TIME_TO],
					state_to_raw = raw_to,
					personnel = personnel,
				))
			end
		end
		@info "Solve complete cleanly" objective=ub.*1e4 bound=lb.*1e4 gap=(ub - lb)/max(1, abs(ub)) time=solve_seconds

		CSV.write("selected_fire_arcs.csv", records)
	end

	return lb, ub  # Return both bounds to the caller.

end
args = get_command_line_args()  # Parse CLI configuration once on startup.


dataset = joinpath(@__DIR__, "..", "ai_wildfire", "fire_models_" * args["date"])
raw_gaccs = split(strip(args["gaccs"]), ',')
target_gaccs = [String(normalize_gacc(strip(g))) for g in raw_gaccs if !isempty(strip(g))]
# target_gaccs = ["Southwest"]
println(dataset)
println(target_gaccs)




# # Ryne added

# dataset = joinpath(@__DIR__, "..", "ai_wildfire", "fire_models_2018_08_01")
# target_gaccs = ["Southwest"]              # any set of GACCs
num_time_periods = 14                     # planning horizon you want
crew_speed = 40.0 * 6.0                        # keep or change depending on study

num_fires = count_selected_fires(target_gaccs, Dict{String,Vector{Int64}}(), dataset)

crew_models, _ = build_crew_models_from_empirical(
	num_fires,
	num_time_periods,
	crew_speed;
	crew_gaccs = target_gaccs,
	fire_gaccs = target_gaccs,
	fire_folder = dataset,
)
num_crews = length(crew_models)

fire_models, fire_meta = build_fire_models_from_empirical(
	num_fires,
	num_crews,
	num_time_periods;
	fire_gaccs = target_gaccs,
	fires_by_gacc = Dict{String,Vector{Int64}}(),
	fire_folder = dataset,
)

raw_state_lookup = Vector{Dict{Int,Int}}(undef, length(fire_models))
for (f, entry) in enumerate(fire_meta.fire_order)
	state_meta = get(entry, "state_metadata", nothing)
	lookup = Dict{Int,Int}()
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
crew_step = hasproperty(fire_meta, :crew_step) ? fire_meta.crew_step : 50  #firefighters per crew

full_network_flow(crew_models, fire_models, verbose = false, integer = true, time_limit = 1200)
