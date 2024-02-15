using DataFrames, CSV, DelimitedFiles

module CrewModelFactory

struct LocationAndRestStatus

	rest_by::Vector{Int64}
	current_fire::Vector{Int64}
	rested_periods::Vector{Int64}

end

struct DistancesAndTravelTimes

	ff_dist::Matrix{Float64}
	bf_dist::Matrix{Float64}
	ff_tau::Matrix{Int64}
	bf_tau::Matrix{Int64}

end

# indices for each crew arc
CREW_NUMBER = 1
FROM_TYPE = 2
LOC_FROM = 3
TO_TYPE = 4
LOC_TO = 5
TIME_FROM = 6
TIME_TO = 7
REST_FROM = 8
REST_TO = 9

# integer lookup for "FIRE" and "BASE"
FIRE_CODE = 1
BASE_CODE = 2

function generate_arcs(
	dists_and_times::DistancesAndTravelTimes,
	crew_status::LocationAndRestStatus,
	num_crews,
	num_fires,
	num_time_periods,
)

	# get fire-to-fire arcs
	ff = [
		[
			c,
			FIRE_CODE,
			f_from,
			FIRE_CODE,
			f_to,
			t_from,
			t_from + dists_and_times.ff_tau[f_to, f_from],
			rest,
			rest,
		]
		for c ∈ 1:num_crews, f_from ∈ 1:num_fires, f_to ∈ 1:num_fires,
		t_from ∈ 1:num_time_periods, rest ∈ 0:1
	]
	ff = copy(reduce(hcat, ff)')

	# get fire-to-fire arcs from start, based on current crew locations
	from_start_ff = [
		[
			c,
			FIRE_CODE,
			crew_status.current_fire[c],
			FIRE_CODE,
			f_to,
			0,
			dists_and_times.ff_tau[f_to, crew_status.current_fire[c]],
			0,
			0,
		]
		for
		c ∈ 1:num_crews, f_to ∈ 1:num_fires if crew_status.current_fire[c] != -1
	]
	from_start_ff = copy(reduce(hcat, from_start_ff)')

	# get base-to-fire arcs
	rf = [
		[
			c,
			BASE_CODE,
			c,
			FIRE_CODE,
			f_to,
			t_from,
			t_from + dists_and_times.bf_tau[c, f_to],
			rest,
			rest,
		]
		for c ∈ 1:num_crews, f_to ∈ 1:num_fires, t_from ∈ 1:num_time_periods,
		rest ∈ 0:1
	]
	rf = copy(reduce(hcat, rf)')

	# get base-to-fire arcs from start
	from_start_rf = [
		[c, BASE_CODE, c, FIRE_CODE, f_to, 0, dists_and_times.bf_tau[c, f_to], 0, 0]
		for
		c ∈ 1:num_crews, f_to ∈ 1:num_fires if crew_status.current_fire[c] == -1
	]
	from_start_rf = copy(reduce(hcat, from_start_rf)')

	# get fire-to-base arcs
	fr = [
		[
			c,
			FIRE_CODE,
			f_from,
			BASE_CODE,
			c,
			t_from,
			t_from + dists_and_times.bf_tau[c, f_from],
			rest,
			rest,
		]
		for c ∈ 1:num_crews, f_from ∈ 1:num_fires, t_from ∈ 1:num_time_periods,
		rest ∈ 0:1
	]
	fr = copy(reduce(hcat, fr)')

	# get fire-to-base arcs from start, based on cs.current crew locations
	from_start_fr = [
		[
			c,
			FIRE_CODE,
			crew_status.current_fire[c],
			BASE_CODE,
			c,
			0,
			dists_and_times.bf_tau[c, crew_status.current_fire[c]],
			0,
			0,
		]
		for c ∈ 1:num_crews if crew_status.current_fire[c] != -1
	]
	from_start_fr = copy(reduce(hcat, from_start_fr)')

	# get base-to-base arcs
	rr = [
		[
			c,
			BASE_CODE,
			c,
			BASE_CODE,
			c,
			t_from,
			t_from + 1 + (BREAK_LENGTH - 1) * rest,
			0,
			rest,
		]
		for c ∈ 1:num_crews, t_from ∈ 1:num_time_periods, rest ∈ 0:1
	]
	rr = copy(reduce(hcat, rr)')
	rr_rested = [
		[c, BASE_CODE, c, BASE_CODE, c, t_from, t_from + 1, 1, 1]
		for c ∈ 1:num_crews, t_from ∈ 1:num_time_periods
	]
	rr_rested = copy(reduce(hcat, rr_rested)')

	# get base-to-base arcs from start, based on cs.current days rested
	from_start_rr = [
		[c, BASE_CODE, c, BASE_CODE, c, 0,
			1 + (BREAK_LENGTH - max(crew_status.rested_periods[c], 0) - 1) * rest, 0,
			rest]
		for c ∈ 1:num_crews, rest ∈ 0:1 if crew_status.current_fire[c] == -1
	]
	from_start_rr = copy(reduce(hcat, from_start_rr)')

	A = vcat(
		ff,
		from_start_ff,
		rf,
		from_start_rf,
		fr,
		from_start_fr,
		rr,
		rr_rested,
		from_start_rr,
	)

	return A
end

function get_distance(from_type, from_ix, to_type, to_ix, fire_fire, base_fire)

	dist = 0

	# if fire to fire
	if (from_type == FIRE_CODE) & (to_type == FIRE_CODE)
		dist = fire_fire[from_ix, to_ix]

		# if fire to base
	elseif (from_type == FIRE_CODE) & (to_type == BASE_CODE)
		dist = base_fire[to_ix, from_ix]

		# if base to fire
	elseif (from_type == BASE_CODE) & (to_type == FIRE_CODE)
		dist = base_fire[from_ix, to_ix]

		# otherwise dist still 0
	end

	return dist
end

function get_static_crew_arc_costs(gd, arcs, cost_param_dict)

	# get number of arcs
	n_arcs = length(arcs[:, 1])

	# initialize costs to 0
	costs = zeros(n_arcs)

	# if there is travel cost per mile
	if "cost_per_mile" in keys(cost_param_dict)

		# find the miles for each arc
		miles_per_arc = [
			get_distance(arcs[i, 2], arcs[i, 3],
				arcs[i, 4], arcs[i, 5],
				gd.ff_dist, gd.bf_dist) for i in 1:n_arcs
		]
		# add to costs
		costs = costs .+ (cost_param_dict["cost_per_mile"] * miles_per_arc)
	end

	# if there are rest violations
	if "rest_violation" in keys(cost_param_dict)

		# find the rest violation scores
		rest_violation_matrix = cost_param_dict["rest_violation"]
		rest_violations = [
			(arcs[i, 8] == 0) & (arcs[i, 6] > 0) ?
			rest_violation_matrix[arcs[i, 1], arcs[i, 6]] : 0
			for i in 1:n_arcs
		]

		# add to costs
		costs = costs .+ rest_violations
	end

	if "fight_fire" in keys(cost_param_dict)
		costs =
			costs .+ [
				(arcs[i, 4] == FIRE_CODE) ? cost_param_dict["fight_fire"] : 0
				for i in 1:n_arcs
			]
	end

	return copy(costs)
end

function crew_data_from_path(path)

	# get distance from fire f to fire g 
	fire_dists = readdlm(path * "/fire_distances.csv", ',')

	# get distance from base c to fire g (NUM_CREWS-by-NUM_FIRES)
	base_fire_dists = readdlm(path * "/base_fire_distances.csv", ',')

	# initialize travel times (number of periods) from fire f to fire g
	tau = convert(Array{Int}, ones(size(fire_dists)))

	# initialize number of periods to travel from base c to fire g (NUM_CREWS-by-NUM_FIRES)
	tau_base_to_fire = convert(Array{Int}, ones((size(base_fire_dists))))

	# read intial crew statuses (location, period by which they must rest)
	# (-1 in current_fire means crew is currently at base)
	# (rested_periods is the amount of time crew has been at base, relevant for completing rest)
	crew_starts = CSV.read(path * "/sample_crew_starts.csv", DataFrame)
	rest_by = crew_starts[!, "rest_by"]
	current_fire = crew_starts[!, "current_fire"]
	rested_periods = crew_starts[!, "rested_periods"]


	return (GlobalData(fire_dists, base_fire_dists, tau, tau_base_to_fire),
		CrewStatus(rest_by, current_fire, rested_periods))
end

end

module FireModelFactory

end
