include("CommonStructs.jl")
include("TSNetworkGeneration.jl")
include("BranchingRules.jl")

function cut_adjust_arc_costs!(costs::Vector{Float64}, cut_arc_lookup, cut_duals)

	# for each cut
	for ix in eachindex(cut_duals)

		# if the dual cost is strictly positive
		if cut_duals[ix] > 0

			# grab the sparse representation of the dual adjustment
			# Note: cost_ix = (arc_ix, coeff)

			# for each arc to be adjusted
			for cost_ix ∈ cut_arc_lookup[ix]

				# add the dual cost 
				costs[cost_ix[1]] += cut_duals[ix] * cost_ix[2]
			end
		end
	end
end

function adjust_fire_sp_arc_costs(
	branching_rule::GlobalFireAllotmentBranchingRule,
	fire::Int64,
	orig_costs::Vector{Float64},
	rule_dual_value::Float64)

	# copy the costs
	adj_costs = copy(orig_costs)

	# for each affected arc, add the dual value
	# TODO verify that >= convention means we are adding dual value either way
	for arc in branching_rule.fire_sp_arc_lookup[fire]
		adj_costs[arc] += rule_dual_value
	end

	return adj_costs

end


function adjust_crew_arc_costs!!(
	costs::Vector{Float64},
	prohibited_arcs::BitVector,
	long_arcs::Matrix{Int64},
	linking_duals::Matrix{Float64},
	linking_dual_arc_lookup::Matrix{Vector{Int64}},
	branching_rules::Vector{CrewAssignmentBranchingRule},
)

	fires, time_periods = size(linking_duals)

	for g ∈ 1:fires
		for t ∈ 1:time_periods
			for i ∈ linking_dual_arc_lookup[g, t]
				costs[i] -= linking_duals[g, t]
			end
		end
	end
			

	# get disallowed arcs due to branching rules
	# TODO refactor to track this info in B-and-B tree, check each rule just once
	for rule ∈ branching_rules
		for arc ∈ 1:size(long_arcs)[1]
			if (long_arcs[arc, CM.CREW_NUMBER] == rule.crew_ix) &&
			   (long_arcs[arc, CM.TIME_TO] == rule.time_ix)
				suppresses_fire =
					(long_arcs[arc, CM.TO_TYPE] == CM.FIRE_CODE) &&
					(long_arcs[arc, CM.LOC_TO] == rule.fire_ix)
				if ~satisfies_branching_rule(rule, suppresses_fire)
					prohibited_arcs[arc] = true
				end
			end
		end
	end
end

function crew_dp_inner_loop!!(
	min_cost::Ref{Float64},
	min_index::Ref{Int64},
	arc::SubArray{Int64, 1, Matrix{Int64}},
	arc_ix::Int64,
	this_arc_cost::Float64,
	path_costs::Array{Float64, 3},
)

	# get the time from which this arc comes
	time_from = arc[CM.TIME_FROM]

	# arcs from time 0 have no past state cost
	past_state_cost = 0.0

	# but for other arcs
	if time_from >= 1

		# get the info of where the arc came from
		from_type = arc[CM.FROM_TYPE]
		loc_from = arc[CM.LOC_FROM]
		rest_from = arc[CM.REST_FROM] + 1

		# if we came from a base, that's the last row in the state matrix
		if from_type == CM.BASE_CODE
			past_state_cost = path_costs[size(path_costs)[1], time_from, rest_from]

			# otherwise it is the fire index row
		else
			past_state_cost = path_costs[loc_from, time_from, rest_from]
		end
	end

	# find the path cost, update min cost and index if needed

	possible_cost = this_arc_cost + past_state_cost
	if possible_cost < min_cost[]
		min_cost[] = possible_cost
		min_index[] = arc_ix
	end
end


function crew_dp_subproblem(
	arcs::Matrix{Int64},
	arc_costs::Vector{Float64},
	prohibited_arcs::BitVector,
	state_in_arcs::Array{Vector{Int64}, 3},
)
	""" Probably this could be refactored so the matrix is state * time
	and then we could generalize code between fire and crew subproblem"""

	# initialize path costs to all states
	path_costs = zeros(size(state_in_arcs)) .+ Inf

	# initialize arc taken to each state
	in_arcs = zeros(Int, size(state_in_arcs))

	# get dimensions
	locs, times, rests = size(state_in_arcs)

	# iterate over times first for algorithm correctness, locs last for performance
	for t in 1:times
		for r in 1:rests
			for l in 1:locs
				min_cost = Ref{Float64}(Inf)
				min_index = Ref{Int64}(-1)

				# for each arc entering this state
				for arc_ix ∈ state_in_arcs[l, t, r]
					if ~prohibited_arcs[arc_ix]

						arc = @view arcs[:, arc_ix]
						this_arc_cost = arc_costs[arc_ix]
						crew_dp_inner_loop!!(
							min_cost,
							min_index,
							arc,
							arc_ix,
							this_arc_cost,
							path_costs,
						)

					end
				end

				# store state shortest path and cost
				path_costs[l, t, r] = min_cost[]
				in_arcs[l, t, r] = min_index[]

			end
		end
	end

	lowest_cost, end_index = findmin(path_costs[:, times, :])

	# if we did not find a path, return
	if lowest_cost == Inf
		return Inf, Int64[]

		# else, find it
	else
		full_min_index = (Tuple(end_index)[1], times, Tuple(end_index)[2])

		# starting at the end state
		current_state = full_min_index
		arcs_used = Int64[]

		# while we are not at the start state
		while (current_state[2] != 0)

			# find the arc we used to get to the current state
			arc_ix = in_arcs[current_state...]
			push!(arcs_used, arc_ix)

			# find where that arc came from
			arc = @view arcs[:, arc_ix]
			time_from = arc[CM.TIME_FROM]
			from_type = arc[CM.FROM_TYPE]
			loc_from = arc[CM.LOC_FROM]
			rest_from = arc[CM.REST_FROM] + 1

			# update to that state
			if from_type == CM.BASE_CODE
				current_state = (locs, time_from, rest_from)
			else
				current_state = (loc_from, time_from, rest_from)
			end
		end
	end

	return lowest_cost, arcs_used

end

function get_fires_fought(
	wide_arcs::Matrix{Int64},
	arcs_used::Vector{Int64},
	(num_fires, num_time_periods)::Tuple{Int64, Int64},
)

	# initialize fires fought bit-matrix
	fires_fought = falses(num_fires, num_time_periods)

	# for each arc used
	for ix in arcs_used
		arc = @view wide_arcs[:, ix]

		# update fires_fought
		if (arc[CM.TO_TYPE] == CM.FIRE_CODE) & (arc[CM.TIME_TO] <= num_time_periods)
			fires_fought[arc[CM.LOC_TO], arc[CM.TIME_TO]] += 1
		end
	end

	return fires_fought

end

function get_adjusted_fire_arc_costs(
	long_arcs,
	linking_duals,
	branching_rules,
)

	# no cost for starting arc
	duals = vcat(0.0, linking_duals)

	# + 1 is because we appended the 0
	rel_costs =
		duals[long_arcs[:, FM.TIME_FROM].+1] .* long_arcs[:, FM.CREWS_PRESENT]

	# get disallowed arcs due to branching rules
	# TODO refactor to track this info in B-and-B tree, check each rule just once
	prohibited_arcs = falses(size(long_arcs)[1])
	for rule ∈ branching_rules
		for arc ∈ 1:size(long_arcs)[1]
			if long_arcs[arc, FM.TIME_FROM] == rule.time_ix
				if ~satisfies_branching_rule(rule, long_arcs[arc, FM.CREWS_PRESENT])
					prohibited_arcs[arc] = true
				end
			end
		end
	end

	return rel_costs, prohibited_arcs

end

function fire_dp_inner_loop!!(
	min_cost::Ref{Float64},
	min_index::Ref{Int64},
	arc::SubArray{Int64, 1, Matrix{Int64}},
	arc_ix::Int64,
	this_arc_cost::Float64,
	path_costs::Matrix{Float64},
)

	# get the time from which this arc comes
	time_from = arc[FM.TIME_FROM]

	# get the past state cost
	past_state_cost =
		time_from == 0 ? 0.0 : path_costs[arc[FM.STATE_FROM], time_from]

	# find the path cost, update min cost and index if needed
	possible_cost = this_arc_cost + past_state_cost
	if possible_cost < min_cost[]
		min_cost[] = possible_cost
		min_index[] = arc_ix
	end

end

function fire_dp_subproblem(arcs::Matrix{Int64},
	arc_costs::Vector{Float64},
	prohibited_arcs::BitVector,
	state_in_arcs::Matrix{Vector{Int64}})

	path_costs = zeros(Float64, size(state_in_arcs)) .+ 1e30
	in_arcs = zeros(Int, size(state_in_arcs))
	states, times = size(state_in_arcs)

	# iterate over times first for algorithm correctness
	for t in 1:times
		for s in 1:states
			min_cost = Ref{Float64}(Inf)
			min_index = Ref{Int64}(-1)

			# for each arc entering this state
			for arc_ix in state_in_arcs[s, t]
				if ~prohibited_arcs[arc_ix]
					arc = @view arcs[:, arc_ix]
					this_arc_cost = arc_costs[arc_ix]
					fire_dp_inner_loop!!(
						min_cost,
						min_index,
						arc,
						arc_ix,
						this_arc_cost,
						path_costs,
					)
				end
			end

			# store state shortest path and cost
			path_costs[s, t] = min_cost[]
			in_arcs[s, t] = min_index[]
		end
	end

	lowest_cost, end_index = findmin(path_costs[:, times, :])

	# if we did not find a path, return
	if lowest_cost == Inf
		return Inf, Int64[]

		# else, find it
	else
		full_min_index = (Tuple(end_index)[1], times, Tuple(end_index)[2])

		current_state = full_min_index
		arcs_used = Int64[]
		while (current_state[2] != 0)
			arc_ix = in_arcs[current_state...]
			push!(arcs_used, arc_ix)
			arc = @view arcs[:, arc_ix]
			state_from = arc[FM.STATE_FROM]
			time_from = arc[FM.TIME_FROM]
			current_state = (state_from, time_from)
		end
		return lowest_cost, arcs_used
	end
end


function get_crew_demands(
	wide_arcs,
	arcs_used,
	num_time_periods,
)

	FM.CREWS_PRESENT = 6
	return reverse(wide_arcs[FM.CREWS_PRESENT, arcs_used][1:num_time_periods])

end
