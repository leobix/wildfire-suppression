include("CommonStructs.jl")

function get_adjusted_crew_arc_costs(
	long_arcs::Matrix{Float64},
	linking_duals::Matrix{Float64},
	branching_rules::Vector{CrewSupplyBranchingRule},
)

	FIRE_CODE_ = 1
	TO_TYPE_ = 4
	LOC_TO_ = 5
	TIME_TO_ = 7

	time_periods = size(linking_duals)[2]
	n_arcs = size(long_arcs)[1]

	costs = [
		(
			(long_arcs[i, TO_TYPE_] == FIRE_CODE_) &
			(long_arcs[i, TIME_TO_] <= time_periods)
		) ?
		-linking_duals[long_arcs[i, LOC_TO_], long_arcs[i, TIME_TO_]] : 0
		for i in 1:n_arcs
	]
	prohibited = []

	for rule in branching_rules
		error("Not implemented")
	end

	return costs, prohibited

end

function crew_dp_inner_loop(arc, arc_ix, this_arc_cost, path_costs, min_cost, min_index) 

    BASE_CODE_ = 2
    FROM_TYPE_ = 2
    LOC_FROM_ = 3
    TIME_FROM_ = 6
    REST_FROM_ = 8

    # get the time from which this arc comes
    time_from = arc[TIME_FROM_]

    # arcs from time 0 have no past state cost
    past_state_cost = 0

    # but for other arcs
    if time_from >= 1

        # get the info of where the arc came from
        from_type = arc[FROM_TYPE_]
        loc_from = arc[LOC_FROM_]
        rest_from = arc[REST_FROM_] + 1

        # if we came from a base, that's the last row in the state matrix
        if from_type == BASE_CODE_
            past_state_cost = path_costs[size(path_costs)[1], time_from, rest_from]

        # otherwise it is the fire index row
        else
            past_state_cost = path_costs[loc_from, time_from, rest_from]   
        end
    end

    # find the path cost, update min cost and index if needed
    
    possible_cost = this_arc_cost + past_state_cost
    if possible_cost < min_cost
        min_cost = possible_cost
        min_index = arc_ix
    end

    return min_cost, min_index
end


function crew_dp_subproblem(
	arcs::Matrix{Int64},
	arc_costs::Vector{Int64},
	prohibited_arcs::Vector{Int64},
	state_in_arcs::Matrix{Vector{Int64}},
)
    """ Probably this should be refactored so the matrix is state * time
    and then we could generalize code between fire and crew subproblem"""

	BASE_CODE_ = 2
	FROM_TYPE_ = 2
	LOC_FROM_ = 3
	TIME_FROM_ = 6
	REST_FROM_ = 8

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
				min_cost = Inf
				min_index = -1

				# for each arc entering this state
				for arc_ix ∈ state_in_arcs[l, t, r]
					if arc_ix ∉ prohibited_arcs

						arc = @view arcs[:, arc_ix]
						this_arc_cost = arc_costs[arc_ix]
						min_cost, min_index = crew_dp_inner_loop(
							arc,
							arc_ix,
							this_arc_cost,
							path_costs,
							min_cost,
							min_index,
						)

					end
				end

				# store state shortest path and cost
				path_costs[l, t, r] = min_cost
				in_arcs[l, t, r] = min_index

			end
		end
	end

	lowest_cost, end_index = findmin(path_costs[:, times, :])

	# if we did not find a path, return
	if lowest_cost == Inf
		return Inf, []

	# else, find it
	else
		full_min_index = (Tuple(end_index)[1], times, Tuple(end_index)[2])

		# starting at the end state
		current_state = full_min_index
		arcs_used = []

		# while we are not at the start state
		while (current_state[2] != 0)

			# find the arc we used to get to the current state
			arc_ix = in_arcs[current_state...]
			push!(arcs_used, arc_ix)

			# find where that arc came from
			arc = @view arcs[:, arc_ix]
			time_from = arc[TIME_FROM_]
			from_type = arc[FROM_TYPE_]
			loc_from = arc[LOC_FROM_]
			rest_from = arc[REST_FROM_] + 1

			# update to that state
			if from_type == BASE_CODE_
				current_state = (locs, time_from, rest_from)
			else
				current_state = (loc_from, time_from, rest_from)
			end
		end
	end

	return lowest_cost, arcs_used

end

function get_fires_fought(
	wide_arcs::Matrix{Float64},
	arcs_used::Vector{Int64},
	(num_fires, num_time_periods)::Tuple{Int64},
)

	FIRE_CODE_ = 1
	TO_TYPE_ = 4
	LOC_TO_ = 5
	TIME_TO_ = 7

	# initialize fires fought matrix
	fires_fought = zeros(Int, num_fires, num_time_periods)

	# for each arc used
	for ix in arcs_used
		arc = @view wide_arcs[:, ix]

		# update fires_fought
		if (arc[TO_TYPE_] == FIRE_CODE_) & (arc[TIME_TO_] <= time_periods)
			fires_fought[arc[LOC_TO_], arc[TIME_TO_]] += 1
		end
	end

	return fires_fought

end

function get_adjusted_fire_arc_costs(
	long_arcs,
	linking_duals,
	branching_rules,
)
	TIME_FROM_ = 3 
	CREWS_PRESENT_ = 6

	# no cost for starting arc
	duals = vcat(0.0, linking_duals)

    # + 1 is because we appended the 0
	rel_costs = duals[long_arcs[:, TIME_FROM_].+1] .* long_arcs[:, CREWS_PRESENT_]
    prohibited_arcs = []

    for rule in branching_rules
        error("not implemented, see DCG.jl#1825")
    end

    return rel_costs, prohibited_arcs

end

function fire_dp_inner_loop(arc, arc_ix, this_arc_cost, path_costs, min_cost, min_index) 

    TIME_FROM_ = 3
    STATE_FROM_ = 2 

    # get the time from which this arc comes
    time_from = arc[TIME_FROM_]

    # arcs from time 0 have no past state cost
    past_state_cost = 0
    
    # but for other arcs
    if time_from >= 1

        # get the info of where the arc came from
        state_from = arc[STATE_FROM_]
        time_from = arc[TIME_FROM_]

        # get the min cost path to the prior state
        past_state_cost = path_costs[state_from, time_from]   
    end

    # find the path cost, update min cost and index if needed
    possible_cost = this_arc_cost + past_state_cost
    if possible_cost < min_cost
        min_cost = possible_cost
        min_index = arc_ix
    end

    return min_cost, min_index
end

function fire_dp_subproblem(arcs::Matrix{Int64},
	arc_costs::Vector{Int64},
	prohibited_arcs::Vector{Int64},
	state_in_arcs::Matrix{Vector{Int64}})

    TIME_FROM_ = 3
    STATE_FROM_ = 2 

    path_costs = zeros(Float64, size(state_in_arcs)) .+ 1e30
    in_arcs = zeros(Int, size(state_in_arcs))
    states, times = size(state_in_arcs)
    
    # iterate over times first for algorithm correctness
    for t in 1:times
        for s in 1:states
            min_cost = Inf
            min_index = -1

            # for each arc entering this state
            for arc_ix in state_in_arcs[s, t]
                if arc_ix ∉ prohibited_arcs
                    arc = @view arcs[:, arc_ix]
                    this_arc_cost = arc_costs[arc_ix]
                    min_cost, min_index = fire_dp_inner_loop(arc, arc_ix, this_arc_cost, path_costs, min_cost, min_index)
                end
            end

            # store state shortest path and cost
            path_costs[s, t] = min_cost
            in_arcs[s, t] = min_index        
        end
    end

    lowest_cost, end_index = findmin(path_costs[:, times, :])

    # if we did not find a path, return
	if lowest_cost == Inf
		return Inf, []

	# else, find it
	else
        full_min_index = (Tuple(end_index)[1], times, Tuple(end_index)[2])
        
        current_state = full_min_index
        arcs_used = []
        while (current_state[2] != 0)
            arc_ix = in_arcs[current_state...]
            push!(arcs_used, arc_ix)
            arc = @view arcs[:, arc_ix]
            state_from = arc[STATE_FROM_]
            time_from = arc[TIME_FROM_]
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

    CREWS_PRESENT_ = 6
    return reverse(wide_arcs[CREWS_PRESENT_, arcs_used][1:num_time_periods])

end
