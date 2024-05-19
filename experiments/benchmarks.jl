include("../BranchingRules.jl")
include("../Subproblems.jl")

using JuMP, Gurobi

const GRB_ENV = Gurobi.Env()

function backward_crew_dp_subproblem(
    arcs::Matrix{Int64},
    prohibited_arcs::BitVector,
    state_in_arcs::Array{Vector{Int64},3},
)
    """ Probably this could be refactored so the matrix is state * time
    and then we could generalize code between fire and crew subproblem"""

    # initialize path costs to all states
    reachable = falses(size(state_in_arcs))

    # get dimensions
    locs, times, rests = size(state_in_arcs)

    # iterate over times first for algorithm correctness, locs last for performance
    for t in reverse(1:times)
        for r in 1:rests
            for l in 1:locs

                if t == times
                    reachable[l, t, r] = true
                end

                if reachable[l, t, r]

                    # for each arc entering this state
                    for arc_ix ∈ state_in_arcs[l, t, r]
                        if ~prohibited_arcs[arc_ix]

                            arc = @view arcs[:, arc_ix]

                            # get the time from which this arc comes
                            time_from = arc[CM.TIME_FROM]
                            loc_from = (arc[CM.FROM_TYPE] == CM.FIRE_CODE) ? arc[CM.LOC_FROM] : locs
                            rest_from = arc[CM.REST_FROM] + 1

                            if time_from > 0
                                reachable[loc_from, time_from, rest_from] = true
                            end

                        end
                    end
                end
            end
        end
    end
    return reachable
end

function no_suppression_fire_damage(
	fire_model::TimeSpaceNetwork,
	num_times::Int,
	branching_rules::Vector{FireDemandBranchingRule},
)

    mock_duals = zeros(Float64, num_times)
    max_time = maximum(vcat([rule.time_ix for rule ∈ branching_rules], 0))
    for i ∈ 1:num_times
        if i > max_time
            mock_duals[i] = 1e30
        end
    end

	rel_costs, prohibited_arcs = get_adjusted_fire_arc_costs(
		fire_model.long_arcs,
		mock_duals,
		branching_rules,
	)
	arc_costs = rel_costs .+ fire_model.arc_costs

	objective, _ = fire_dp_subproblem(
		fire_model.wide_arcs,
		arc_costs,
		prohibited_arcs,
		fire_model.state_in_arcs,
	)

	return objective
end

function triage_then_route_by_time_period(
	crew_models::Vector{TimeSpaceNetwork},
	fire_models::Vector{TimeSpaceNetwork},
	λ::Float64,
	β::Float64,
)

	# get problem size
	num_crews = length(crew_models)
	num_fires = length(fire_models)
	_, num_times, _ = size(crew_models[1].state_in_arcs)

	# initialize arcs traversed
	crew_arcs_traversed = [Int[] for crew ∈ 1:num_crews]
	fire_arcs_traversed = [Int[] for fire ∈ 1:num_fires]

	# intitalize crew states
	crew_time_from = Vector{Int}(undef, num_crews)
	crew_from_type = Vector{Int}(undef, num_crews)
	crew_loc_from = Vector{Int}(undef, num_crews)
	crew_rest_from = Vector{Int}(undef, num_crews)

	for crew ∈ 1:num_crews
		crew_model = crew_models[crew]
		long_arcs = crew_model.long_arcs
		first_arcs = findall(
			(long_arcs[:, CM.TIME_FROM] .== 0) .&&
			(long_arcs[:, CM.CREW_NUMBER] .== crew),
		)
		arc = long_arcs[first_arcs[1], :]
		crew_time_from[crew] = arc[CM.TIME_FROM]
		crew_from_type[crew] = arc[CM.FROM_TYPE]
		crew_loc_from[crew] = arc[CM.LOC_FROM]
		crew_rest_from[crew] = arc[CM.REST_FROM]
	end

	# initialize fire states
	# fire_states = Vector{Int}(undef, num_fires)
	# for fire ∈ 1:num_fires
	# 	fire_long_arcs = fire_models[fire].long_arcs
	# 	start_arcs = findall(fire_long_arcs[:, FM.TIME_FROM] .== 0)
	# 	@assert length(unique(fire_long_arcs[start_arcs, FM.STATE_FROM]))[1] == 1
	# 	fire_states[fire] = fire_long_arcs[start_arcs[1], FM.STATE_FROM]
	# end

	# intitialize fire branching rules
	branching_rules = [FireDemandBranchingRule[] for fire ∈ 1:num_fires]

	# myopic assignment can send crews off too far to make their rest
	# so we precalculate the states that should be under consideration

    # TODO FIX THIS
	reachable_states = [backward_crew_dp_subproblem(
			crew_models[crew].wide_arcs,
			crew_models[crew].arc_costs .> 1e3,
			crew_models[crew].state_in_arcs,
		) for crew ∈ 1:num_crews]

    println(reachable_states[9])

	# for each time period
	for curr_time ∈ 0:num_times

        println(crew_time_from)
        println(crew_from_type)
        println(crew_rest_from)
        println(crew_loc_from)

		# update fire branching rules to encode suppression at this time
		if curr_time > 0
			for fire ∈ 1:num_fires
				crews_present = sum([
					1 for i ∈ 1:num_crews if
					(crew_from_type[i] == CM.FIRE_CODE) & (crew_loc_from[i] == fire)
				])
				push!(
					branching_rules[fire],
					FireDemandBranchingRule(
						fire,
						curr_time,
						crews_present,
						"less_than_or_equal",
					),
				)
			end
		end

        println(branching_rules)

		# get the fire triage scores
		no_suppression_fire_damages = [
			no_suppression_fire_damage(
				fire_models[fire],
				num_time_periods,
				branching_rules[fire],
			) for fire ∈ 1:num_fires
		]

		fire_triage_scores =
			no_suppression_fire_damages ./ sum(no_suppression_fire_damages)
		
        println(fire_triage_scores)

		# get the possible assignments

		# formulate assignment problem
		m = direct_model(Gurobi.Optimizer(GRB_ENV))
		set_optimizer_attribute(m, "OutputFlag", 0) # put this first so others don't print
		set_optimizer_attribute(m, "OptimalityTol", 1e-9)
		set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
		fire_totals = @variable(m, s[g = 1:num_fires, t = curr_time+1:num_times+1] >= 0)
		@objective(
			m,
			Max,
			sum(
				λ^(t - curr_time) * fire_totals[g, t] for g ∈ 1:num_fires,
				t ∈ curr_time+1:num_times+1
			)
		)
		arc_variables =
			@variable(m, A[crew ∈ 1:1, loc ∈ 1:1, time ∈ 1:1, rest ∈ 0:1; 0 != 0])
		arc_lookup = Dict()
		triage_rules = @constraint(
			m,
			[
				fire_1 ∈ 1:num_fires,
				fire_2 ∈ 1:num_fires;
				fire_triage_scores[fire_1] > fire_triage_scores[fire_2],
			],
			sum(
				λ^(t - curr_time) * fire_totals[fire_1, t] for
				t ∈ curr_time+1:num_times+1
			) -
			(fire_triage_scores[fire_1] / fire_triage_scores[fire_2])^β * sum(
				λ^(t - curr_time) * fire_totals[fire_2, t] for
				t ∈ curr_time+1:num_times+1
			) >= 0
		)

		linking = @constraint(
			m,
			[g = 1:num_fires, t = curr_time+1:num_times+1],
			fire_totals[g, t] <= 0
		)

		for crew ∈ 1:num_crews
			if crew_time_from[crew] == curr_time
				crew_long_arcs = crew_models[crew].long_arcs
				potential_arc_ixs =
					findall(
                        (crew_long_arcs[:, CM.CREW_NUMBER] .== crew) .&&
						(crew_long_arcs[:, CM.TIME_FROM] .== curr_time) .&&
						(crew_long_arcs[:, CM.LOC_FROM] .== crew_loc_from[crew]) .&&
						(
							crew_long_arcs[:, CM.FROM_TYPE] .== crew_from_type[crew]
						) .&&
						(crew_long_arcs[:, CM.REST_FROM] .== crew_rest_from[crew]),
					)
				for ix ∈ potential_arc_ixs
					loc_to =
						(crew_long_arcs[ix, CM.TO_TYPE] == CM.BASE_CODE) ?
						num_fires + 1 : crew_long_arcs[ix, CM.LOC_TO]
					time_to = crew_long_arcs[ix, CM.TIME_TO]
					rest_to = crew_long_arcs[ix, CM.REST_TO]
					if reachable_states[crew][loc_to, time_to, rest_to + 1]
						arc_variables[crew, loc_to, time_to, rest_to] =
							@variable(m, binary = true)
						arc_lookup[(crew, loc_to, time_to, rest_to)] = ix
						if loc_to < num_fires + 1
							set_normalized_coefficient(
								linking[loc_to, time_to],
								arc_variables[crew, loc_to, time_to, rest_to],
								-1,
							)
						end
					end
				end
                print(crew)
				c = @constraint(m, sum(arc_variables[crew, :, :, :]) == 1)
			end
		end



		# solve assignment problem
		optimize!(m)

		# extract decisions
		chosen_variables = [i for i ∈ eachindex(arc_variables) if value(arc_variables[i]) > 0.5]
		println(chosen_variables)

		# update arcs traversed and crew state
		for (crew, loc_to, time_to, rest_to) ∈ chosen_variables
			crew_long_arcs = crew_models[crew].long_arcs
			chosen_arc = arc_lookup[(crew, loc_to, time_to, rest_to)]
			push!(crew_arcs_traversed[crew], chosen_arc)
			crew_time_from[crew] = crew_long_arcs[chosen_arc, CM.TIME_TO]
			crew_from_type[crew] = crew_long_arcs[chosen_arc, CM.TO_TYPE]
			crew_loc_from[crew] = crew_long_arcs[chosen_arc, CM.LOC_TO]
			crew_rest_from[crew] = crew_long_arcs[chosen_arc, CM.REST_TO]
		end

		# calculate cost from arcs traversed
	end

end

function triage_then_route_full_horizon(
	crew_models::Vector{TimeSpaceNetwork},
	fire_models::Vector{TimeSpaceNetwork},
)

end

num_crews = 10
num_fires = 3
num_time_periods = 14

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
triage_then_route_by_time_period(crew_models, fire_models, 0.9, 1.0)
