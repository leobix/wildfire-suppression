include("CommonStructs.jl")

struct CrewSupplyBranchingRule

	crew_ix::Int64
	fire_ix::Int64
	time_ix::Int64
	visits::Bool
end

struct FireDemandBranchingRule

	fire_ix::Int64
	time_ix::Int64
	allotment::Int64
	branch_direction::String
end

function satisfies_branching_rule(
	b_rule::CrewSupplyBranchingRule,
	fires_fought::BitArray{2},
)

	return fires_fought[b_rule.fire_ix, b_rule.time_ix] ==
		   b_rule.visits
end

function satisfies_branching_rule(
	b_rule::CrewSupplyBranchingRule,
	fire_fought::Bool,
)

	return fire_fought == b_rule.visits
end

function satisfies_branching_rule(
	b_rule::FireDemandBranchingRule,
	crews_present::Vector{Int64},
)
	satisfies::Bool = false
	if b_rule.branch_direction == "less_than_or_equal"
		satisfies = crews_present[b_rule.time_ix] <= b_rule.allotment
	elseif b_rule.branch_direction == "equals"
		satisfies = crews_present[b_rule.time_ix] == b_rule.allotment
	elseif b_rule.branch_direction == "greater_than_or_equal"
		satisfies = crews_present[b_rule.time_ix] >= b_rule.allotment
	else
		error("Branch direction not implemented")
	end

	return satisfies
end

function satisfies_branching_rule(
	b_rule::FireDemandBranchingRule,
	crews_present::Int64,
)
	satisfies::Bool = false
	if b_rule.branch_direction == "less_than_or_equal"
		satisfies = crews_present <= b_rule.allotment
	elseif b_rule.branch_direction == "equals"
		satisfies = crews_present == b_rule.allotment
	elseif b_rule.branch_direction == "greater_than_or_equal"
		satisfies = crews_present >= b_rule.allotment
	else
		error("Branch direction not implemented")
	end

	return satisfies
end

struct GlobalFireAllotmentBranchingRule

	allotment_matrix::Matrix{Int}
	geq_flag::Bool
	fire_sp_arc_lookup::Dict{Int64, Vector{Int64}} # fire sp arcs that exceed allotment
	mp_lookup::Dict{Tuple{Int64, Int64}, Int64} # total amount by which fire plans exceed allotment
end


function add_fire_plan_to_mp_lookup!(
	branching_rule::GlobalFireAllotmentBranchingRule,
	fire::Int64,
	plan_ix::Int64,
	fire_plans::FirePlanData,
)

	coeff = sum(
		branching_rule.allotment_matrix[fire, :] .<
		fire_plans.crews_present[fire, plan_ix, :],
	)
	branching_rule.mp_lookup[(fire, plan_ix)] = coeff

end


function find_fire_sp_arc_lookup(
	allotment_by_time::Vector{Int},
	fire_sp::TimeSpaceNetwork,
)

	TIME_FROM_ = 3
	CREWS_PRESENT_ = 6

	exceed_arcs = Int64[]
	for time_ix ∈ eachindex(allotment_by_time)
		for arc in 1:length(fire_sp.arc_costs)
			if (fire_sp.long_arcs[arc, TIME_FROM_] == time_ix) &
			   (fire_sp.long_arcs[arc, CREWS_PRESENT_] > allotment_by_time[time_ix])
				push!(exceed_arcs, arc)
			end
		end
	end

	return exceed_arcs

end

function GlobalFireAllotmentBranchingRule(
	allotment_matrix::Matrix{Int},
	geq_flag::Bool,
	fire_plans::FirePlanData,
	fire_subproblems::Vector{TimeSpaceNetwork},
)

	fire_sp_arc_lookup = Dict{Int64, Vector{Int64}}()
	mp_lookup = Dict{Tuple{Int64, Int64}, Int64}()

	for fire in eachindex(fire_subproblems)
		fire_allots = allotment_matrix[fire, :]
		fire_sp_arc_lookup[fire] = find_fire_sp_arc_lookup(fire_allots, fire_subproblems[fire])

		for ix in 1:fire_plans.plans_per_fire[fire]
			coeff = sum(fire_allots .< fire_plans.crews_present[fire, ix, :])
			mp_lookup[(fire, ix)] = coeff
		end
	end

	return GlobalFireAllotmentBranchingRule(
		allotment_matrix,
		geq_flag,
		fire_sp_arc_lookup,
		mp_lookup,
	)
end
