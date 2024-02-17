struct TimeSpaceNetwork # TODO always make time the first index

	arc_costs::Vector{Float64}
	state_in_arcs::Array{Vector{Int64}}
	model_type::String

	# not quite sure how much this helps performance, but storing 2 copies
	# of array data allows always column-wise access
	long_arcs::Matrix{Int64}
	wide_arcs::Matrix{Int64}

end

mutable struct CrewRouteData

	const n_crews::Int64
	routes_per_crew::Vector{Int64}
	route_costs::Matrix{Float64}
	fires_fought::BitArray{4}

end

function CrewRouteData_init(
	max_routes::Int64,
	num_fires::Int64,
	num_crews::Int64,
	num_time_periods::Int64,
)

	return CrewRouteData(num_crews, zeros(num_crews),
		Matrix{Float64}(undef, num_crews, max_routes),
		BitArray(undef, num_crews, max_routes, num_fires, num_time_periods) .> 2)
end

mutable struct FirePlanData

	const n_fires::Int64
	plans_per_fire::Vector{Int64}
	plan_costs::Matrix{Float64}
	crews_present::Array{Int16, 3}

end

function FirePlanData_init(
	max_supp_plans::Int64,
	num_fires::Int64,
	num_time_periods::Int64,
)

	return FirePlanData(num_fires, zeros(num_fires),
		Matrix{Float64}(undef, num_fires, max_supp_plans),
		zeros(Int16, (num_fires, max_supp_plans, num_time_periods)),
	)
end

struct FireTriageAndRouteInstance

	num_crews::Int64
	num_fires::Int64
	num_time_periods::Int64

end


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


@kwdef mutable struct BranchAndBoundNode

	const ix::Int64
	const parent_ix::Int64
	child_ixs::Vector{Int64} = []
	l_bound::Float64 = Inf
	feasible::Union{Nothing, Bool} = nothing
	integer::Union{Nothing, Bool} = nothing
	fire_plan_ixs_allowed::Vector{Int64} = []
	crew_plan_ixs_allowed::Vector{Int64} = []
end
