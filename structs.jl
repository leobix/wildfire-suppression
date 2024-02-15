using JuMP

mutable struct CrewRouteData

    n_crews::Int64
    routes_per_crew::Vector{Int64} 
    route_costs::Matrix{Float64}
    fires_fought::BitArray{4}

end

mutable struct FirePlanData
    
    n_fires::Int64
    plans_per_fire::Vector{Int64}
    plan_costs::Matrix{Float64}
    crews_present::Array{Int8, 3}
    
end

struct FireTriageAndRouteInstance

    num_crews::Int64
    num_fires::Int64
    num_time_periods::Int64

end


struct TimeSpaceNetwork # always make time the first index

	arc_costs::Vector{Float64}
	state_in_arcs::Matrix{Vector{Int64}}

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

struct RelaxedMasterProblem

    model::JuMP.Model

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