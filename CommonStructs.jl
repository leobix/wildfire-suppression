using JuMP

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

function CrewRouteData(
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
	crews_present::Array{Int64, 3}

end

function FirePlanData(
	max_supp_plans::Int64,
	num_fires::Int64,
	num_time_periods::Int64,
)

	return FirePlanData(num_fires, zeros(num_fires),
		Matrix{Float64}(undef, num_fires, max_supp_plans),
		zeros(Int64, (num_fires, max_supp_plans, num_time_periods)),
	)
end

struct GUBCoverCut
	time_ix::Int64
	fire_lower_bounds::Dict{Int64, Tuple{Int64, Int8}} # fire : (allotment : coeff)
	inactive_crews::Vector{Int64}
	rhs::Int64
end


function GUBCoverCut(
	time_ix::Int64,
	fire_lower_bounds::Dict{Int64, Int64},
	inactive_crews::Vector{Int64},
	rhs::Int64,
)

	return GUBCoverCut(
		time_ix,
		Dict(key => (fire_lower_bounds[key], 1) for key in keys(fire_lower_bounds)),
		inactive_crews,
		rhs,
	)
end

struct GUBCoverCutData

	cut_dict::Dict{Any, GUBCoverCut}
	cuts_per_time::Vector{Int64}
	fire_sp_lookup::Dict{Int64, Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}}
	crew_sp_lookup::Dict{Int64, Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}}
end

# constructor
function GUBCoverCutData(num_crews::Int, num_fires::Int, num_time_periods::Int)

	cut_objects = Dict{Any, GUBCoverCut}()
	cuts_per_time = [0 for t ∈ 1:num_time_periods]
	fire_cut_sp_lookup = Dict{Int64, Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}}()
	crew_cut_sp_lookup = Dict{Int64, Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}}()
	for fire ∈ 1:num_fires
		fire_cut_sp_lookup[fire] = Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}()
	end
	for crew ∈ 1:num_crews
		crew_cut_sp_lookup[crew] = Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}()
	end

	return GUBCoverCutData(
		cut_objects,
		cuts_per_time,
		fire_cut_sp_lookup,
		crew_cut_sp_lookup,
	)
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

mutable struct RestrictedMasterProblem # TODO can make some JuMP things const?

	# model
	model::JuMP.Model

	# vars
	routes::JuMP.Containers.SparseAxisArray # could speed up?
	plans::JuMP.Containers.SparseAxisArray # could speed up?

	# constraints
	route_per_crew::Vector{ConstraintRef}
	plan_per_fire::Vector{ConstraintRef}
	supply_demand_linking::Matrix{ConstraintRef}
	gub_cover_cuts::JuMP.Containers.SparseAxisArray
	# linking_perturbation::Matrix{ConstraintRef}

	# termination status
	termination_status::MOI.TerminationStatusCode

end

mutable struct BranchAndBoundNode

	const ix::Int64
	const parent::Union{Nothing, BranchAndBoundNode}
	const new_crew_branching_rules::Vector{CrewSupplyBranchingRule}
	const new_fire_branching_rules::Vector{FireDemandBranchingRule}
	children::Vector{BranchAndBoundNode}
	l_bound::Float64
	master_problem::Union{Nothing, RestrictedMasterProblem}
	feasible::Union{Nothing, Bool}
	integer::Union{Nothing, Bool}
end

function BranchAndBoundNode(
	;
	ix::Int64,
	parent::Union{Nothing, BranchAndBoundNode},
	new_crew_branching_rules::Vector{CrewSupplyBranchingRule} = CrewSupplyBranchingRule[],
	new_fire_branching_rules::Vector{FireDemandBranchingRule} = FireDemandBranchingRule[],
	children::Vector{BranchAndBoundNode} = BranchAndBoundNode[],
	l_bound::Float64 = -Inf,
	master_problem::Union{Nothing, RestrictedMasterProblem} = nothing,
	feasible::Union{Nothing, Bool} = nothing,
	integer::Union{Nothing, Bool} = nothing)

	return BranchAndBoundNode(
		ix,
		parent,
		new_crew_branching_rules,
		new_fire_branching_rules,
		children,
		l_bound,
		master_problem,
		feasible,
		integer,
	)
end
