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

struct CrewRouteData

	n_crews::Int64
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

struct FirePlanData

	n_fires::Int64
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
	fire_mp_lookup::Dict{Tuple{Int64, Int64}, Dict{Tuple{Int64, Int64}, Int8}}
	crew_mp_lookup::Dict{Tuple{Int64, Int64}, Dict{Tuple{Int64, Int64}, Int8}}
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
		Dict{Tuple{Int64, Int64}, Dict{Tuple{Int64, Int64}, Int64}}(),
		Dict{Tuple{Int64, Int64}, Dict{Tuple{Int64, Int64}, Int64}}(),
	)
end

function restrict_GUBCoverCutData(orig::GUBCoverCutData, ixs)
	cut_objects = Dict(a => b for (a, b) in orig.cut_objects if a in ixs)
	cuts_per_time = orig.cuts_per_time
	fire_sp_lookup = Dict(
		a => Dict(b => c for (b, c) in orig.fire_sp_lookup[a] if b ∈ ixs) for
		a ∈ keys(orig.fire_sp_lookup)
	)
	crew_sp_lookup = Dict(
		a => Dict(b => c for (b, c) in orig.crew_sp_lookup[a] if b ∈ ixs) for
		a ∈ keys(orig.crew_sp_lookup)
	)
	fire_mp_lookup = Dict(a => b for (a, b) in orig.fire_mp_lookup if a ∈ ixs)
	crew_mp_lookup = Dict(a => b for (a, b) in orig.crew_mp_lookup if a ∈ ixs)
	return GUBCoverCutData(
		cut_objects,
		cuts_per_time,
		fire_cut_sp_lookup,
		crew_cut_sp_lookup,
		fire_mp_lookup,
		crew_mp_lookup)
end

mutable struct DualRegularizer

	const strategy::String
	epsilon::Float64
	constraint_ref::Union{Nothing, ConstraintRef}

end

mutable struct RestrictedMasterProblem # TODO can make some JuMP things const?

	# model
	const model::JuMP.Model

	# vars
	const routes::JuMP.Containers.SparseAxisArray # could speed up?
	const plans::JuMP.Containers.SparseAxisArray # could speed up?

	# constraints
	const route_per_crew::Vector{ConstraintRef}
	const plan_per_fire::Vector{ConstraintRef}
	const supply_demand_linking::Matrix{ConstraintRef}
	const gub_cover_cuts::JuMP.Containers.SparseAxisArray
	const fire_allotment_branches::Vector{ConstraintRef}
	const dual_regularizer::Union{Nothing, DualRegularizer}

	# termination status
	termination_status::MOI.TerminationStatusCode

end
