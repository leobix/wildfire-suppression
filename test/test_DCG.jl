include("../BranchAndPrice.jl")
include("../DoubleColumnGeneration.jl")
include("../GUBKnapsackCoverCuts.jl")

using JuMP, Gurobi, Profile
const GRB_ENV = Gurobi.Env()

function initialize_data_structures(
	num_fires::Int64,
	num_crews::Int64,
	num_time_periods::Int64,
)
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


	crew_routes = CrewRouteData(100000, num_fires, num_crews, num_time_periods)
	fire_plans = FirePlanData(100000, num_fires, num_time_periods)
	cut_data = CutData(num_crews, num_fires, num_time_periods)

	rmp = define_restricted_master_problem(
		GRB_ENV,
		crew_routes,
		[Int64[] for i ∈ 1:num_crews],
		fire_plans,
		[Int64[] for i ∈ 1:num_fires],
		cut_data
	)
	return crew_routes, fire_plans, crew_models, fire_models, rmp, cut_data

end


s = @elapsed crew_routes, fire_plans, crew_models, fire_models, rmp, cut_data =
    initialize_data_structures(3, 10, 14)


t = @elapsed double_column_generation!(
    rmp,
    crew_models,
    fire_models,
    CrewAssignmentBranchingRule[],
    FireDemandBranchingRule[],
    crew_routes,
    fire_plans,
	cut_data
)

println(s)
println(t)

s = @elapsed crew_routes, fire_plans, crew_models, fire_models, rmp, cut_data =
    initialize_data_structures(3, 10, 14)

Profile.init()
t = @profile @elapsed double_column_generation!(
    rmp,
    crew_models,
    fire_models,
    CrewAssignmentBranchingRule[],
    FireDemandBranchingRule[],
    crew_routes,
    fire_plans,
	cut_data
)

println(s)
println(t)
@assert abs(objective_value(rmp.model) - 1.1388005347712322e6) <= 1e-6
println("Matched expected objective value")
Profile.print()