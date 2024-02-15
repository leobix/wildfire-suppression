include("structs.jl")

using Gurobi
# const GRB_ENV = Gurobi.Env()

function branch_and_price()

    

end


## data processing

function preprocess_fire_data()

end

function preprocess_crew_data()

end


### cg

function intialize_column_generation()

end


function master_problem(gurobi_env, crew_route_data, fire_plan_data, dual_warm_start_config=nothing)

    m = Model(() -> Gurobi.Optimizer(gurobi_env))

    set_optimizer_attribute(m, "OptimalityTol", 1e-9)
    set_optimizer_attribute(m, "FeasibilityTol", 1e-9)


    # decision variables
    @variable(m, route[c=1:NUM_CREWS, r=1:crew_route_data.routes_per_crew[c]] >= 0)
    @variable(m, plan[g=1:NUM_FIRES, p=1:fire_plan_data.plans_per_fire[g]] >= 0)

    if dual_warm_start_config != nothing

        # dual stabilization variables
        @variable(m, delta_plus[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS] >= 0)
        @variable(m, delta_minus[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS] >= 0)

    end


    # constraints that you must choose a plan per crew and per fire
    @constraint(m, route_per_crew[c=1:NUM_CREWS],
        sum(route[c, r] for r = 1:route_data.routes_per_crew[c]) == 1)
    @constraint(m, plan_per_fire[g=1:NUM_FIRES],
        sum(plan[g, p] for p = 1:supp_plan_data.plans_per_fire[g]) >= 1)

    # linking constraint


    if dual_warm_start_config == nothing
  
        @constraint(m, linking[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS],

            # crews at fire
            sum(route[c, r] * route_data.fires_fought[c, r, g, t]
                for c = 1:NUM_CREWS, r = 1:route_data.routes_per_crew[c])
            >=

            # crews suppressing
            sum(plan[g, p] * supp_plan_data.crews_present[g, p, t]
                for p = 1:supp_plan_data.plans_per_fire[g]))

    elseif dual_warm_start_config.strategy == "global"

        # get expected dual value ratios
        ratios = dual_warm_start_config.warm_start_values
        ratios = ratios / sum(ratios)

        # get secondary dual objective epsilon
        secondary_eps = dual_warm_start_config.epsilon

        @constraint(m, linking[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS],

            # crews at fire
            sum(route[c, r] * route_data.fires_fought[c, r, g, t]
                for c = 1:NUM_CREWS, r = 1:route_data.routes_per_crew[c])
            +

            # perturbation
            delta_plus[g, t] - delta_minus[g, t] -
            sum(ratios .* delta_plus) + sum(ratios .* delta_minus)
            >=

            # crews suppressing
            sum(plan[g, p] * supp_plan_data.crews_present[g, p, t]
                for p = 1:supp_plan_data.plans_per_fire[g]))
                    
        @constraint(m, perturb[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS],
            delta_plus[g, t] + delta_minus[g, t] <= secondary_eps)
    else
        error("Dual stabilization type not implemented")
    end


    @objective(m, Min,

        # route costs
        sum(route[c, r] * route_data.route_costs[c, r]
            for c = 1:NUM_CREWS, r = 1:route_data.routes_per_crew[c])
        +

        # suppression plan costs
        sum(plan[g, p] * supp_plan_data.plan_costs[g, p]
            for g = 1:NUM_FIRES, p = 1:supp_plan_data.plans_per_fire[g])
    )

    return Dict("m" => m, "sigma" => route_per_crew, "pi" => plan_per_fire,
        "rho" => linking, "route" => route, "plan" => plan)

end

function solve_crew_subproblem(
	ts_network::TimeSpaceNetwork,
	linking_constraint_duals,
	branching_rules,
)

end

function solve_fire_subproblem(
	ts_network::TimeSpaceNetwork,
	linking_constraint_duals,
	branching_rules,
)

end

function choose_natural_variable_for_branching(
	branch_and_bound_node::BranchAndBoundNode,
)

end

function update_master_problem()

end

function double_column_generation(
	relaxed_master_problem::RelaxedMasterProblem,
	crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork},
    crew_branching_rules::Vector{CrewSupplyBranchingRule},
    fire_branching_rules::Vector{FireDemandBranchingRule},
    crew_routes::CrewRouteData,
    fire_plans::FirePlanData)

    # gather global information
    num_crews = 1 #master_problem...;
    num_fires = 1 #master_problem...;
    num_time_periods = 1 #master_problem...;

    # all zeros is a feasible dual solution, initialize there
    fire_duals = zeros(num_fires)
    crew_duals = zeros(num_crews)
    linking_duals = zeros(num_fires, num_time_periods)

    # initialize column generation loop
    terminate::Bool = false
    while ~terminate

        # default is to terminate unless we find a new variable
        terminate = true

        # for each crew

            # generate the local costs of the arcs

            # solve the subproblem

            # if there is an improving route

                # update the master problem

                # add the route to the routes

                terminate = false

        # for each fire

            # generate the local costs of the arcs

            # solve the subproblem

            # if there is an improving plan

                # update the master problem

                # add the plan to the plans

                terminate = false
    end

end


function explore_node(branch_and_bound_node::BranchAndBoundNode,
	all_nodes::Vector{BranchAndBoundNode},
    crew_subproblems::Vector{TimeSpaceNetwork},
	fire_subproblems::Vector{TimeSpaceNetwork})

    # get active branching rules by following tree upward

    # get dual values at solution of parent node for dual warm start

	# find available crew and fire plans based on branching rule to warm-start DCG

	# define master problem
	mp = define_relaxed_master_problem(crew_plans, fire_plans)

	# run DCG, adding columns as needed
	double_column_generation!(mp, crew_subproblems, fire_subproblems)

	# update the branch-and-bound node to be feasible or not

	# decide the branching rule


end




function test_BranchAndBoundNode()

	bb_node = BranchAndBoundNode(ix=1, parent_ix=-1)
    println(bb_node)

end

test_BranchAndBoundNode()