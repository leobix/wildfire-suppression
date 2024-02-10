include("DCG.jl")

function preprocess(in_path)

    # get inital fire perimeters and no-suppression progression parameters
    M = readdlm(in_path * "/sample_growth_patterns.csv", ',')
    start_perims = M[:, 1]
    progressions = M[:, 2:15]
    global NUM_TIME_PERIODS = size(M)[2] - 1

    g_data, crew_status = load_data(in_path)
    r_data = RegionData(convert.(Int, ones(NUM_CREWS)), convert.(Int, ones(100)))

    rotation_order = get_rotation_orders(r_data.crew_regions)
    A = generate_arcs(g_data, r_data, crew_status)

    rest_pen = get_rest_penalties(crew_status.rest_by, 1e10, positive)
    cost_params = Dict("cost_per_mile" => 1, "rest_violation" => rest_pen, "fight_fire" => ALPHA)
    arc_costs = get_static_arc_costs(g_data, A, cost_params)

    c_data = define_network_constraint_data(A)

    fire_configs = []
    for fire in 1:NUM_FIRES
        model_config = Dict("model_type" => "simple_linear", "progressions" => progressions[fire, :],
            "start_perim" => start_perims[fire], "line_per_crew" => LINE_PER_CREW,
            "beta" => BETA)
        push!(fire_configs, model_config)
    end

    return A, arc_costs, r_data, c_data, g_data, rotation_order, fire_configs
end

function single_DCG_node(test_features, data)

    ts = Dict{String,Float64}()

    # read in data 
    A, arc_costs, r_data, c_data, g_data, rotation_order, fire_configs = data


    rhos = []
    sigmas = []
    pis = []
    allotment_history = []
    reduced_costs = []

    global AGG_PREC = test_features["agg_prec"]
    global PASSIVE_STATES = test_features["passive_states"]

    gamma = test_features["gamma"]


    best_sols = Dict{String,Float64}()
    allotments = Dict{String,Any}()

    fire_solver_configs = [Dict{String,Any}("solver_type" => test_features["fire_solver_type"]) for fire in 1:NUM_FIRES]
    crew_solver_configs = [Dict{String,Any}("solver_type" => test_features["crew_solver_type"]) for crew in 1:NUM_CREWS]

    max_plans = 1000
    ts["init_cg"] = @elapsed col_gen_data = initialize_column_generation(A, arc_costs, c_data, fire_configs,
        fire_solver_configs, crew_solver_configs, max_plans)


    if (test_features["apply_warm_start"]) | (test_features["dual_warm_start"] == "calculate")

        global AGG_PREC = 30
        global PASSIVE_STATES = 5

        ts["ws"] = @elapsed ws, ws_dict = warm_start_suppression_plans(10, fire_configs, r_data, c_data, rotation_order, arc_costs,
            gamma, false, false, 3600)


        if test_features["dual_warm_start"] == "calculate"

            # recreate linking constraints in full network flow formulation
            test_features["dual_warm_start"] = dual.(ws_dict["f_linking"])
            # test_features["dual_warm_start"] += dual.(ws_dict["r_linking"])
            # state_bounds = vcat([0], ws_dict["fire_state_lookup"])
            # test_features["dual_warm_start"] += dropdims(sum(dual.(ws_dict["crew_flow"]), dims=(1, 4)), dims=(1, 4))
            # fire_flow_duals = reduce(hcat, [dropdims(sum(dual.(ws_dict["fire_flow"])[state_bounds[g]+1:state_bounds[g+1], :], dims=1), dims=1) for g in 1:NUM_FIRES])'
            # test_features["dual_warm_start"] += fire_flow_duals

            test_features["warm_start_allotments"] = value.(ws_dict["d"])
        end

        # if test_features["apply_warm_start"]
        #     for fire in 1:NUM_FIRES
        #         ws_fire = ws[fire]
        #         for plan in ws_fire
        #             cost = get_fire_cost(plan, fire_configs[fire])
        #             update_available_supp_plans(fire, cost, plan, col_gen_data.suppression_plans)
        #         end
        #     end
        # end

        global AGG_PREC = test_features["agg_prec"]
        global PASSIVE_STATES = test_features["passive_states"]

    end

    objs = []
    max_iters = test_features["max_iters"]
    n_iters = 0
    opt = false
    ts["cg"] = 0
    max_time = test_features["cg_time_limit"]
    col_gen_config = Dict{String,Any}("ws_dual" => false)
    col_gen_config["master_problem"] = Dict{String,Any}()

    col_gen_config["master_problem"]["dual_stabilization"] = test_features["dual_stab_type"]
    col_gen_config["master_problem"]["dual_secondary_eps"] = test_features["dual_stab_eps"]
    col_gen_config["master_problem"]["dual_warm_start"] = test_features["dual_warm_start"]

    if test_features["dual_stab_anchor"] == "all_ones"
        col_gen_config["master_problem"]["dual_warm_start"] = ones(size(col_gen_config["master_problem"]["dual_warm_start"]))
    end
    mp = master_problem(col_gen_config, col_gen_data.routes, col_gen_data.suppression_plans,
        r_data, rotation_order, 0, false)

    current_num_routes = copy(col_gen_data.routes.routes_per_crew)
    current_num_plans = copy(col_gen_data.suppression_plans.plans_per_fire)

    while (~opt) & (n_iters < max_iters) & (ts["cg"] < max_time)

        n_iters += 1

        for fire in 1:NUM_FIRES
            fire_solver_configs[fire]["int_aware_adjustment_pattern"] = test_features["int_aware_price_adjustment"][n_iters]
            if "int_aware_capacities" in keys(test_features)
                fire_solver_configs[fire]["int_aware_capacities"] = test_features["int_aware_capacities"]
            end
        end

        t = @elapsed mp, a, r_costs, r, s, p, c =  run_CG_step(col_gen_data, collect(A'), arc_costs, g_data, r_data, fire_configs,
            fire_solver_configs, crew_solver_configs, col_gen_config, rotation_order, gamma,
            test_features["restore_cost"], mp)
        println(t)
        println()
        ts["cg"] += t

        push!(objs, a)
        push!(reduced_costs, r_costs)
        push!(rhos, r)
        push!(sigmas, s)
        push!(pis, p)
        push!(allotment_history, c)

        next_num_routes = col_gen_data.routes.routes_per_crew
        next_num_plans = col_gen_data.suppression_plans.plans_per_fire

        if (sum(next_num_routes) == sum(current_num_routes)) & (sum(next_num_plans) == sum(current_num_plans))
            opt = true
        end

        current_num_routes = copy(next_num_routes)
        current_num_plans = copy(next_num_plans)

    end

    # no more stabilization for assessment
    col_gen_config["master_problem"]["dual_stabilization"] = false

    # for new JuMP, re-optimize before querying
    optimize!(mp["m"])
    best_sols["master_problem"] = objective_value(mp["m"])
    allotments["master_problem"] = get_fire_allotments(mp, col_gen_data)

    ts["mp_reformulate"] = @elapsed mp = master_problem(col_gen_config, col_gen_data.routes, col_gen_data.suppression_plans,
        r_data, rotation_order, gamma, false)
    ts["mp_resolve"] = @elapsed optimize!(mp["m"])
    best_sols["master_problem_reconstructed"] = objective_value(mp["m"])
    allotments["master_problem_reconstructed"] = get_fire_allotments(mp, col_gen_data)

    if test_features["solve_pb_time_limit"] > 0

        ts["price_and_branch_formulate"] = @elapsed pb = master_problem(col_gen_config,
            col_gen_data.routes, col_gen_data.suppression_plans, r_data, rotation_order, gamma, true)
        set_optimizer_attribute(pb["m"], "TimeLimit", test_features["solve_pb_time_limit"])
        ts["price_and_branch_solve"] = @elapsed optimize!(pb["m"])

        best_sols["price_and_branch"] = objective_value(pb["m"])

        plans_chosen = [i for i in eachindex(pb["plan"]) if value(pb["plan"][i]) > 0.5]

        pb_allotment = convert.(Int8, zeros(size(col_gen_data.suppression_plans.crews_present[:, 1, :])))
        for plan in plans_chosen
            pb_allotment[plan[1], :] = col_gen_data.suppression_plans.crews_present[plan[1], plan[2], :]
        end

        allotments["price_and_branch"] = pb_allotment
    end

    if test_features["solve_explicit_int_time_limit"] > 0

        tl = test_features["solve_explicit_int_time_limit"]
        m, p, l, z, q, oor, linking = full_formulation(true, r_data, c_data, rotation_order, arc_costs,
            progressions, start_perims, BETA, gamma, false, tl)
        ts["explicit_int"] = @elapsed optimize!(m)

        if has_values(m)
            best_sols["explicit_int"] = objective_value(m)
            allotments["explicit_int"] = convert.(Int8, ceil.(value.(l) / LINE_PER_CREW .- 0.001))
        else
            best_sols["explicit_int"] = false
        end

    end


    if test_features["solve_explicit_lin_time_limit"] > 0

        tl = test_features["solve_explicit_lin_time_limit"]

        m2, p, l, z, q, oor, linking = full_formulation(false, r_data, c_data, rotation_order, arc_costs,
            progressions, start_perims, BETA, gamma, false, tl)
        ts["explicit_lr"] = @elapsed optimize!(m2)

        if has_values(m2)
            best_sols["explicit_lr"] = objective_value(m2)
            allotments["explicit_lr"] = value.(l) / LINE_PER_CREW
        else
            best_sols["explicit_lr"] = false
        end

    end

    if test_features["solve_net_flow_time_limit"] > 0
        tl = test_features["solve_net_flow_time_limit"]
        ts["net_flow_int"] = @elapsed ws, ws_dict = warm_start_suppression_plans(10, fire_configs, r_data, c_data, rotation_order, arc_costs, gamma, true, false, tl)
        if has_values(ws_dict["m"])
            best_sols["net_flow_int"] = objective_value(ws_dict["m"])
            best_sols["net_flow_int_bound"] = objective_bound(ws_dict["m"])
            allotments["net_flow_int"] = value.(ws_dict["d"])
            allotments["net_flow_int_full_data"] = ws_dict
        else
            best_sols["net_flow_int"] = false
        end
    end

    d = Dict{String,Any}()
    d["iterations"] = n_iters
    d["times"] = ts
    d["solutions"] = best_sols
    d["rhos"] = rhos
    d["objs"] = objs
    d["sigmas"] = sigmas
    d["pis"] = pis
    d["mp"] = mp
    d["allotment_history"] = allotment_history
    allotments["warm_start"] = test_features["warm_start_allotments"]
    d["allotments"] = allotments
    d["reduced_costs"] = reduced_costs
    d["dual_warm_start"] = test_features["dual_warm_start"]
    return d, col_gen_data
end

function default_params()

    params = Dict{String,Any}()
    params["gamma"] = 0

    params["in_path"] = "data/raw/big_fire"
    params["agg_prec"] = 10
    params["passive_states"] = 30
    params["num_fires"] = 6
    params["num_crews"] = 20
    params["line_per_crew"] = 17
    params["fire_solver_type"] = "dp_fast"
    params["crew_solver_type"] = "dp_fast"

    params["apply_warm_start"] = false
    params["restore_cost"] = false

    params["solve_pb_time_limit"] = 0
    params["solve_explicit_lin_time_limit"] = 0
    params["solve_explicit_int_time_limit"] = 0
    params["cg_time_limit"] = 300
    params["solve_net_flow_time_limit"] = 0

    params["max_iters"] = 1000

    params["dual_stab_type"] = "global"
    params["dual_stab_anchor"] = false
    params["dual_stab_eps"] = 0.01
    params["dual_warm_start"] = "calculate"

    pattern = [(0, 0)]
    params["int_aware_price_adjustment"] = [pattern for i in 1:params["max_iters"]]

    return params
end

function generate_new_plans(cg_params, prepped_data, cg_data, excesses, crew_thresh, fire_thresh)

    iters = []
    times = []
    existing_plans = [[cg_data.suppression_plans.crews_present[i, j, :]
                       for j = 1:cg_data.suppression_plans.plans_per_fire[i]
    ]
                      for i = 1:NUM_FIRES
    ]

    existing_routes = [[cg_data.routes.fires_fought[i, j, :, :]
                        for j = 1:cg_data.routes.routes_per_crew[i]
    ]
                       for i = 1:NUM_CREWS
    ]

    for excess in excesses

        pattern = [(excess, 1e30)]
        cg_params["int_aware_price_adjustment"] = [pattern for i in 1:cg_params["max_iters"]]

        t = @elapsed d2, cg_data2 = single_DCG_node(cg_params, deepcopy(prepped_data))
        push!(iters, d2["iterations"])
        push!(times, t)


        for fire in 1:NUM_FIRES

            for plan_ix in 1:cg_data2.suppression_plans.plans_per_fire[fire]

                if value(d2["mp"]["plan"][fire, plan_ix]) > fire_thresh

                    plan = cg_data2.suppression_plans.crews_present[fire, plan_ix, :]
                    cost = cg_data2.suppression_plans.plan_costs[fire, plan_ix]

                    if ~(plan in existing_plans[fire])

                        push!(existing_plans[fire], plan)
                        update_available_supp_plans(fire, cost, plan, cg_data.suppression_plans)

                    end
                end
            end
        end


        for crew in 1:NUM_CREWS

            for route_ix in 1:cg_data2.routes.routes_per_crew[crew]

                if value(d2["mp"]["route"][crew, route_ix]) > crew_thresh

                    route = cg_data2.routes.fires_fought[crew, route_ix, :, :]
                    cost = cg_data2.routes.route_costs[crew, route_ix]
                    oo_region = cg_data2.routes.out_of_reg[crew, route_ix, :]

                    if ~(route in existing_routes[crew])

                        push!(existing_routes[crew], route)
                        update_available_crew_routes(crew, cost, route, oo_region, cg_data.routes)
                    end

                end
            end
        end
    end

    return iters, times, cg_data
end

function restore_integrality(cg_data, time_limit)

    config = Dict{String,Any}("ws_dual" => false)
    config["master_problem"] = Dict{String,Any}()

    config["master_problem"]["dual_stabilization"] = false
    config["master_problem"]["dual_secondary_eps"] = false
    config["master_problem"]["dual_warm_start"] = false

    r_data = RegionData(convert.(Int, ones(NUM_CREWS)), convert.(Int, ones(100)))
    rotation_order = get_rotation_orders(r_data.crew_regions)

    form_time = @elapsed pb = master_problem(config, cg_data.routes, cg_data.suppression_plans,
        r_data, rotation_order, 0, true)
    set_optimizer_attribute(pb["m"], "TimeLimit", max(1, time_limit - form_time))
    set_optimizer_attribute(pb["m"], "OutputFlag", 0)
    sol_time = @elapsed optimize!(pb["m"])

    return form_time, sol_time, pb
end
