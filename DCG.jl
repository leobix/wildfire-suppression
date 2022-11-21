include("CGStructs.jl")

using DataFrames, CSV, DelimitedFiles, JuMP, Gurobi


NUM_CREWS = 10                
BREAK_LENGTH = 2       # how long at base to be considered "rested"

# tradeoffs
BETA = 100             # cost of one area unit burned / cost of mile traveled
ALPHA = 200            # cost of crew-day of suppression / cost of mile traveled
LINE_PER_CREW = 17     # how much perimeter prevented per crew per time period

FIRE_CODE = 1
BASE_CODE = 2

AGG_PREC = 0
PASSIVE_STATES = 0

function get_rotation_orders(crew_regions)
    
    # initialize output
    out = Dict()
    
    # get the unique regions where there are crews
    regions = unique(crew_regions)
    
    # for each region
    for region in regions
        
        # initialize dictionary corresponding to the order
        out[region] = Dict() 
        crews_in_region = 0
        
        # for each crew in the region
        for crew in 1:NUM_CREWS
            
            if crew_regions[crew] == region
                
                # update crew count, log rotation order 
                crews_in_region += 1
                out[region][crew] = crews_in_region
            end
        end
    end
    
    return out
end


function arc_exits_region(crew, from_type, from_ix, to_type, to_ix, region_data)
    
    # get the region where the arc originates
    from_region = 0
    if from_type == FIRE_CODE
        from_region = region_data.fire_regions[from_ix]
    elseif from_type == BASE_CODE
        from_region = region_data.crew_regions[from_ix]
    else
        throw(DomainError(from_type, "from_type invalid"))
    end
    
    # get the region where the arc terminates
    to_region = 0
    if to_type == FIRE_CODE
        to_region = region_data.fire_regions[to_ix]
    elseif to_type == BASE_CODE
        to_region = region_data.crew_regions[to_ix]
    else
        throw(DomainError(from_type, "to_type invalid"))
    end
    
    # if these are different regions
    if from_region != to_region
        
        # if the crew is leaving its home region
        if region_data.crew_regions[crew] == from_region
        
            # return the region that the arc exited
            return from_region
        
        end
        
    end
    
    # otherwise
    return 0
    
end    

# crew, from_type, from_ix, to_type, to_ix, from_time, to_time, from_rested, to_rested, exited_region

function generate_arcs(gd, rd, cs)
    
    # get fire-to-fire arcs
    ff = [[c, FIRE_CODE, f_from, FIRE_CODE, f_to, t_from, t_from + gd.ff_tau[f_to, f_from], rest, rest]
          for c=1:NUM_CREWS, f_from=1:NUM_FIRES, f_to=1:NUM_FIRES, t_from=1:NUM_TIME_PERIODS, rest=0:1]
    ff = copy(reduce(hcat, ff)')

    # get fire-to-fire arcs from start, based on cs.current crew locations
    from_start_ff = [[c, FIRE_CODE, cs.current_fire[c], FIRE_CODE, f_to, 0, gd.ff_tau[f_to, cs.current_fire[c]], 0, 0]
                      for c=1:NUM_CREWS, f_to=1:NUM_FIRES if cs.current_fire[c] != -1]
    from_start_ff = copy(reduce(hcat, from_start_ff)')

    # get base-to-fire arcs
    rf = [[c, BASE_CODE, c, FIRE_CODE, f_to, t_from, t_from + gd.bf_tau[c, f_to], rest, rest]
           for c=1:NUM_CREWS, f_to=1:NUM_FIRES, t_from=1:NUM_TIME_PERIODS, rest=0:1]
    rf = copy(reduce(hcat, rf)')

    # get base-to-fire arcs from start
    from_start_rf = [[c, BASE_CODE, c, FIRE_CODE, f_to, 0, gd.bf_tau[c, f_to], 0, 0]
                      for c=1:NUM_CREWS, f_to=1:NUM_FIRES if cs.current_fire[c] == -1]
    from_start_rf = copy(reduce(hcat, from_start_rf)')

    # get fire-to-base arcs
    fr = [[c, FIRE_CODE, f_from, BASE_CODE, c, t_from, t_from + gd.bf_tau[c, f_from], rest, rest]
           for c=1:NUM_CREWS, f_from=1:NUM_FIRES, t_from=1:NUM_TIME_PERIODS, rest=0:1]
    fr = copy(reduce(hcat, fr)')

    # get fire-to-base arcs from start, based on cs.current crew locations
    from_start_fr = [[c, FIRE_CODE, cs.current_fire[c], BASE_CODE, c, 0, gd.bf_tau[c, cs.current_fire[c]], 0, 0]
                      for c=1:NUM_CREWS if cs.current_fire[c] != -1]
    from_start_fr = copy(reduce(hcat, from_start_fr)')

    # get base-to-base arcs
    rr = [[c, BASE_CODE, c, BASE_CODE, c, t_from, t_from + 1 + (BREAK_LENGTH - 1) * rest, 0, rest]
          for c=1:NUM_CREWS, t_from=1:NUM_TIME_PERIODS, rest=0:1]
    rr = copy(reduce(hcat, rr)')
    rr_rested = [[c, BASE_CODE, c, BASE_CODE, c, t_from, t_from + 1, 1, 1]
          for c=1:NUM_CREWS, t_from=1:NUM_TIME_PERIODS]
    rr_rested  = copy(reduce(hcat, rr_rested)')

    # get base-to-base arcs from start, based on cs.current days rested
    from_start_rr = [[c, BASE_CODE, c, BASE_CODE, c, 0, 
                      1 + (BREAK_LENGTH - max(cs.rested_periods[c], 0) - 1) * rest, 0, rest] 
                      for c=1:NUM_CREWS, rest=0:1 if cs.current_fire[c] == -1]
    from_start_rr = copy(reduce(hcat, from_start_rr)')

    A = vcat(ff, from_start_ff, rf, from_start_rf, fr, from_start_fr, rr, rr_rested, from_start_rr)

    out_of_region = [arc_exits_region(A[i, 1], A[i, 2], A[i, 3], A[i, 4], A[i, 5], rd) 
                     for i in 1:length(A[:, 1])]
    A = hcat(A, out_of_region)
    
    return A
end

function get_distance(from_type, from_ix, to_type, to_ix, fire_fire, base_fire)
    
    dist = 0
    
    # if fire to fire
    if (from_type == FIRE_CODE) & (to_type == FIRE_CODE)
        dist = fire_fire[from_ix, to_ix]
    
    # if fire to base
    elseif (from_type == FIRE_CODE) & (to_type == BASE_CODE)
        dist = base_fire[to_ix, from_ix]
    
    # if base to fire
    elseif (from_type == BASE_CODE) & (to_type == FIRE_CODE)
        dist = base_fire[from_ix, to_ix]
        
    # otherwise dist still 0
    end
    
    return dist
end 

function get_arc_costs(gd, arcs, cost_param_dict)
    
    # get number of arcs
    n_arcs = length(arcs[:, 1])
    
    # initialize costs to 0
    costs = zeros(n_arcs)
    
    # if there is travel cost per mile
    if "cost_per_mile" in keys(cost_param_dict)
        
        # find the miles for each arc
        miles_per_arc =  [get_distance(arcs[i, 2], arcs[i, 3], 
                                       arcs[i, 4], arcs[i, 5], 
                                       gd.ff_dist, gd.bf_dist) for i in 1:n_arcs]
        # add to costs
        costs = costs .+ (cost_param_dict["cost_per_mile"] * miles_per_arc)
    end
    
    # if there are rest violations
    if "rest_violation" in keys(cost_param_dict)
        
        # find the rest violation scores
        rest_violation_matrix = cost_param_dict["rest_violation"]
        rest_violations = [(arcs[i, 8] == 0) & (arcs[i, 6] > 0) ? 
                           rest_violation_matrix[arcs[i, 1], arcs[i, 6]] : 0
                           for i in 1:n_arcs]
        
        # add to costs
        costs = costs .+ rest_violations
    end
    
    if "fight_fire" in keys(cost_param_dict)
        costs = costs .+ [(arcs[i, 4] == FIRE_CODE) ? cost_param_dict["fight_fire"] : 0
                          for i in 1:n_arcs]
    end
    
    # if we have to adjust for linking dual constraints
    if "linking_dual" in keys(cost_param_dict)
        
        # get the dual variables
        rho = cost_param_dict["linking_dual"]
        
        # get linking costs (really benefits) if arc goes to a fire
        linking_costs = [((arcs[i, 4] == FIRE_CODE) & (arcs[i, 7] <= NUM_TIME_PERIODS)) ? 
                          - rho[arcs[i, 5], arcs[i, 7]] : 0
                          for i in 1:n_arcs]
        
        # add to costs
        costs = costs .+ linking_costs
        
    end
    
    # if we have to adjust for linking dual constraints
    if "out_of_region_dual" in keys(cost_param_dict)
        
        # get needed regional info
        regs = cost_param_dict["region_data"].crew_regions
        rot_order = cost_param_dict["rotation_order"]
        
        # get the dual variables
        eta = cost_param_dict["out_of_region_dual"]

        # get adjustment for crew allotment
        c1 = [(arcs[i, 10] > 0) ? sum(eta[arcs[i, 1], t_0]
                                        for t_0=arcs[i, 6]:NUM_TIME_PERIODS
                                      ) : 0
                                                   
               for i in 1:n_arcs
             ]
        
        # get adjustment for region average allotment
        c2 = [(arcs[i, 10] > 0) ? sum(eta[c, t_0]
                                            for c in keys(rot_order[regs[arcs[i, 1]]]),
                                                t_0=arcs[i, 6]:NUM_TIME_PERIODS) /
                                        length(keys(rot_order[regs[arcs[i, 1]]])) : 0
                                                   
               for i in 1:n_arcs
             ]
        
        # get adjustment for big-M constraint
        c3 = [(arcs[i, 10] > 0) ? NUM_TIME_PERIODS * eta[arcs[i, 1], arcs[i, 6]] : 0
               for i in 1:n_arcs
             ]
            
        # add to costs
        costs = costs .+ c1 .- c2 .+ c3
        
    end   
    
    return costs
end

function positive(x)
    
    if x > 0
        return 1
    end
    
    return 0
end

function is_one(x)
    
    if x == 1
        return 1
    end
    
    return 0
end

# should return matrix indexed by crew, time, 
function get_rest_penalties(rest_by_periods, lambda, accounting_func)
    
    penalties = zeros(NUM_CREWS, NUM_TIME_PERIODS)
    
    for c in 1:NUM_CREWS
        penalties[c, :] = [lambda * accounting_func(t - rest_by_periods[c]) 
                           for t in 1:NUM_TIME_PERIODS]
    end
    
    return penalties    
end

function define_network_constraint_data(arcs)
    
    # shorten some global variable names
    C = NUM_CREWS
    G = NUM_FIRES
    T = NUM_TIME_PERIODS
    
    # get number of arcs
    n_arcs = length(arcs[:, 1])
      
    ## flow balance ##
    
    # initialize arrays of vectors for flow balance
    f_out = Array{Vector{Int64}}(undef, C, G, T, 2)
    f_in = Array{Vector{Int64}}(undef, C, G, T, 2)
    b_out = Array{Vector{Int64}}(undef, C, T, 2)
    b_in = Array{Vector{Int64}}(undef, C, T, 2)
    start = Array{Vector{Int64}}(undef, C)
    out_of_region = Array{Vector{Int64}}(undef, C, T+1)
    
    # for each crew
    for crew in 1:C
        
        # get indices of this crew's arcs only
        crew_ixs = [i for i in 1:n_arcs if arcs[i, 1] == crew]
        
        # get time 0 indices
        start[crew] = [i for i in crew_ixs if arcs[i, 6] == 0]
        
        # for each time period (including start)
        for tm in 0:T
        
            # get indices for out of region assignments
            out_of_region[crew, tm+1] = [i for i in crew_ixs if
                                           (arcs[i, 6] == tm) &
                                           (arcs[i, 10] > 0)
                                        ]
        end
        
        # for each time period
        for tm in 1:T
            
            # for each rest state
            for rest in 1:2
                
                # get arcs leaving crew base at this time with this rest
                b_out[crew, tm, rest] = [i for i in crew_ixs if
                                         (arcs[i, 2] == BASE_CODE) &
                                         (arcs[i, 6] == tm) &
                                         (arcs[i, 8] == rest-1)
                                        ]
                
                # get arcs entering crew base at this time with this rest
                b_in[crew, tm, rest] = [i for i in crew_ixs if
                                        (arcs[i, 4] == BASE_CODE) &
                                        (arcs[i, 7] == tm) &
                                        (arcs[i, 9] == rest-1)
                                       ]
                # for each fire
                for fire in 1:G
                    
                    # get arcs where this crew leaves this fire at this time
                    # with this rest state
                    f_out[crew, fire, tm, rest] = [i for i in crew_ixs if
                                                   (arcs[i, 2] == FIRE_CODE) &
                                                   (arcs[i, 3] == fire) &
                                                   (arcs[i, 6] == tm) &
                                                   (arcs[i, 8] == rest-1)
                                                   ]
                    
                    # get arcs where this crew enters this fire at this time
                    # with this rest state
                    f_in[crew, fire, tm, rest] = [i for i in crew_ixs if
                                                  (arcs[i, 4] == FIRE_CODE) &
                                                  (arcs[i, 5] == fire) &
                                                  (arcs[i, 7] == tm) &
                                                  (arcs[i, 9] == rest-1)
                                                  ]
                end
            end
        end
    end
    
    ## linking constraints ##
    linking = Array{Vector{Int64}}(undef, G, T)
    for fire in 1:G
        for tm in 1:T
            
            # we count the crew as working *where they arrived* during this timestep
            linking[fire, tm] = [i for i in 1:n_arcs if (arcs[i, 4] == FIRE_CODE) &
                                                        (arcs[i, 5] == fire) &
                                                        (arcs[i, 7] == tm)]
        end
    end
    
    
    return KeyArcIndices(f_out, f_in, b_out, b_in, linking, start, out_of_region)
end

function get_route_stats(arc_ixs_used, arcs, costs)
    
    # get total cost
    route_cost = sum(costs[arc_ixs_used])
    
    # initialize fires fought matrix
    fires_fought =  falses(NUM_FIRES, NUM_TIME_PERIODS)
    
    # initialize out of region matrix
    out_of_region = falses(NUM_TIME_PERIODS + 1)
    
    # for each arc used
    for ix in arc_ixs_used
        arc = arcs[ix, :]
        
        # update fires_fought
        if (arc[4] == FIRE_CODE) & (arc[7] <= NUM_TIME_PERIODS)
            @assert ~fires_fought[arc[5], arc[7]] "Visited fire twice at same time"
            fires_fought[arc[5], arc[7]] = true
        end
        
        # update out_of_region
        if arc[10] > 0
            @assert ~out_of_region[arc[6] + 1] "Left region twice at same time"
            out_of_region[arc[6] + 1] = true
        end
    end
    
    return route_cost, fires_fought, out_of_region
end

function initialize_route_data(max_routes)
    
    return RouteData(zeros(NUM_CREWS), Matrix{Float64}(undef, NUM_CREWS, max_routes),
                     BitArray(undef, NUM_CREWS, max_routes, NUM_FIRES, NUM_TIME_PERIODS) .> 2,
                     BitArray(undef, NUM_CREWS, max_routes, NUM_TIME_PERIODS + 1) .> 2)
end

function update_available_routes(crew, route_ixs, arcs, costs, route_data)
    
    # get the required information from the arcs used
    route_cost, fires_fought, out_of_region = get_route_stats(route_ixs, arcs, costs)
    
    ## store this information to the route_data ##
    return update_available_crew_routes(crew, route_cost, fires_fought, out_of_region, route_data)

end

function update_available_crew_routes(crew, cost, fires_fought, oor, route_data)
    
    # add 1 to number of routes for this crew, store the index
    route_data.routes_per_crew[crew] += 1
    ix = route_data.routes_per_crew[crew]
    
    # append the route cost
    route_data.route_costs[crew, ix] = cost
    
    # append the fires fought
    route_data.fires_fought[crew, ix, :, :] = fires_fought
    
    # append the out-of-region assignments
    route_data.out_of_reg[crew, ix, :] = oor
    
    return 1

end

function get_supp_plan_stats(var_p, var_d, beta, tolerance=0.0001)
    
    # get the cost based on the perimeter progression
    cost = beta * (sum(value.(var_p)) - value(var_p[1])/2 - value(var_p[NUM_TIME_PERIODS+1]/2))
    
    # get the number of crews present each time period from line constructed
    crew_vector = value.(var_d)
    int_crew_vector = convert.(Int64, round.(crew_vector))
    @assert maximum(abs.(crew_vector - int_crew_vector)) < tolerance "Not an integer plan"
    
    return cost, int_crew_vector

end

function initialize_supp_plan_data(max_supp_plans)
    
    return SuppressionPlanData(zeros(NUM_FIRES), 
                               Matrix{Float64}(undef, NUM_FIRES, max_supp_plans),
                               zeros(Int8, (NUM_FIRES, max_supp_plans, NUM_TIME_PERIODS))
                              )
end

function update_available_supp_plans(fire, cost, allotment, plan_data)
    
    # add 1 to number of plans for this fire, store the index
    plan_data.plans_per_fire[fire] += 1
    ix = plan_data.plans_per_fire[fire]
    
    # append the route cost
    plan_data.plan_costs[fire, ix] = cost
    
    # append the fires fought
    plan_data.crews_present[fire, ix, :] = allotment
    
    return 1

end

function full_formulation(integer_routes, region_data, constraint_data, rotation_order, 
                          costs, progs, perims, beta, gamma, verbose, time_limit)
    
    # get number of arcs
    n_arcs = length(costs)
    
    # shorten some global variable names
    C = NUM_CREWS
    G = NUM_FIRES
    T = NUM_TIME_PERIODS
    regs = region_data.crew_regions
    
    # intialize model
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    
    if ~verbose
        set_optimizer_attribute(m, "OutputFlag", 0)
    end
    
    if time_limit != false
        set_optimizer_attribute(m, "TimeLimit", time_limit)
    end
        

    # fire suppression plan section
    @variable(m, p[g=1:G, t=1:T+1] >= 0)
    @variable(m, l[g=1:G, t=1:T])
    @variable(m, d[g=1:G, t=1:T] >= 0)
    @constraint(m, perim_growth[g=1:G, t=1:T], p[g, t+1] >= progs[g, t] * 
                                                           (p[g, t] - l[g, t] / 2) - l[g, t] / 2)
    @constraint(m, perim_start[g=1:G], p[g, 1] == perims[g])

    # routing plan section
    if integer_routes
        @variable(m, z[1:n_arcs] >= 0, Int)
    else
        @variable(m, z[1:n_arcs] >= 0)
    end
    
    if integer_routes
        @variable(m, q[1:C, 0:T] >= 0, Int)
    else
        @variable(m, q[1:C, 0:T] >= 0)
    end
    
    # build out_of_region constraints
    if gamma > 0
        @constraint(m, out_of_region[c=1:C, t=0:T],

            # out of region penalty is at least
            q[c, t] >=

                # this crew's cumulative rotations
                sum(z[i] for t_0=0:t, i in constraint_data.out_of_region[c, t_0+1]) 

            - 

                # average cumulative rotations among all crews in same region
                sum(z[i] for c_0 in keys(rotation_order[regs[c]]), t_0=0:t, 
                    i in constraint_data.out_of_region[c_0, t_0+1]) /
                length(keys(rotation_order[regs[c]]))

            -

                # normalizing factor for specific crew rotation order
                (1 - rotation_order[regs[c]][c] / length(keys(rotation_order[regs[c]])))

            -
                # big-M for if crew goes not leave region at this time
                T * (1 - sum(z[i] for i in constraint_data.out_of_region[c, t+1]))

        )
    else
        @constraint(m, out_of_region[c=1:C, t=0:T], 0 == 0)
    end


    @constraint(m, fire_flow[c=1:C, g=1:G, t=1:T, rest=1:2],

            sum(z[constraint_data.f_out[c, g, t, rest]]) ==
            sum(z[constraint_data.f_in[c, g, t, rest]])
    
    )
    
    @constraint(m, base_flow[c=1:C, t=1:T, rest=1:2],

            sum(z[constraint_data.b_out[c, t, rest]]) ==
            sum(z[constraint_data.b_in[c, t, rest]])
    
    )


    @constraint(m, linking[g=1:G, t=1:T],

        LINE_PER_CREW * sum(z[constraint_data.supp_fire[g, t]]) >= l[g, t]
    )
    

    # build start constraint
    @constraint(m, start[c=1:C], 

        sum(z[constraint_data.start[c]]) == 1
    )
    
    
    

    @objective(m, Min, 
        beta * (sum(p) - sum(p[1:G, 1])/2 - sum(p[1:G, T+1])/2) + 
        sum(z .* costs) + sum(q) * gamma
    )
    
    return m, p, l, z, q, out_of_region, linking
    
end

function load_data(path)
    
    # get distance from fire f to fire g 
    fire_dists =  readdlm(path * "/fire_distances.csv", ',')

    # get distance from base c to fire g (NUM_CREWS-by-NUM_FIRES)
    base_fire_dists =  readdlm(path * "/base_fire_distances.csv", ',')

    # initialize travel times (number of periods) from fire f to fire g
    tau = convert(Array{Int}, ones(size(fire_dists)))

    # initialize number of periods to travel from base c to fire g (NUM_CREWS-by-NUM_FIRES)
    tau_base_to_fire = convert(Array{Int}, ones((size(base_fire_dists))))

    # read intial crew statuses (location, period by which they must rest)
    # (-1 in current_fire means crew is currently at base)
    # (rested_periods is the amount of time crew has been at base, relevant for completing rest)
    crew_starts = CSV.read(path * "/sample_crew_starts.csv", DataFrame)
    rest_by = crew_starts[!, "rest_by"]
    current_fire = crew_starts[!, "current_fire"]
    rested_periods = crew_starts[!, "rested_periods"]


    return (GlobalData(fire_dists, base_fire_dists, tau, tau_base_to_fire), 
            CrewStatus(rest_by, current_fire, rested_periods))
end

function update_master_problem(mp, route_data, supp_plan_data, new_routes, new_plans)
    
    # for each crew where we found a new route
    for crew in new_routes
        
        # push it to the model
        ix = route_data.routes_per_crew[crew]
        mp["route"][crew, ix] = @variable(mp["m"], base_name="route[$crew,$ix]", lower_bound=0)
        
        # update coefficient in objective
        set_objective_coefficient(mp["m"], mp["route"][crew, ix], route_data.route_costs[crew, ix])
        
        # update coefficient in constraints
        set_normalized_coefficient(mp["sigma"][crew], mp["route"][crew, ix], 1)
        set_normalized_coefficient.(mp["rho"], mp["route"][crew, ix], route_data.fires_fought[crew, ix, :, :])
        
        ## TODO out_of_region ##
        
    end
        
    # for each fire where we found a new plan
    for fire in new_plans
        
        # push it to the model
        ix = supp_plan_data.plans_per_fire[fire]
        mp["plan"][fire, ix] = @variable(mp["m"], base_name="plan[$fire,$ix]", lower_bound=0)
        
        # update coefficient in objective
        set_objective_coefficient(mp["m"], mp["plan"][fire, ix], supp_plan_data.plan_costs[fire, ix])
        
        # update coefficient in constraints
        set_normalized_coefficient(mp["pi"][fire], mp["plan"][fire, ix], 1)
        set_normalized_coefficient.(mp["rho"][fire, :], mp["plan"][fire, ix], 
                                    -supp_plan_data.crews_present[fire, ix, :])
        
    end 
    

    return mp
    
end

function master_problem(col_gen_config, route_data, supp_plan_data, region_data, rotation_order, gamma, price_branch)
    
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(m, "OutputFlag", 0)
    set_optimizer_attribute(m, "OptimalityTol", 1e-9)
    set_optimizer_attribute(m, "FeasibilityTol", 1e-9)
    
    regs = region_data.crew_regions
    
    # decision variables
    if price_branch
        @variable(m, route[c=1:NUM_CREWS, r=1:route_data.routes_per_crew[c]] >= 0, Int)
        @variable(m, plan[g=1:NUM_FIRES, p=1:supp_plan_data.plans_per_fire[g]] >= 0, Int)
        @variable(m, q[c=1:NUM_CREWS, t=0:NUM_TIME_PERIODS] >= 0, Int)
    else
        @variable(m, route[c=1:NUM_CREWS, r=1:route_data.routes_per_crew[c]] >= 0)
        @variable(m, plan[g=1:NUM_FIRES, p=1:supp_plan_data.plans_per_fire[g]] >= 0)
        @variable(m, q[c=1:NUM_CREWS, t=0:NUM_TIME_PERIODS] >= 0)
    end
    
    # if we are doing a dual_stabilization
    if col_gen_config["master_problem"]["dual_stabilization"] != false
        
        # should not do this in price-and-branch framework
        @assert price_branch == false
        
        # linking constraint perturbations
        @variable(m, delta_plus[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS] >= 0)
        @variable(m, delta_minus[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS] >= 0)
        
    end
    
    # constraints that you must choose a plan per crew and per fire
    @constraint(m, route_per_crew[c=1:NUM_CREWS], 
                sum(route[c, r] for r=1:route_data.routes_per_crew[c]) == 1)
    @constraint(m, plan_per_fire[g=1:NUM_FIRES], 
                sum(plan[g, p] for p=1:supp_plan_data.plans_per_fire[g]) >= 1)
    
    # linking constraint
    
    if col_gen_config["master_problem"]["dual_stabilization"] == false
        @constraint(m, linking[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS],

                        # crews at fire
                        sum(route[c, r] * route_data.fires_fought[c, r, g, t] 
                            for c=1:NUM_CREWS, r=1:route_data.routes_per_crew[c]) 

                        >=

                        # crews suppressing
                        sum(plan[g, p] * supp_plan_data.crews_present[g, p, t] 
                            for p=1:supp_plan_data.plans_per_fire[g]) 

                    )
    elseif col_gen_config["master_problem"]["dual_stabilization"] == "global"
       
        # get expected dual value ratios
        ratios = col_gen_config["master_problem"]["dual_warm_start"]
        ratios = ratios / sum(ratios)
        
        # get secondary dual objective epsilon
        secondary_eps = col_gen_config["master_problem"]["dual_secondary_eps"]
        
        @constraint(m, linking[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS],

                # crews at fire
                sum(route[c, r] * route_data.fires_fought[c, r, g, t] 
                    for c=1:NUM_CREWS, r=1:route_data.routes_per_crew[c]) 
            
                +
            
                # perturbation
                delta_plus[g, t] - delta_minus[g, t] - 
                sum(ratios .* delta_plus) + sum(ratios .* delta_minus)
                
                >=

                # crews suppressing
                sum(plan[g, p] * supp_plan_data.crews_present[g, p, t] 
                    for p=1:supp_plan_data.plans_per_fire[g]) 

            )
        @constraint(m, perturb[g=1:NUM_FIRES, t=1:NUM_TIME_PERIODS], 
                    delta_plus[g, t] + delta_minus[g, t] <= secondary_eps)
    else
        error("Dual stabilization type not implemented")
    end
    

    # out_of_region constraint
    @constraint(m, out_of_region[c=1:NUM_CREWS, t=0:NUM_TIME_PERIODS],
    
        # out of region penalty is at least
        q[c, t] >=
        
            # this crew's cumulative rotations
            sum(route[c, r] * route_data.out_of_reg[c, r, t_0 + 1] 
            for r=1:route_data.routes_per_crew[c], t_0=0:t)
        
        - 
        
            # average cumulative rotations among all crews in same region
            sum(route[c_0, r] * route_data.out_of_reg[c_0, r, t_0 + 1] 
                for c_0 in keys(rotation_order[regs[c]]), r=1:route_data.routes_per_crew[c_0],
                t_0=0:t) /
            length(keys(rotation_order[regs[c]]))
        
        -
        
            # normalizing factor for specific crew rotation order
            (1 - rotation_order[regs[c]][c] / length(keys(rotation_order[regs[c]])))
        
        -
            # big-M for if crew goes not leave region at this time
            NUM_TIME_PERIODS * (1 - sum(route[c, r] * route_data.out_of_reg[c, r, t+1]
                                        for r=1:route_data.routes_per_crew[c])
                               )
        
    )
    
    @objective(m, Min, 
        
                  # route costs
                  sum(route[c, r] * route_data.route_costs[c, r] 
                        for c=1:NUM_CREWS, r=1:route_data.routes_per_crew[c])
        
                  +
                     
                  # suppression plan costs
                  sum(plan[g, p] * supp_plan_data.plan_costs[g, p] 
                     for g=1:NUM_FIRES, p=1:supp_plan_data.plans_per_fire[g]) 
        
                  +
        
                  # rotational queueing violations cost
                  sum(q) * gamma
               )
    
    return Dict("m" => m, "q" => q, "sigma" => route_per_crew, "pi" => plan_per_fire, 
                "rho" => linking, "eta" => out_of_region, "route" => route, "plan" => plan)
end 

function init_route_subproblem(crew_ixs, crew, constraint_data, integer_routes=false)
    
    # shorten some global variable names
    C = NUM_CREWS
    G = NUM_FIRES
    T = NUM_TIME_PERIODS
    
    # intialize model
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(m, "OutputFlag", 0)

    # routing plan section
    if integer_routes
        @variable(m, z[crew_ixs] >= 0, Int)
    else
        @variable(m, z[crew_ixs] >= 0)
    end


    @constraint(m, fire_flow[g=1:G, t=1:T, rest=1:2],

            sum(z[constraint_data.f_out[crew, g, t, rest]]) ==
            sum(z[constraint_data.f_in[crew, g, t, rest]])
    
    )
    
    @constraint(m, base_flow[t=1:T, rest=1:2],

            sum(z[constraint_data.b_out[crew, t, rest]]) ==
            sum(z[constraint_data.b_in[crew, t, rest]])
    
    )

    # build start constraint
    @constraint(m, start, 

        sum(z[constraint_data.start[crew]]) == 1
    )
    
    return Dict("m" => m, "z" => z, "ff" => fire_flow)
    
end

function init_suppression_plan_subproblem(progs, perims, fire, beta)
    
    T = NUM_TIME_PERIODS
    
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(m, "OutputFlag", 0)

    # fire suppression plan section
    @variable(m, p[t=1:T+1] >= 0)
    @variable(m, l[t=1:T] >= 0)
    @variable(m, NUM_CREWS >= d[t=1:T] >= 0, Int)
    @constraint(m, suppression_per_crew[t=1:T], l[t] <= d[t] * LINE_PER_CREW)
    @constraint(m, perim_growth[t=1:T], p[t+1] >= progs[fire, t] * (p[t] - l[t] / 2) - l[t] / 2)
    @constraint(m, perim_start, p[1] == perims[fire])
    
#    
    
    return Dict("m" => m, "p" => p, "d" => d, "beta" => beta)
end

function full_network_flow(region_data, constraint_data, rotation_order, 
                               r_arc_costs, f_arc_costs, crews_needed, fire_linking_arcs, 
                               fire_start_arcs, f_out_arcs, f_in_arcs, gamma, integer=false, verbose=false)
    
    # get number of arcs
    n_arcs_route = length(r_arc_costs)
    n_arcs_fire = length(f_arc_costs)
    
    # shorten some global variable names
    C = NUM_CREWS
    G = NUM_FIRES
    T = NUM_TIME_PERIODS
    regs = region_data.crew_regions
    
    # intialize model
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    
    if ~verbose
        set_optimizer_attribute(m, "OutputFlag", 0)
    end

    # fire suppression plan section
    
    if integer
        @variable(m, y[i=1:n_arcs_fire] >= 0, Int)
    else
        @variable(m, y[i=1:n_arcs_fire] >= 0)
    end
        
    @variable(m, d[fire=1:NUM_FIRES, t=1:NUM_TIME_PERIODS] >= 0)

    @constraint(m, fire_linking[fire=1:NUM_FIRES, t=1:NUM_TIME_PERIODS],

        d[fire, t] >= sum(crews_needed[i] * y[i] for i in fire_linking_arcs[fire][t])

                )
    
    @constraint(m, fire_start[fire=1:NUM_FIRES], y[fire_start_arcs[fire]] == 1)
    
    n_states = size(f_out_arcs)[1]
    @constraint(m, state_flow[s=1:n_states, t=1:NUM_TIME_PERIODS], sum(y[f_out_arcs[s, t]]) == sum(y[f_in_arcs[s, t]]))

    if integer
        @variable(m, z[1:n_arcs_route] >= 0, Int)
    else
        @variable(m, z[1:n_arcs_route] >= 0)
    end
    
    @variable(m, q[1:C, 0:T] >= 0)
    
    # build out_of_region constraints
    @constraint(m, out_of_region[c=1:C, t=0:T],
    
        # out of region penalty is at least
        q[c, t] >=
        
            # this crew's cumulative rotations
            sum(z[i] for t_0=0:t, i in constraint_data.out_of_region[c, t_0+1]) 
        
        - 
        
            # average cumulative rotations among all crews in same region
            sum(z[i] for c_0 in keys(rotation_order[regs[c]]), t_0=0:t, 
                i in constraint_data.out_of_region[c_0, t_0+1]) /
            length(keys(rotation_order[regs[c]]))
        
        -
        
            # normalizing factor for specific crew rotation order
            (1 - rotation_order[regs[c]][c] / length(keys(rotation_order[regs[c]])))
        
        -
            # big-M for if crew goes not leave region at this time
            T * (1 - sum(z[i] for i in constraint_data.out_of_region[c, t+1]))
        
    )


    @constraint(m, fire_flow[c=1:C, g=1:G, t=1:T, rest=1:2],

            sum(z[constraint_data.f_out[c, g, t, rest]]) ==
            sum(z[constraint_data.f_in[c, g, t, rest]])
    
    )
    
    @constraint(m, base_flow[c=1:C, t=1:T, rest=1:2],

            sum(z[constraint_data.b_out[c, t, rest]]) ==
            sum(z[constraint_data.b_in[c, t, rest]])
    
    )


    @constraint(m, route_linking[g=1:G, t=1:T],

        sum(z[constraint_data.supp_fire[g, t]]) >= d[g, t] 
    )
    
    # build start constraint
    @constraint(m, start[c=1:C], 

        sum(z[constraint_data.start[c]]) == 1
    )
    
    
    
    

    @objective(m, Min, 
        sum(y .* f_arc_costs) + 
        sum(z .* r_arc_costs) + sum(q) * gamma
    )
    
    return Dict("m" => m, "d" => d, "z" => z, "y" => y, "q" => q, "oor" => out_of_region, 
                "r_linking" => route_linking, "f_linking" => fire_linking)
    
end

function warm_start_suppression_plans(plans_per_fire, fire_model_configs, region_data, constraint_data, rotation_order, 
                               r_arc_costs, gamma, integer, verbose, time_limit)
    
    all_fire_configs = [Dict{String,Any}("solver_type" => "network_flow_gurobi") for fire in 1:NUM_FIRES]
    
    fire_sps = []
    for fire in 1:NUM_FIRES
        config = all_fire_configs[fire]
        config["model_data"] = fire_model_configs[fire]
        config["warm_start"] = "lp_relax"
        d = init_suppression_plan_subproblem(config)
        push!(fire_sps, d)
    end
    
    round_type = "ceiling"
    arcs_per_sp = [size(fire_sps[fire][round_type]["arcs"])[1] for fire in 1:NUM_FIRES]
    
    for fire in 1:NUM_FIRES
        fire_arcs = fire_sps[fire][round_type]["arcs"]
        fire_arcs[:, 1] .= fire
    end
    concat_fire_arcs = reduce(vcat, [fire_sps[fire][round_type]["arcs"] for fire in 1:NUM_FIRES])
    n_fire_arc = size(concat_fire_arcs)[1]


    fire_linking = []
    for fire in 1:NUM_FIRES
        fire_arcs = fire_sps[fire][round_type]["arcs"]
        push!(fire_linking, [[j for j in 1:size(fire_arcs)[1] if fire_arcs[j, 3] == i] for i in 1:NUM_TIME_PERIODS])
    end


    offset = 0
    for fire in 1:NUM_FIRES
        for i in 1:length(fire_linking[fire])
            fire_linking[fire][i] = fire_linking[fire][i] .+ offset
        end
        offset = offset + arcs_per_sp[fire]
    end
    
    all_out_arcs = reduce(vcat, [fire_sps[fire][round_type]["out_arcs"] for fire in 1:NUM_FIRES])
    all_in_arcs = reduce(vcat, [fire_sps[fire][round_type]["in_arcs"] for fire in 1:NUM_FIRES])

    state_offsets = [size(fire_sps[fire][round_type]["out_arcs"])[1] for fire in 1:NUM_FIRES]
    state_offsets  = cumsum(state_offsets)

    for all_arcs in [all_out_arcs, all_in_arcs]
        for i in 1:size(all_arcs)[1]
            for fire in 1:NUM_FIRES
                if i > state_offsets[fire]
                    for j in 1:NUM_TIME_PERIODS
                        all_arcs[i, j] = all_arcs[i, j] .+ arcs_per_sp[fire]
                    end
                end
            end
        end
    end
    
    crew_requirements = concat_fire_arcs[:, 6]
    fire_arc_costs = reduce(vcat, [fire_sps[fire][round_type]["costs"] for fire in 1:NUM_FIRES])
    start_arcs = reduce(vcat, all_in_arcs[:, 1])
    
    network_flow_model = full_network_flow(region_data, constraint_data, rotation_order, 
                               r_arc_costs, fire_arc_costs, crew_requirements, fire_linking, 
                               start_arcs, all_out_arcs, all_in_arcs, gamma, integer, verbose)
    
    if integer
        set_optimizer_attribute(network_flow_model["m"], "TimeLimit", time_limit)
        set_optimizer_attribute(network_flow_model["m"], "MIPFocus", 1)
    end
        
    optimize!(network_flow_model["m"])
    
    d = network_flow_model["d"]
    network_flow_model["arcs"] = concat_fire_arcs
    network_flow_model["arc_costs"] = fire_arc_costs
    
    return [suppression_plan_perturbations(ceil.(value.(d[fire, :]) .- 0.00001), plans_per_fire) for fire in 1:NUM_FIRES], network_flow_model
end

function arcs_from_state_graph(graph)
    
    visitable = [(i,j) for i in 1:size(graph)[1], j in 1:size(graph)[2] if isassigned(graph, i, j)]
    edges = []

    for (i, j) in visitable
        push!(edges, copy(reduce(hcat, [[i, j, j+1, a[1], a[2]] for a in graph[i, j]])'))
    end

    return vcat([size(graph)[1], 0, 1, size(graph)[1], 0]', reduce(vcat, edges))
end

function get_fire_cost(crew_allocations, params)
    
    if params["model_type"] == "simple_linear"
        
        perims = zeros(NUM_TIME_PERIODS + 1)
        perims[1] = params["start_perim"]
        for i in 1:NUM_TIME_PERIODS
            perims[i+1] = update_fire_stats(perims[i], i, crew_allocations[i], params)
        end
        
        return params["beta"] * (sum(perims) - 0.5 * perims[1] - 0.5 * perims[end])
        
    else
        error("Model type not implemented")
    end
end

function update_fire_stats(curr_stats, curr_time, crew_allocation, params)
    
    if params["model_type"] == "simple_linear"
        
        line_per_crew = params["line_per_crew"]
        prog = params["progressions"][curr_time]
        line = line_per_crew * crew_allocation
        next_stats = (curr_stats - line/2) * prog - line/2
        next_stats = max(0, next_stats)
        
    else
        error("Not implemented")
    end
    
    return next_stats
end

function inverse_update_fire_stats(stats_from, stats_to, time_from, time_to, params)
    """ Returns number of crews needed for given fire transition """
    
    if params["model_type"] == "simple_linear"

        line_per_crew = params["line_per_crew"]
        prog = params["progressions"][time_from]
        
        crews = 2 / line_per_crew * (prog * stats_from - stats_to) / (1 + prog)
        
    else
        error("Not implemented")
    end
    
    return crews
    
end

function create_discrete_fire_states(params)
    
    if params["model_type"] == "simple_linear"
        
        # get the no-suppression progression of this fire
        progs = params["progressions"]
        start_perim = params["start_perim"]
        no_supp = [start_perim]
        for i in 1:NUM_TIME_PERIODS
            push!(no_supp, no_supp[i] * progs[i])
        end
        
        # generalize this later
        aggressive_precision = AGG_PREC
        num_aggressive_states = convert(Int, round(start_perim * 2 / aggressive_precision))
        num_passive_states = PASSIVE_STATES

        aggressive_states = LinRange(0, num_aggressive_states * aggressive_precision, num_aggressive_states)
        passive_states = exp.(LinRange(log(num_aggressive_states * aggressive_precision + 1), 
                                       maximum(log.(no_supp .+ 1)), num_passive_states + 1)
                             )
        passive_states = passive_states[2:num_passive_states+1] .- 1
        all_states = vcat(aggressive_states, passive_states)
        all_states = vcat(all_states, 9999999)
        
        push!(all_states, start_perim)
        
    else
        error("Not implemented")
    end
    
    return all_states
end

function generate_state_transition_crew_reqs(curr_stats, curr_time, sorted_states, params, round_types)
    
    if params["model_type"] == "simple_linear"
        
        # get all possibly feasible states
        min_state_val = update_fire_stats(curr_stats, curr_time, NUM_CREWS + 0.5, params)
        max_state_val = update_fire_stats(curr_stats, curr_time, 0, params)
        
        # get min and max index of the possibly feasible states
        min_state_ix = searchsorted(sorted_states, min_state_val).start
        max_state_ix = searchsorted(sorted_states, max_state_val).start
        
        # inititalize output 
        output = Dict()
        for round_type in round_types
            output[round_type] = []
        end
        
        # inititalize minimum crews needed so far, for trimming (explained below)
        min_used = Dict()
        for round_type in round_types
            min_used[round_type] = NUM_CREWS + 1
        end
        
        # for each feasible state
        for state_ix in min_state_ix:max_state_ix
            crews_needed = inverse_update_fire_stats(curr_stats, sorted_states[state_ix], curr_time, 
                                                     curr_time + 1, params)
            
            # for each round type
            for round_type in round_types
                
                # round the number of crews
                if round_type == "ceiling"
                    crews = max(0, convert(Int, ceil(crews_needed - 0.0001)))
                elseif round_type == "floor"
                    crews = max(0, convert(Int, floor(crews_needed + 0.0001)))
                elseif round_type == "nearest"
                    crews = max(0, convert(Int, round(crews_needed)))
                else
                    error("Round type not implemented")
                end
                
                # since the states are sorted in increasing level of cost and risk
                # we can trim arcs that are dominated
                
                # if this is a feasible number of crews and we do not have a dominating arc
                if (crews <= NUM_CREWS) & (crews < min_used[round_type])
                    
                    # update the minimum crews used for this round type
                    min_used[round_type] = crews
                    
                    # push the crew requirement for this transition to the corresponding list
                    push!(output[round_type], (state_ix, crews))
                end
            end
        end
        
        return output

    
    else
        error("Model type not implemented")
    end
    
end

function generate_graphs(states, params, round_types)
    
    # separate out non start states
    non_start_states = states[1:length(states)-1]
    
    # initializations
    n_states = length(non_start_states)
    crews_needed = Dict()
    
    # for each round type
    for round_type in round_types
        
        # init graph as array of vectors (each index is a state*time vertex, each vector holds edges out)
        crews_needed[round_type] = Array{Vector{}}(undef, n_states + 1, NUM_TIME_PERIODS + 1)
        
        # initialize time, state indices to check next
        curr_time = 1
        next_to_check = [n_states + 1]
        
        # while we have not yet reached the end of the horizon
        while curr_time < NUM_TIME_PERIODS + 1
            
            # copy the indices to check and reset next_to_check
            to_check = copy(next_to_check)
            next_to_check = []

            # for each state index feasible to reach at this time period
            for check in to_check

                # if it is not the last (no return) state
                if check != n_states
                    
                    # find the crew requirements for feasible edges
                    edges = generate_state_transition_crew_reqs(states[check], curr_time, non_start_states, 
                                                                params, [round_type])[round_type]
                
                # otherwise stay at this state
                else
                    edges = [(n_states, 0)]
                end

                # update the crews needed for each edge coming out of this state
                crews_needed[round_type][check, curr_time] = edges

                # append the neighbors to next_to_check
                next_to_check = vcat(next_to_check, [edges[i][1] for i in 1:length(edges) 
                                                     if ~(edges[i][1] in next_to_check)])
            end
            
            # add one to the time
            curr_time += 1
        end
    end
    
    return crews_needed
end

function create_local_edge_price_func(rho, edge_caps, break_cap_penalty)

    function local_edge_price(curr_time, num_crews)
    
        out = rho[curr_time] * num_crews
        if num_crews > edge_caps[curr_time]
            out += break_cap_penalty * (num_crews - edge_caps[curr_time])
        end
        
        return out
    end
    
    return local_edge_price
end 

function solve_fire_dp(graph, avail_nodes, state_enter_costs, local_edge_price_func)
    
    n_states = size(state_enter_costs)[1]
    
    curr_time = NUM_TIME_PERIODS
    costs = zeros(n_states, NUM_TIME_PERIODS + 1) .+ 1.0e12
    costs[:, NUM_TIME_PERIODS + 1] = state_enter_costs[:, NUM_TIME_PERIODS + 1]
    bests = Dict()
    
    # backward induction
    while curr_time > 0
    
        state_time_costs = state_enter_costs[:, curr_time]
        nodes_to_check = avail_nodes[curr_time]
        
        for node in nodes_to_check
            current_cost = 1e30
            current_best = -1
            current_allot = -1
            edges = graph[node, curr_time]
            for edge in edges
                edge_cost = costs[edge[1], curr_time + 1] + local_edge_price_func(curr_time, edge[2])
                if edge_cost < current_cost
                    current_best = edge[1]
                    current_allot = edge[2]
                    current_cost = edge_cost
                end
            end
            costs[node, curr_time] = current_cost + state_time_costs[node]
            bests[(node, curr_time)] = (current_best, current_allot);
        end

        curr_time -= 1
    end
    
    # recreate path
    curr_state = n_states
    curr_time = 1
    curr_cost = state_enter_costs[curr_state, curr_time]
    allotments = convert.(Int, zeros(NUM_TIME_PERIODS))
    all_states = convert.(Int, zeros(NUM_TIME_PERIODS + 1))
    all_states[1] = curr_state

    while curr_time < NUM_TIME_PERIODS + 1

        next_state = bests[(curr_state, curr_time)]
        curr_state = next_state[1]
        allotments[curr_time] = next_state[2]
        curr_time += 1
        all_states[curr_time] = curr_state
        curr_cost += state_enter_costs[curr_state, curr_time]
    end
    
    return allotments, all_states, curr_cost, costs[n_states, 1]
end

function get_state_entrance_cost(state, enter_time, params)
    
    if params["model_type"] == "simple_linear"
        
        if (enter_time == 1) | (enter_time == NUM_TIME_PERIODS + 1)
            return params["beta"] * state / 2
        elseif (enter_time > 1) & (enter_time < NUM_TIME_PERIODS + 1)
            return params["beta"] * state
        else
            error("Invalid time")
        end
    else
        error("Not implemented")
    end
end

function init_suppression_plan_subproblem(config)
             
    if config["solver_type"] == "dp"

        states = create_discrete_fire_states(config["model_data"])
        graphs = generate_graphs(states, config["model_data"], ["ceiling", "floor"])
        nodes_avail = Dict()
        nodes_avail["floor"] = Dict(i => [j for j in 1:size(graphs["floor"])[1] 
                                          if isassigned(graphs["floor"], j, i)] 
                                    for i in 1:size(graphs["floor"])[2])
        nodes_avail["ceiling"] = Dict(i => [j for j in 1:size(graphs["ceiling"])[1] 
                                            if isassigned(graphs["ceiling"], j, i)] 
                                      for i in 1:size(graphs["ceiling"])[2])
        state_costs = zeros(length(states), NUM_TIME_PERIODS + 1)
        for i = 1:length(states)
            for j = 1:NUM_TIME_PERIODS + 1
                state_costs[i, j] = get_state_entrance_cost(states[i], j, config["model_data"])
            end
        end
        
        return Dict("graphs" => graphs, "avail_nodes" => nodes_avail, "state_costs" => state_costs)
        
    elseif config["solver_type"] == "network_flow_gurobi"
        
        # get graph formulation
        d = copy(config)
        d["solver_type"] = "dp"
        sp_dp = init_suppression_plan_subproblem(d)
        
        output = Dict{String, Any}()
        
        for strategy in keys(sp_dp["graphs"])
            
            graph = sp_dp["graphs"][strategy]

            arc_array = arcs_from_state_graph(graph)
            n_arcs = length(arc_array[:, 1])

            # add crew to front (actually, this is deprecated, adding just -1's to keep column order)
            arc_array = hcat(convert.(Int, zeros(length(arc_array[:, 1]))) .- 1, arc_array)

            state_costs = sp_dp["state_costs"]
            arc_costs = [state_costs[arc_array[i, 5], arc_array[i, 4]] for i in 1:n_arcs]


            n_states = size(state_costs)[1]  

            out_arcs = Array{Vector{Int64}}(undef, n_states, NUM_TIME_PERIODS)
            in_arcs = Array{Vector{Int64}}(undef, n_states, NUM_TIME_PERIODS)

            for s in 1:n_states
                state_arcs = [i for i in 1:n_arcs if arc_array[i, 2] == s]
                for t in 1:NUM_TIME_PERIODS
                    out_arcs[s, t] = [i for i in state_arcs if arc_array[i, 3] == t]
                end
            end

            for s in 1:n_states
                state_arcs = [i for i in 1:n_arcs if arc_array[i, 5] == s]
                for t in 1:NUM_TIME_PERIODS
                    in_arcs[s, t] = [i for i in state_arcs if arc_array[i, 4] == t]
                end
            end

            if config["warm_start"] == "lp_relax"
                output[strategy] = Dict("arcs" => arc_array, "costs" => arc_costs, 
                                        "out_arcs" => out_arcs, "in_arcs" => in_arcs)
            else
                # intialize model
                m = Model(() -> Gurobi.Optimizer(GRB_ENV))
                set_optimizer_attribute(m, "OutputFlag", 0)

                @variable(m, z[1:n_arcs] >= 0, Int)
                @constraint(m, flow[s=1:n_states, t=1:NUM_TIME_PERIODS], sum(z[out_arcs[s, t]]) == sum(z[in_arcs[s, t]]))
                @constraint(m, start, z[1] == 1)


                output[strategy] = Dict("arcs" => arc_array, "costs" => arc_costs, "m" => m, "z" => z, 
                                        "out_arcs" => out_arcs, "in_arcs" => in_arcs)
            end
        end
        
        return output
        
    elseif config["solver_type"] == "gurobi"
        
        if config["model_data"]["model_type"] == "simple_linear"

            progs = config["model_data"]["progressions"]
            perim = config["model_data"]["start_perim"]
            beta = config["model_data"]["beta"]

            T = NUM_TIME_PERIODS

            m = Model(() -> Gurobi.Optimizer(GRB_ENV))
            set_optimizer_attribute(m, "OutputFlag", 0)

            # fire suppression plan section
            @variable(m, p[t=1:T+1] >= 0)
            @variable(m, l[t=1:T] >= 0)
            @variable(m, NUM_CREWS >= d[t=1:T] >= 0, Int)
            @constraint(m, suppression_per_crew[t=1:T], l[t] <= d[t] * LINE_PER_CREW)
            @constraint(m, perim_growth[t=1:T], p[t+1] >= progs[t] * (p[t] - l[t] / 2) - l[t] / 2)
            @constraint(m, perim_start, p[1] == perim)


            return Dict("m" => m, "p" => p, "d" => d, "beta" => beta)
        else
            error("Model type not implemented for Gurobi")
        end
    else
        error("Solver type not implemented")
    end
end
      

function run_fire_subproblem(sp, config, rho)
    
    if config["solver_type"] == "gurobi"
        
        if config["model_data"]["model_type"] == "simple_linear"
            
            m = sp["m"]
            p = sp["p"]
            d = sp["d"]
            beta = sp["beta"]
            
            if config["warm_start"] == "dummy"
                @objective(m, Min, sum(d))
            elseif config["warm_start"] == false
                @objective(m, Min, beta * (sum(p) - p[1]/2 - p[NUM_TIME_PERIODS + 1]/2) + 
                                   sum(d .* rho) + 0.0001 * sum(d))
            else
                error("warm start type not implemented")
            end
            
            optimize!(m)
            
            # get the required information from the model decision variables
            cost, crew_vector = get_supp_plan_stats(p, d, beta)
            
            return cost, objective_value(m), crew_vector
        else
            error("model type not implemented")
        end
        
    elseif config["solver_type"] == "dp"
        strategy = config["solver_strategy"]
        graph = sp["graphs"][strategy]
        avail_nodes = sp["avail_nodes"][strategy]
        state_enter_costs = sp["state_costs"]
        
        if config["warm_start"] == "dummy"
            duals = zeros(length(rho)) .+ 1e30
        elseif config["warm_start"] == false
            duals = rho
        else
            error("warm start type not implemented")
        end
        
        # get capacities
        capacities = config["capacities"]
        break_capacity_penalty = config["break_capacity_penalty"]
    
        # make edge price function
        edge_prices = create_local_edge_price_func(duals, capacities, break_capacity_penalty)
        
        # solve dp
        allotments, states, cost, rel_cost = solve_fire_dp(graph, avail_nodes, state_enter_costs, edge_prices)
        
        return cost, rel_cost, allotments
        
    elseif config["solver_type"] == "network_flow_gurobi"
        strategy = config["solver_strategy"]
        model_data = sp[strategy]
        m = model_data["m"]
        z = model_data["z"]
        arc_costs = model_data["costs"]
        arc_array = model_data["arcs"]
        
        if config["warm_start"] == "dummy"
            duals = zeros(length(rho)) .+ 1e30
        elseif config["warm_start"] == false
            duals = rho
        else
            error("warm start type not implemented")
        end
        
        # no cost to start edge
        duals = vcat([0], duals)
        
        adj_costs = arc_costs .+ [arc_array[i, 6] * duals[arc_array[i, 3] + 1] for i in 1:length(arc_costs)]
        @objective(m, Min, sum(adj_costs .* z))
        optimize!(m)
        
        rel_cost = objective_value(m)
        arcs_used = [i for i in 1:length(adj_costs) if value(z[i]) > 0.9]
        cost = sum(arc_costs[arcs_used])
        allotments = arc_array[arcs_used[2:length(arcs_used)], 6]

        return cost, rel_cost, allotments
    else
        error("solver type not implemented")
    end
end

function initialize_column_generation(arcs, costs, constraint_data, fire_model_configs, solver_configs, max_plans)
    
    # initialize subproblems
    route_sps = []
    for crew in 1:NUM_CREWS
        ixs = [i for i in 1:length(arcs[:, 1]) if arcs[i, 1] == crew]
        d = init_route_subproblem(ixs, crew, constraint_data)
        d["arc_ixs"] = ixs
        push!(route_sps, d)
    end
    
    plan_sps = []
    for fire in 1:NUM_FIRES
        config = copy(solver_configs[fire])
        config["model_data"] = fire_model_configs[fire]
        config["warm_start"] = false
        d = init_suppression_plan_subproblem(config)
        push!(plan_sps, d)
    end
    
    # initialize routes and suppression plans to populate
    routes = initialize_route_data(max_plans)
    suppression_plans = initialize_supp_plan_data(max_plans)
    
    ## generate dummy plans (no suppression) to ensure feasibility at first step ##
    
    # for each crew
    for crew in 1:NUM_CREWS
        
        # get the crew's subproblem instance
        crew_sp = route_sps[crew]
        m = crew_sp["m"]
        z = crew_sp["z"]
        crew_ixs = crew_sp["arc_ixs"]

        # set objective in light of dual variables
        @objective(m, Min, sum(z[ix] * (costs[ix]) for ix in crew_ixs))

        # optimize
        optimize!(m)

        # update crew routes
        crew_arcs = [i for i in crew_ixs if (value(z[i]) > 0.5)]
        update_available_routes(crew, crew_arcs, arcs, costs, routes)
    
    end
    
    # for each fire
    for fire in 1:NUM_FIRES
        
        sp_config = copy(solver_configs[fire])
        sp_config["model_data"] = fire_model_configs[fire]
        sp_config["solver_strategy"] = "ceiling"
        sp_config["warm_start"] = "dummy"
        
        if sp_config["solver_type"] == "dp"
            sp_config["capacities"] = zeros(NUM_TIME_PERIODS)
            sp_config["break_capacity_penalty"] = 0
        end
        
        cost, rel_cost, allotment = run_fire_subproblem(plan_sps[fire], sp_config, zeros(NUM_TIME_PERIODS))

        # update suppression plans
        update_available_supp_plans(fire, cost, allotment, suppression_plans)

    end
    
    return ColumnGeneration(route_sps, plan_sps, routes, suppression_plans)
    
end

function run_crew_subproblem(sps, crew, costs, local_costs)
    
    # get the crew's subproblem instance
    crew_sp = sps[crew]
    m = crew_sp["m"]
    z = crew_sp["z"]
    crew_ixs = crew_sp["arc_ixs"]
    
    # set objective in light of dual variables
    @objective(m, Min, sum(z[ix] * (local_costs[ix] + costs[ix]) for ix in crew_ixs))
        
    # optimize
    optimize!(m)
    
    return objective_value(m), z
end

function run_CG_step(cg, arcs, costs, global_data, region_data, fire_model_configs, solver_configs, cg_config,
                     rot_order, gamma, recover_fire_sp_cost, mp)
    
    last_num_routes = copy(cg.routes.routes_per_crew)
    last_num_plans = copy(cg.suppression_plans.plans_per_fire) 
    
    t = @elapsed optimize!(mp["m"])
    obj = objective_value(mp["m"])
    rho = dual.(mp["rho"])
    allotments = get_fire_allotments(mp, cg)
    # println("solve")
    # println(t)

    # grab the dual variables
    sigma = dual.(mp["sigma"])
    rho = dual.(mp["rho"])
    eta = dual.(mp["eta"])
    pie = dual.(mp["pi"]) # lol can't overwrite "pi" in Julia
    

    
    true_rho = copy(rho)
    
    if "ws_dual_weight" in keys(cg_config)
        if cg_config["ws_dual_weight"] > 0
            error("Wrong type of dual warm start (useless stabilization)")
            # lambda = cg_config["ws_dual_weight"]
            # ws_dual_values = cg_config["ws_dual"]
            # rho = (lambda * (ws_dual_values)) .+ ((1 - lambda) * rho)
        end
    end
    
    

    # using the dual variables, get the local adjustments to the arc costs in the route subproblems
    d = Dict("out_of_region_dual" => eta, "region_data"=> region_data, "rotation_order" => rot_order, "linking_dual" => rho)
    local_costs = get_arc_costs(global_data, arcs, d)

    ## run subproblems ##

    # for each fire
    for fire in 1:NUM_FIRES

        # run the subproblem
        sp_config = copy(solver_configs[fire])
        sp_config["model_data"] = fire_model_configs[fire]
        sp_config["solver_strategy"] = "ceiling"
        sp_config["warm_start"] = false
        
        found_plan = false
        for adjust_price in sp_config["int_aware_adjustment_pattern"]
            
            if ~found_plan
                if "int_aware_capacities" in keys(sp_config)
                    allotments = sp_config["int_aware_capacities"]
                else
                    # grab the allotments
                    allotments = get_fire_allotments(mp, cg) 
                end
                
                sp_config["capacities"] = (allotments .+ adjust_price[1])[fire, :]
                sp_config["capacities"] = max.(sp_config["capacities"], 0)
                sp_config["break_capacity_penalty"] = adjust_price[2]

        
                cost, modified_rel_cost, allotment = run_fire_subproblem(cg.plan_sps[fire], sp_config, rho[fire, :])
                 
                rel_cost = cost + sum(allotment .* true_rho[fire, :])

                # if there is an improving plan
                if rel_cost < pie[fire] - 0.0001

                    found_plan = true

                    # adjust cost to true value if needed
                    if recover_fire_sp_cost
                        cost = get_fire_cost(allotment, fire_model_configs[fire])
                    end

                    # add the plan
                    update_available_supp_plans(fire, cost, allotment, cg.suppression_plans)
                end
            end
        end
    end

    # for each crew
    for crew in 1:NUM_CREWS

        # run the crew subproblem
        obj, assignments = run_crew_subproblem(cg.route_sps, crew, costs, local_costs)
                 

        # if there is an improving route
        if obj < sigma[crew] - 0.0001

            # add it
            crew_arcs = [i for i in cg.route_sps[crew]["arc_ixs"] if (value(assignments[i]) > 0.5)]
            update_available_routes(crew, crew_arcs, arcs, costs, cg.routes)

        end

    end 
    
    # formulate and solve the master problem
    t = 0
    t += @elapsed plans = findall(last_num_plans .< cg.suppression_plans.plans_per_fire)
    t += @elapsed routes = findall(last_num_routes .< cg.routes.routes_per_crew)
    t += @elapsed mp = update_master_problem(mp, cg.routes, cg.suppression_plans, routes, plans)
    # println("formulate")
    # println(t)
    return mp, obj, rho, allotments
end

function get_fire_allotments(solved_mp, cg_data_object)
    
    mp_allotment = zeros(size(cg_data_object.suppression_plans.crews_present[:, 1, :]))
    
    for plan in eachindex(solved_mp["plan"])
        new_allot = cg_data_object.suppression_plans.crews_present[plan[1], plan[2], :] * value(solved_mp["plan"][plan])
        mp_allotment[plan[1], :] += new_allot
    end
    
    return mp_allotment
end


function suppression_plan_perturbations(start_plan, count)
    
    # get the indices we may perturb, chosen to be anything at most one index away
    # from a time period when suppression was >0 for the start_plan
    ixs = [i for i in 1:NUM_TIME_PERIODS if start_plan[i] > 0]
    ixs_1 = [i+1 for i in 1:NUM_TIME_PERIODS-1 if start_plan[i] > 0]
    ixs_2 = [i-1 for i in 2:NUM_TIME_PERIODS if start_plan[i] > 0]
    ixs_to_perturb = sort(unique(vcat(ixs, ixs_1, ixs_2)))
    
    start_plan_copy = copy(start_plan)
    
    # hack if no fire suppression
    if length(ixs) == 0
        start_plan_copy[1] += 1
        ixs_to_perturb = 1:NUM_TIME_PERIODS
        ixs = [1]
    end
    
    found = []
    
    # perturbing the total number of crews by 0, -1, or 1
    for perturb in [0, 1, -1]
        
        # make sure we don't explode
        curr_length = length(found)
        
        # push a new possible plan
        new_plan = copy(start_plan_copy)
        new_plan[ixs[1]] = new_plan[ixs[1]] + perturb
        push!(found, copy(new_plan))
        
        # find all ways to perturb the valid indices while saying within
        # prescribed crew bounds and keeping total crews unchanged
        for current_arr in found
            if length(found) < curr_length + 500
                nonzero = [ix for ix in ixs_to_perturb if current_arr[ix] > 0]
                nonfull = [ix for ix in ixs_to_perturb if current_arr[ix] < NUM_CREWS]
                for ix in nonzero
                    for ix2 in nonfull
                        if ix2 != ix
                            next_arr = copy(current_arr)
                            next_arr[ix] -= 1
                            next_arr[ix2] += 1
                            if !(next_arr in found)
                                push!(found, copy(next_arr))
                            end
                        end
                    end
                end
            end
        end
    end
    
    # get the closest plans to the original, using L1 norm
    ixs_to_keep = sortperm([sum(abs.(i - start_plan)) for i in found])[1:min(count, length(found))]
    
    return found[ixs_to_keep]

end 