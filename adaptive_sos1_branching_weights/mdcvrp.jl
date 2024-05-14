""" 
This code is adapted from Alex Schmid's 15.095 recitation 
"""

using LinearAlgebra, JuMP, Gurobi, LightGraphs, Printf

const GRB_ENV = Gurobi.Env()

# THIS COULD BE MORE GENERAL TO ALLOW DIFFERENT CONSTRAINTS TO DETERMINE THE WEIGHTS 
# IN DIFFERENT VRP VARIANTS
struct SOSBranchingRule
    vehicle::Int
    customer::Int
    visits::Bool
end

function satisfies(route::Vector{Int}, vehicle::Int, branch::SOSBranchingRule)
    if vehicle != branch.vehicle
        return true
    end

    if branch.visits
        return branch.customer ∈ route
    end
    return branch.customer ∉ route
end

function satisfies(route::Vector{Int}, vehicle::Int, branches::Vector{SOSBranchingRule})
    for branch ∈ eachindex(branches)
        if ~satisfies(route, vehicle, branch)
            return false
        end
    end
    return true
end

""""
Create initial set of routes and associated costs
"""
function initial_routes(distances::Matrix)
    n = size(distances, 1) - 1 # number of customers
    return [[i, n + 1] for i in 1:n], distances[1:n, n+1]
end

"""
    Given the induced graph as an adjacency list (i.e., next[i] is the next node to visit after node i),
        compute all subtours.
    Return them as a list of lists of nodes in the same component
"""
function find_subtours(next::Vector{Int})
    n = length(next)
    g = DiGraph(n)
    for i = 1:n
        if next[i] != 0
            add_edge!(g, i, next[i])
        end
    end
    components = strongly_connected_components(g)
    return sort(components, by=length)
end

"Find next column to add, returns route, cost, and whether or not the reduced cost is nonnegative"
function new_route_mdcvrp(distances::Matrix,
    dual_variables::Vector,
    demands::Vector,
    branches::Vector{SOSBranchingRule},
    Q::Int)

    n = size(distances, 1) - 1 # number of customers
    @assert length(dual_variables) == n
    @assert length(demands) == n
    duals = copy(dual_variables)
    push!(duals, 0) # depot

    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "InfUnbdInfo", 1)
    @variable(model, x[i=1:(n+1), j=1:(n+1)], Bin)
    @constraint(model, flow_conservation[i=1:(n+1)],
        sum(x[j, i] for j = 1:(n+1)) == sum(x[i, j] for j = 1:(n+1)))
    @constraint(model, max_flow[i=1:n],
        sum(x[i, j] for j = 1:(n+1)) <= 1)
    @constraint(model, no_self_flow[i=1:n], x[i, i] == 0) # allow depot self-flow so there is a naive solution
    @constraint(model, into_depot, sum(x[j, n+1] for j = 1:(n+1)) == 1)

    @constraint(model, branching_rules[b=1:length(branches)],
        sum(x[j, branches[b].customer] for j = 1:(n+1)) == Int(branches[b].visits))

    @objective(model, Min, sum(x[i, j] * (distances[i, j] - duals[i]) for i = 1:(n+1), j = 1:(n+1)))

    # CAPACITY CONSTRAINT
    @constraint(model, sum(demands[j] * x[i, j] for i = 1:(n+1), j = 1:n) <= Q)

    "Define the callback function"
    function eliminate_subtours(cb_data)
        status = callback_node_status(cb_data, model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            # get value of current solution
            next = zeros(Int, n + 1)
            for i = 1:(n+1), j = 1:(n+1)
                if callback_value(cb_data, x[i, j]) > 0.5
                    next[i] = j
                end
            end
            subtours = find_subtours(next)
            # solve dual subproblems
            for subtour in subtours
                if length(subtour) == 1 || n + 1 in subtour
                    continue
                else
                    cut = @build_constraint(sum(x[subtour[i], subtour[i+1]] for i = 1:(length(subtour)-1)) +
                                            x[subtour[length(subtour)], subtour[1]] <= length(subtour) - 1)
                    MOI.submit(model, MOI.LazyConstraint(cb_data), cut)
                end
            end
        end
    end

    # set callback function and attach to model
    MOI.set(model, MOI.LazyConstraintCallback(), eliminate_subtours)
    optimize!(model)

    next = zeros(Int, n + 1)
    for i = 1:(n+1), j = 1:(n+1)
        if value(x[i, j]) > 0.5
            next[i] = j
        end
    end

    route = [n + 1]
    while (length(route) == 1) || (route[end] != route[1])
        push!(route, argmax(value.(x[route[end], :])))
    end
    reduced_cost = sum(distances[route[i], route[i+1]] - duals[route[i]] for i = 1:(length(route)-1))

    return route, sum(distances[route[i], route[i+1]] for i = 1:(length(route)-1)), reduced_cost
end

"Solve OVRP relaxation using column generation"
function rmp_mdcrvp(distances::Matrix, depot_distances::Matrix, branches::Vector{SOSBranchingRule}, capacities::Vector, demands::Vector; T::Int=10, verbose::Bool=false)

    n = size(distances, 1)  # number of customers
    d = size(depot_distances, 1) # number of depots
    @assert length(capacities) == d
    @assert length(demands) == n
    routes = Containers.@container(x[i=1:d, j=1:1000; 0 != 0], Int[])
    costs = Containers.@container(x[i=1:d, j=1:1000; 0 != 0], 0.0)
    vehicle = 1

    t = 0
    upper_bounds = []
    lower_bounds = []
    valid_routes = [r for r ∈ eachindex(routes) if satisfies(routes[r], r[1], branches)]
    while true

        t += 1

        # TODO refactor model formulation outside of loop
        model = Model(() -> Gurobi.Optimizer(GRB_ENV))
        set_optimizer_attribute(model, "OutputFlag", 0)
        set_optimizer_attribute(model, "InfUnbdInfo", 1)
        @variable(model, y[i=1:d, j=1:1000; (i, j) ∈ valid_routes] >= 0)
        @variable(model, 0 <= visit_customer[i=1:n] <= 1)
        @constraint(model, sos1[i=1:d], sum(y[r] for r ∈ valid_routes if r[1] == i) <= 1)
        @constraint(model, linking[i=1:n],
            sum(y[r] for r = valid_routes if i ∈ routes[r] && r ∈ valid_routes) >= visit_customer[i])
        @objective(model, Min, sum(costs[r] * y[r] for r ∈ eachindex(routes)) - 1000 * sum(visit_customer) + n * 1000)

        optimize!(model)

        push!(upper_bounds, objective_value(model))
        p = dual.(linking)
        σ = dual.(sos1)
        if termination_status(model) != MOI.OPTIMAL # farkas vector
            p = p ./ sum(p) .* 1000
            println("infeasibility certificate fa vector")
            println(p)
        end

        rc_sum = 0
        found_col = false
        for vehicle ∈ 1:d
            dists = vcat(distances, depot_distances[vehicle, :]')
            dists = hcat(dists, vcat(depot_distances[vehicle, :], 0))
            route, cost, reduced_cost = new_route_mdcvrp(dists, p, demands, [b for b in branches if b.vehicle == vehicle], capacities[vehicle])
            reduced_cost -= σ[vehicle]
            rc_sum += reduced_cost

            if reduced_cost <= -1e-6
                found_col = true
                verbose && @printf("Found column with reduced cost: %.2g\n", reduced_cost)
                ix = length(routes[vehicle, :]) + 1
                routes[vehicle, ix] = route
                costs[vehicle, ix] = cost
                push!(valid_routes, (vehicle, ix))
            end
        end
        push!(lower_bounds, objective_value(model) + rc_sum)

        if (t > 1) && (~found_col || t > T)
            visit_allots = zeros(Float64, n, d)
            for i ∈ 1:n
                for j ∈ 1:d
                    for r ∈ eachindex(y)
                        if r[1] == j
                            visit_allots[i, j] += (value(y[r]) * normalized_coefficient(linking[i], y[r]))
                        end
                    end
                end
            end
            println(visit_allots)
            println(valid_routes)
            println(value.(y))
            println(routes)
            println(objective_value(model))
            println(p)
            println(σ)
            return [routes[r] for r = eachindex(y) if value(y[r]) >= 0.1], objective_value(model),
            upper_bounds, lower_bounds
        end
    end
end


locations = [[1, 1], [1, 20], [3, 5], [1, 2], [1, 9], [3, 6]]
N = length(locations)
demands = [4, 3, 2, 1, 2, 3]
vehicle_bases = [[10, 50], [2, 7]]
M = length(vehicle_bases)
capacities = [8, 8]

dist = [norm(locations[i, :] .- locations[j, :]) for i = 1:N, j = 1:N]
base_dist = [norm(locations[i, :] .- vehicle_bases[j, :]) for j = 1:M, i = 1:N]

b = SOSBranchingRule(2, 2, false)
routes, total_cost, upper, lower = rmp_mdcrvp(dist, base_dist, [b], capacities, demands, T=20, verbose=true)
@time routes, total_cost, upper, lower = rmp_mdcrvp(dist, base_dist, [b], capacities, demands, T=20, verbose=true)
