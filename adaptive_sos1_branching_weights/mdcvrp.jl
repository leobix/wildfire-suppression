""" 
This code is adapted from Alex Schmid's 15.095 recitation 
"""

using LinearAlgebra, JuMP, Gurobi, LightGraphs, Printf

const GRB_ENV = Gurobi.Env()

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
function new_route_mdcvrp(distances::Matrix, dual_variables::Vector, Q::Int)

    n = size(distances, 1) - 1 # number of customers
    @assert length(dual_variables) == n
    push!(dual_variables, 0) # depot

    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attributes(model, "TimeLimit" => 10, "OutputFlag" => 0)
    @variable(model, x[i=1:(n+1), j=1:(n+1)], Bin)
    @constraint(model, flow_conservation[i=1:(n+1)],
        sum(x[j, i] for j = 1:(n+1)) == sum(x[i, j] for j = 1:(n+1)))
    @constraint(model, max_flow[i=1:n],
        sum(x[i, j] for j = 1:(n+1)) <= 1)
    @constraint(model, no_self_flow[i=1:(n+1)], x[i, i] == 0)
    @constraint(model, into_depot, sum(x[j, n+1] for j = 1:(n+1)) == 1)

    # TODO apply branching rules just like above

    @objective(model, Min, sum(x[i, j] * (distances[i, j] - dual_variables[i]) for i = 1:(n+1), j = 1:(n+1)))

    # CAPACITY CONSTRAINT
    @constraint(model, sum(x[i, j] for i = 1:(n+1), j = 1:(n+1)) <= Q)

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

    route = [n+1]
    while (length(route) == 1) || (route[end] != route[1])
        push!(route, argmax(value.(x[route[end], :])))
    end
    reduced_cost = sum(distances[route[i], route[i+1]] - dual_variables[route[i]] for i = 1:(length(route)-1))

    return route, sum(distances[route[i], route[i+1]] for i = 1:(length(route)-1)), reduced_cost
end

"Solve OVRP relaxation using column generation"
function rmp_mdcrvp(distances::Matrix, Q::Int; T::Int=10, verbose::Bool=false)

    n = size(distances, 1) - 1 # number of customers
    routes, costs = [], []

    t = 0
    upper_bounds = []
    lower_bounds = []

    while true

        t += 1

        model = Model(() -> Gurobi.Optimizer(GRB_ENV))
        set_optimizer_attributes(model, "InfUnbdInfo" => 1, "TimeLimit" => 10, "OutputFlag" => 0)
        @variable(model, y[r=eachindex(routes)] >= 0)
        @constraint(model, visit_customer[i=1:n],
            sum(y[r] for r = eachindex(routes) if i in routes[r]) >= 1)
        @objective(model, Min, sum(costs[r] * y[r] for r = eachindex(routes)))

        optimize!(model)

        push!(upper_bounds, objective_value(model))
        p = [dual(visit_customer[i]) for i = 1:n]
        if termination_status(model) != MOI.OPTIMAL # farkas vector
            p = p .* 1000
        end
        println("here")
        route, cost, reduced_cost = new_route_mdcvrp(distances, p, Q)

        verbose && @printf("Found column with reduced cost: %.2g\n", reduced_cost)
        verbose && print(route)
        push!(lower_bounds, objective_value(model) + n * reduced_cost)

        if (t > 1) && (reduced_cost > -1e-6 || t > T)
            return [routes[r] for r = eachindex(routes) if value(y[r]) >= 0.1], objective_value(model),
            upper_bounds, lower_bounds
        end

        push!(routes, route)
        push!(costs, cost)

    end
end


locations = [[1, 1], [1, 10], [3, 5]]
N = length(locations)
demands = [4, 5, 6]
vehicle_bases = [[1, 3], [2, 7]]
M = length(vehicle_bases)
capacities = [9, 6]

dist = [norm(locations[i, :] .- locations[j, :]) for i = 1:N, j = 1:N]
# base_dist = [norm(locations[i, :] .- vehicle_bases[j, :]) for i = 1:N, j = 1:M]
# @printf base_dist

@time routes, total_cost, upper, lower = rmp_mdcrvp(dist, 4, T=10, verbose = true)
