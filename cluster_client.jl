include("cluster_test_helper.jl")
include("DCG.jl")

using Profile
using JSON
const GRB_ENV = Gurobi.Env()
const GRAPH_TYPE = "ceiling"

function run_experiment(date, json_name, write_output)
    
    # get path to input parameters
    json_path = "data/experiment_inputs/" * date * "/" * json_name * ".json"
    
    # read in parameters
    external_params = JSON.parse(open(json_path))

    # create full dict of parameters
    dcg_params = default_params()
    for key in keys(external_params["dcg"])
        dcg_params[key] = external_params["dcg"][key]
    end
    for key in keys(external_params)
        if key in keys(dcg_params)
            dcg_params[key] = external_params[key]
        end
    end


    # update global variables
    global NUM_FIRES = dcg_params["num_fires"]
    global NUM_CREWS = dcg_params["num_crews"]
    global LINE_PER_CREW = dcg_params["line_per_crew"]

    # process all the relevant location, fire data
    t = @elapsed preprocessed_data = preprocess(dcg_params["in_path"], dcg_params["agg_prec"], dcg_params["passive_states"])
    println(t)

    # run DCG at root node, saving dual warm start for future iterations
    d, cg_data = single_DCG_node(dcg_params, deepcopy(preprocessed_data))
    @assert dcg_params["dual_warm_start"] != "calculate"

    plan_usages = []
    for plan in eachindex(d["mp"]["plan"])
        if value(d["mp"]["plan"][plan]) > 0
            push!(plan_usages, 
            (plan[1], cg_data.suppression_plans.crews_present[plan[1], plan[2], :], value(d["mp"]["plan"][plan])))
        end
    end

    route_usages = []
    for route in eachindex(d["mp"]["route"])
        if value(d["mp"]["route"][route]) > 0
            push!(route_usages, 
            (route[1], cg_data.routes.fires_fought[route[1], route[2], :, :], value(d["mp"]["route"][route])))
        end
    end

    d["plan_usages"] = plan_usages
    d["route_usages"] = route_usages
    d["routes"] = deepcopy(cg_data.routes.routes_per_crew)
    d["plans"] = deepcopy(cg_data.suppression_plans.plans_per_fire)
    
    # update int-aware capacities from weighted average primal solution
    dcg_params["int_aware_capacities"] = d["allotments"]["master_problem_reconstructed"]
    
    # generate int-aware plans
    # set perturbations in int-aware-plan generating phase
    capacity_perturbations = [-2, -1, 0, 1, 2] .* (NUM_CREWS / 20)
    iterations, timings, cg_data = generate_new_plans(dcg_params, preprocessed_data, cg_data, capacity_perturbations, -1, 0)
    
    # restore integrality
    form_time, sol_time, pb = restore_integrality(cg_data, 300)
    pb_allot = convert.(Int, round.(100 .* get_fire_allotments(pb, cg_data)) ./ 100);
    
    # write output to JSON
    if write_output
        
        out_dir = "data/experiment_outputs/" * date * "/"
        if (~isdir(out_dir))
            mkpath(out_dir)
        end
        
        outputs = Dict{String, Any}()
        delete!(d, "mp")
        outputs["initial_DCG"] = d
        outputs["generate_additional_plans"] = Dict{String, Any}()
        outputs["generate_additional_plans"]["iterations"] = iterations
        outputs["generate_additional_plans"]["timings"] = timings
        outputs["cg_data"] = Dict{String, Any}()
        outputs["cg_data"]["routes_per_crew"] = cg_data.routes.routes_per_crew
        outputs["cg_data"]["plans_per_fire"] = cg_data.suppression_plans.plans_per_fire
        outputs["restore_integrality"] = Dict{String, Any}()
        outputs["restore_integrality"]["formulation_time"] = form_time
        outputs["restore_integrality"]["solve_time"] = sol_time
        outputs["restore_integrality"]["pb_objective"] = objective_value(pb["m"])
        outputs["restore_integrality"]["pb_objective_bound"] = objective_bound(pb["m"])
        outputs["restore_integrality"]["allotment"] = pb_allot
        
        
        
        open(out_dir * json_name * ".json", "w") do f
            JSON.print(f, outputs, 4)
        end
  
    end
    
    return cg_data
end
    
# run_experiment("20221123", "1", false)
date = ARGS[1]
number = ARGS[2]
run_experiment(string(date), "precompile", false)
run_experiment(string(date), string(number), true)

