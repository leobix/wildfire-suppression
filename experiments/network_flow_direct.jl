include("../Subproblems.jl")

using JuMP, Gurobi, JSON, ArgParse
const GRB_ENV = Gurobi.Env()

function get_command_line_args()
    arg_parse_settings = ArgParseSettings()
    @add_arg_table arg_parse_settings begin
        "--debug"
        help = "run in debug mode, exposing all logging that uses @debug macro"
        action = :store_true
        "--directory_output", "-d"
        help = "directory to write outputs, must exist"
        arg_type = String
        default = "data/experiment_outputs/network_flow_direct/"
    end
    return parse_args(arg_parse_settings)
end

function full_network_flow(
	crew_models::Vector{TimeSpaceNetwork},
	fire_models::Vector{TimeSpaceNetwork};
	integer = false,
	verbose = false,
	time_limit = 180.0)

	ub = Inf
	lb = 0
	num_crews = length(crew_models)
	num_fires = length(fire_models)
	_, num_times, _ = size(crew_models[1].state_in_arcs)


	# intialize model
	m = Model(() -> Gurobi.Optimizer(GRB_ENV))
	if ~verbose
		set_optimizer_attribute(m, "OutputFlag", 0)
	end
	set_optimizer_attribute(m, "TimeLimit", time_limit)

	fire_vars = []
	for fire ∈ 1:num_fires
		fire_model = fire_models[fire]
		y = @variable(
			m,
			[1:size(fire_model.long_arcs)[1]],
			lower_bound = 0,
			upper_bound = 1,
		)
		# println(size(fire_model.long_arcs)[1])
		# for i ∈ 1:size(fire_model.long_arcs)[1]
		#     set_objective_coefficient(m, y[i], fire_model.arc_costs[i])
		# end
		if integer
			set_binary.(y)
		end
		push!(fire_vars, y)
	end

	crew_vars = []
	for crew ∈ 1:num_crews
		crew_model = crew_models[crew]
		ixs = findall(crew_model.long_arcs[:, CM.CREW_NUMBER] .== crew)
		z = @variable(
			m,
			[ixs],
			lower_bound = 0,
			upper_bound = 1,
		)
		# println(length(crew_model.arc_costs[ixs]))
		# set_objective_coefficient.(m, z, crew_model.arc_costs[ixs])
		if integer
			set_binary.(z)
		end
		push!(crew_vars, z)
	end

	@objective(
		m,
		Min,
		sum(
			fire_models[fire].arc_costs[ix] * fire_vars[fire][ix] for
			fire ∈ 1:num_fires,
			ix ∈ 1:size(fire_models[fire].long_arcs)[1]
		) + sum(
			crew_models[crew].arc_costs[ix] * crew_vars[crew][ix] for
			crew ∈ 1:num_crews,
			ix ∈ findall(crew_models[crew].long_arcs[:, CM.CREW_NUMBER] .== crew)
		)
	)

	# fire network flow
	@constraint(
		m,
		fire_flow[
			fire = 1:num_fires,
			t = 1:num_times,
			s = 1:size(fire_models[fire].state_out_arcs)[1],
		],
		sum(
			fire_vars[fire][ix] for ix ∈ fire_models[fire].state_out_arcs[s, t]
		) ==
		sum(fire_vars[fire][ix] for ix ∈ fire_models[fire].state_in_arcs[s, t])
	)

	# fire start
	@constraint(
		m,
		fire_start[fire = 1:num_fires],
		fire_vars[fire][1] == 1
	)

	# crew network flow
	locs, times, rests = size(crew_models[1].state_in_arcs)
	@constraint(
		m,
		crew_flow[crew = 1:num_crews, l = 1:locs, t = 1:times, r = 1:rests],
		sum(crew_vars[crew][ix] for ix ∈ crew_models[crew].state_in_arcs[l, t, r]) ==
		sum(crew_vars[crew][ix] for ix ∈ crew_models[crew].state_out_arcs[l, t, r])
	)
	# crew start
	@constraint(
		m,
		crew_start[crew = 1:num_crews],
		sum(
			crew_vars[crew][ix] for
			ix ∈ findall(crew_models[crew].long_arcs[:, CM.TIME_FROM] .== 0) if
			crew_models[crew].long_arcs[ix, CM.CREW_NUMBER] == crew)
		==
		1
	)

	# linking
	@constraint(
		m,
		linking[fire = 1:num_fires, t = 1:num_times],
		sum(
			crew_vars[crew][ix] for crew ∈ 1:num_crews,
			ix ∈ vcat(
				crew_models[crew].state_in_arcs[fire, t, 1],
				crew_models[crew].state_in_arcs[fire, t, 2],
			)
		) >=
		sum(
			fire_models[fire].long_arcs[ix, FM.CREWS_PRESENT] .* fire_vars[fire][ix]
			for ix ∈ eachindex(fire_vars[fire]) if
			fire_models[fire].long_arcs[ix, FM.TIME_FROM] == t
		)
	)

	optimize!(m)

	if has_values(m)
		ub = objective_value(m)
		lb = objective_bound(m)
	end

	return lb, ub
	# supplies = zeros(num_fires, num_time_periods)
	# for fire ∈ 1:num_fires
	# 	for t ∈ 1:num_time_periods
	# 		supp = 0
	# 		for crew ∈ 1:num_crews
	# 			for ix ∈ crew_models[crew].state_in_arcs[fire, t, 1]
	# 				supp += value(crew_vars[crew][ix])
	# 			end
	# 			for ix ∈ crew_models[crew].state_in_arcs[fire, t, 2]
	# 				supp += value(crew_vars[crew][ix])
	# 			end
	# 		end
	# 		supplies[fire, t] = supp
	# 	end
	# end

end
args = get_command_line_args()
linear_outputs = Dict()
integer_outputs = Dict()

num_crews = 10
num_fires = 3
num_time_periods = 14

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

# for model ∈ crew_models
# 	println(size(model.long_arcs))
# end

# for model ∈ fire_models
# 	println(size(model.long_arcs))
# end

full_network_flow(crew_models, fire_models, verbose = false)


sizes = [(3, 10, 14), (6, 20, 14), (9, 30, 14), (12, 40, 14), (15, 50, 14)]
sizes = [(15, 50, 14)]

for (num_fires, num_crews, num_time_periods) ∈ sizes

	local_crew_models = build_crew_models(
		"data/raw/big_fire",
		num_fires,
		num_crews,
		num_time_periods,
	)

	local_fire_models = build_fire_models(
		"data/raw/big_fire",
		num_fires,
		num_crews,
		num_time_periods,
	)



	t = @elapsed l, u = full_network_flow(
		local_crew_models,
		local_fire_models,
		verbose = false,
		integer = false,
		time_limit = 60.0,
	)
	linear_outputs[num_crews] = Dict("ub" => u, "lb" => l, "time" => t)
	t = @elapsed l, u = full_network_flow(
		local_crew_models,
		local_fire_models,
		verbose = false,
		integer = true,
		time_limit = 120.0,
	)
	integer_outputs[num_crews] = Dict("ub" => u, "lb" => l, "time" => t)
end

open(args["directory_output"] * "output.json", "w") do f
	JSON.print(f, Dict("linear" => linear_outputs, "integer" => integer_outputs), 4)
end


