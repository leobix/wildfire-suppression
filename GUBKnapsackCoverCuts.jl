using IterTools

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

struct GUBCutsWithSubproblemInfo

    cuts::JuMP.Containers.SparseAxisArray
    cut_objects::Dict{Any, Any}
    cuts_per_time::Vector{Int64}
    fire_cut_sp_lookup::Dict{Int64, Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}}
    crew_cut_sp_lookup::Dict{Int64, Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}}

end

function get_gub_fire_relevant_suppression_data(f, t, fire_plans, rmp, used_plans)
	sample = fire_plans.crews_present[f, used_plans[f], :]
	values = [value(rmp.plans[f, i]) for i in used_plans[f]]
	coeff_order = sortperm(sample[:, t])
	allot = 0
	weight = 0
	pairs = Tuple[]
	for coeff in coeff_order
		if sample[coeff, t] > allot
			push!(pairs, (sample[coeff, t], weight))
			allot = sample[coeff, t]
			# weight = 0
		end
		weight += values[coeff]
	end
	return pairs
end

function get_gub_crew_relevant_routing_data(c, t, crew_routes, used_routes, rmp)
	sample =
		.~dropdims(
			maximum(crew_routes.fires_fought[c, used_routes[c], :, :], dims = 2),
			dims = 2,
		)
	values = [value(rmp.routes[c, i]) for i in used_routes[c]]
	coeff_order = sortperm(sample[:, t])
	allot = 0
	weight = 0
	pairs = Tuple[]
	for coeff in coeff_order
		if sample[coeff, t] > allot
			push!(pairs, (sample[coeff, t], weight))
			allot = sample[coeff, t]
			# weight = 0
		end
		weight += values[coeff]
	end
	return pairs
end

function enumerate_minimal_cuts(crew_allots, fire_allots)

	all_cuts = []

	# find the crews that are inactive at time t in all used routes
	inactive_crews = Int64[]
	for crew in eachindex(crew_allots)
		if (length(crew_allots[crew]) > 0)
			if (crew_allots[crew][1][1] == true) & (crew_allots[crew][1][2] == 0)
				push!(inactive_crews, crew)
			end
		end
	end

	# for now, assume that these crews were inactive for a "good" reason and force those to be part of any cover
	available_for_fires = length(crew_allots) - length(inactive_crews)

	### get all minimal GUB covers under this assumption ###

	# enumerate indices of all possible GUB covers 
	for i in eachindex(fire_allots)
		fire_allots[i] = vcat([(0, 0.0)], fire_allots[i])
	end
	allotment_option_ixs = []
	for i in eachindex(fire_allots)
		push!(allotment_option_ixs, [])
		for j in eachindex(fire_allots[i])
			push!(allotment_option_ixs[i], j)
		end
	end

	loop_count = 0
	for i in product(allotment_option_ixs...)
		loop_count += 1
		if loop_count >= 1000
			@info "Broke cover enumeration early" loop_count
			break
		end
		allot = 0
		cost = 0
		min_allot = 10000
		for j in eachindex(i)
			cur_allot = fire_allots[j][i[j]][1]
			allot += cur_allot
			cost += fire_allots[j][i[j]][2]
			if i[j] > 1
				min_allot = min(min_allot, cur_allot)
			end
		end
		if (cost < 1 - 0.02) & (allot > available_for_fires) &
		   (allot - min_allot <= available_for_fires)
			cut_allotments = [fire_allots[j][i[j]][1] for j in eachindex(i)]
			@debug "cut found" i,
			available_for_fires,
			allot,
			min_allot,
			cost,
			cut_allotments
			push!(
				all_cuts,
				(available_for_fires, cost, cut_allotments, inactive_crews),
			)
		end
	end

	return all_cuts
end

function get_fire_and_crew_incumbent_weighted_average(rmp, crew_routes, fire_plans)

	# get problem dimensions
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	fire_allotment = zeros(num_fires, num_time_periods)
	for ix in eachindex(rmp.plans)
		if value(rmp.plans[ix]) > 0
			fire_allotment[ix[1], :] +=
				fire_plans.crews_present[ix..., :] * value(rmp.plans[ix])
		end
	end

	crew_allotment = zeros(num_crews, num_fires, num_time_periods)
	for ix in eachindex(rmp.routes)
		if value(rmp.routes[ix]) > 0
			crew_allotment[ix[1], :, :] +=
				crew_routes.fires_fought[ix..., :, :] * value(rmp.routes[ix])
		end
	end

	fire_allotment, crew_allotment
end

function adjust_cut_fire_allotment(
	current_allotment,
	incumbent_weighted_average,
	budget,
)

	adjusted_allotment = copy(current_allotment)
	while sum(adjusted_allotment) > budget + 1
		max_diff = -1000.0
		max_diff_ix = -1
		for ix in eachindex(adjusted_allotment)
			if (adjusted_allotment[ix] > 0) &
			   (adjusted_allotment[ix] - incumbent_weighted_average[ix] > max_diff)
				max_diff_ix = ix
				max_diff = adjusted_allotment[ix] - incumbent_weighted_average[ix]
			end
		end
		adjusted_allotment[max_diff_ix] -= 1
	end
	return adjusted_allotment
end

function generate_knapsack_cuts(crew_routes, fire_plans, rmp)

	# get problem dimensions
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	all_fire_allots = Array{Vector{Tuple}, 2}(undef, num_fires, num_time_periods)
	all_crew_allots = Array{Vector{Tuple}, 2}(undef, num_crews, num_time_periods)

	# get the plans that are used in the current RMP solution
	used_plans = []
	for g in 1:num_fires
		plans = [
			i[1] for
			i in eachindex(rmp.plans[g, :]) if value(rmp.plans[g, i...]) > 1e-4
		]
		push!(used_plans, plans)
	end

	# get the routes that are used in the current RMP solution
	used_routes = []
	for c in 1:num_crews
		routes = [
			i[1] for
			i in eachindex(rmp.routes[c, :]) if value(rmp.routes[c, i...]) > 1e-4
		]
		push!(used_routes, routes)
	end

	# get the current fire suppression relaxed solution
	fire_alloc, _ =
		get_fire_and_crew_incumbent_weighted_average(rmp, crew_routes, fire_plans)

	# intitialize output 
	knapsack_gub_cuts = []

	# for each time period
	@debug fire_alloc
	for t in 1:num_time_periods

		# process these routes and plans into GUB-relevant allotments (needed for dominated cuts due to GUB constraints)
		for f ∈ 1:num_fires
			all_fire_allots[f, t] = get_gub_fire_relevant_suppression_data(
				f,
				t,
				fire_plans,
				rmp,
				used_plans,
			)
		end
		for c ∈ 1:num_crews
			all_crew_allots[c, t] = get_gub_crew_relevant_routing_data(
				c,
				t,
				crew_routes,
				used_routes,
				rmp,
			)
		end

		minimal_cuts =
			enumerate_minimal_cuts(all_crew_allots[:, t], all_fire_allots[:, t])

		for cut in minimal_cuts
			orig_fire_allots = cut[3]
			unused_crews = cut[4]
			adjusted_cut =
				adjust_cut_fire_allotment(orig_fire_allots, fire_alloc[:, t], cut[1])
			@debug cut, adjusted_cut
			fire_allot_dict = Dict(
				ix => adjusted_cut[ix] for
				ix in eachindex(orig_fire_allots) if orig_fire_allots[ix] > 0
			)
			gub_cut = GUBCoverCut(
				t,
				fire_allot_dict,
				unused_crews,
				length(keys(fire_allot_dict)) + length(unused_crews) - 1,
			)
			push!(knapsack_gub_cuts, gub_cut)
		end

	end

	return knapsack_gub_cuts
end


function incorporate_gub_cover_cuts_into_fire_subproblem!(
	out_dict,
	cuts,
	fire,
	submodel,
)

	TIME_FROM_ = 3
	CREWS_PRESENT_ = 6

	for ix in eachindex(cuts)
		if ix ∉ keys(out_dict)
			cut = cuts[ix]
			costs = Dict{Int64, Int8}()
			if fire in keys(cut.fire_lower_bounds)
				for i in 1:length(submodel.arc_costs)
					if (submodel.long_arcs[i, TIME_FROM_] == cut.time_ix) & (
						submodel.long_arcs[i, CREWS_PRESENT_] >=
						cut.fire_lower_bounds[fire][1]
					)
						costs[i] = Float64(cut.fire_lower_bounds[fire][2])
					end
				end
			end
			out_dict[ix] = costs
		end
	end
	return out_dict
end


function incorporate_gub_cover_cuts_into_crew_subproblem!(
	out_dict,
	cuts,
	crew,
	submodel,
)

	for ix in eachindex(cuts)
		if ix ∉ keys(out_dict)
			cut = cuts[ix]
			costs = Dict{Int64, Int8}()
			if crew in cut.inactive_crews
				for i in 1:length(submodel.arc_costs)
					# -1 for any that do suppress at the time
					if (submodel.long_arcs[i, CM.TIME_TO] == cut.time_ix) &
					   (submodel.long_arcs[i, CM.TO_TYPE] == CM.FIRE_CODE)

						# the route needs to suppress a fire in the cut to avoid the penalty
						fires = [i for i in keys(cut.fire_lower_bounds)]
						if (submodel.long_arcs[i, CM.LOC_TO] ∈ fires)
							costs[i] = -1
						end
					end
				end
			end
			out_dict[ix] = costs
		end
	end
end

function push_cut_to_rmp!!!(cuts_per_time, cuts, cut_objects, cut, rmp)
	cut_time = cut.time_ix
	cut_fire_allots = cut.fire_lower_bounds
	cut_crew_allots = cut.inactive_crews
	cut_rhs = cut.rhs

	ix = cuts_per_time[cut_time] + 1

	lhs = []
	for f in keys(cut_fire_allots)
		for i in eachindex(rmp.plans[f, :])
			if fire_plans.crews_present[f, i..., cut_time] >= cut_fire_allots[f][1]
				push!(lhs, cut_fire_allots[f][2] * rmp.plans[f, i...])
			end
		end
	end

	# TODO allow coeffs for this cut too
	for c in cut_crew_allots
		for i in eachindex(rmp.routes[c, :])

			# crew has to suppress one of the fires in the cut to avoid entering this cut
			fires = [i for i in keys(cut_fire_allots)]
			if maximum(crew_routes.fires_fought[c, i..., fires, cut_time]) == 0
				push!(lhs, rmp.routes[c, i...])
			end
		end
	end

	cuts[cut_time, ix] = @constraint(rmp.model, -sum(lhs) >= -cut_rhs)
	cut_objects[(cut_time, ix)] = cut
	cuts_per_time[cut_time] = ix
end

function generate_and_add_knapsack_cuts_to_rmp!(
	rmp::RestrictedMasterProblem,
	crew_routes,
	fire_plans,
	crew_models,
	fire_models,
)

	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	@constraint(rmp.model, cuts[t = 1:14, u = 1:14; false], 0 >= 0)
	cut_objects = Dict()
	cuts_per_time = [0 for t ∈ 1:14]

	keep_iterating = true
	while keep_iterating
		knapsack_cuts = generate_knapsack_cuts(crew_routes, fire_plans, rmp)

		## lots of choices for how many times to introduce cuts and re-optimize
		## can look for more after reoptimizing to cut off new solution
		## but maybe don't care if the new incumbent is nowhere near the eventual solution
		## for now let's do one round of minimal cuts but could do all of them or tune based on change in objective

		keep_iterating = false
		@info length(knapsack_cuts)

		# if length(knapsack_cuts) == 0
		#     keep_iterating = false
		# end

		for cut in knapsack_cuts
			push_cut_to_rmp!!!(cuts_per_time, cuts, cut_objects, cut, rmp)
		end

		optimize!(rmp.model)
		@info objective_value(rmp.model)
	end

	fire_cut_sp_lookup = Dict{Int64, Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}}()
	for fire ∈ 1:num_fires
		mdl = fire_models[fire]
		d = Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}()
		incorporate_gub_cover_cuts_into_fire_subproblem!(d, cut_objects, fire, mdl)
		fire_cut_sp_lookup[fire] = deepcopy(d)
	end

	crew_cut_sp_lookup = Dict{Int64, Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}}()
	for crew ∈ 1:num_crews
		mdl2 = crew_models[crew]
		e = Dict{Tuple{Int64, Int64}, Dict{Int64, Int8}}()
		incorporate_gub_cover_cuts_into_crew_subproblem!(e, cut_objects, crew, mdl2)
		crew_cut_sp_lookup[crew] = deepcopy(e)
	end

	return cuts, cut_objects, cuts_per_time, fire_cut_sp_lookup, crew_cut_sp_lookup
end
