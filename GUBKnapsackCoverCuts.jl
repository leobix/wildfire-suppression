
include("DoubleColumnGeneration.jl")

using IterTools

function get_gub_fire_relevant_suppression_data(
	fire::Int64,
	time::Int64,
	fire_plans::FirePlanData,
	rmp::RestrictedMasterProblem,
	used_plans,
)
	sample = fire_plans.crews_present[fire, used_plans[fire], :]
	values = [value(rmp.plans[fire, i]) for i in used_plans[fire]]
	coeff_order = sortperm(sample[:, time])
	allot = 0
	weight = 0
	pairs = Tuple[]
	for coeff in coeff_order
		if sample[coeff, time] > allot
			push!(pairs, (sample[coeff, time], weight))
			allot = sample[coeff, time]
			# weight = 0
		end
		weight += values[coeff]
	end
	return pairs
end

function get_crew_suppression_cdf_by_fire_and_time(
	crew_routes::CrewRouteData,
	t::Int64,
	used_routes,
	rmp::RestrictedMasterProblem,
)
	num_crews, _, num_fires, _ = size(crew_routes.fires_fought)

	cdf = zeros(num_crews, num_fires)
	for c in 1:num_crews
		all_fires_fought = crew_routes.fires_fought[c, used_routes[c], :, t]
		coeffs = [value(rmp.routes[c, i]) for i in used_routes[c]]
		fires_fought = all_fires_fought' * coeffs
		cdf[c, :] = fires_fought
	end

	@debug "after populate" cdf

	for col in eachcol(cdf)
		sort!(col)
	end

	@debug "after sort" cdf

	cdf = cumsum(1 .- cdf, dims = 1)

	return cdf
end


function get_gub_crew_relevant_routing_data(
	c::Int64,
	t::Int64,
	crew_routes::CrewRouteData,
	used_routes,
	rmp::RestrictedMasterProblem,
)
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

function cut_generating_LP(gurobi_env,
	time_ix::Int64,
	crew_allots::Matrix{Float64},
	fire_allots::Vector{Vector{Tuple}})

	crew_usages = vec(mapslices(sum, crew_allots, dims = 2))

	inactive_crews = vec(findall(x -> x < 1e-20, crew_usages))

	# for now, assume that these crews were inactive for a "good" reason and force those to be part of any cover
	# TODO could be smarter about inactive crews, not necessarily all 0 only
	num_crews = length(crew_usages)
	available_for_fires = num_crews - length(inactive_crews)

	# define model
	m = direct_model(Gurobi.Optimizer(gurobi_env))
	@objective(m, Max, 0)
	set_optimizer_attribute(m, "OutputFlag", 0)

	# define decision variables
	vars = @variable(m, coeff[fire ∈ 1:1, allot ∈ 1:1; 0 != 0] >= 0)
	costs = []
	for fire in eachindex(fire_allots)
		specific_fire_allots = fire_allots[fire]
		push!(costs, [j[1] for j in specific_fire_allots])

		# turn survival function into CDFs
		current_weights = [j[2] for j in specific_fire_allots]
		for i in eachindex(current_weights)
			if i == length(current_weights)
				current_weights[i] = 1 - current_weights[i]
			else
				current_weights[i] = current_weights[i+1] - current_weights[i]
			end
		end
		for ix in eachindex(costs[fire])
			allot = costs[fire][ix]
			vars[fire, allot] = @variable(
				m,
				base_name = "coeff[$fire,$allot]",
				lower_bound = 0
			)
			set_objective_coefficient(
				m,
				vars[fire, allot],
				current_weights[ix],
			)
		end
	end

	# find constraints
	fires_to_consider =
		[fire for fire in eachindex(fire_allots) if length(fire_allots[fire]) > 0]
	if length(fires_to_consider) < 2
		# For now only consider cuts with >= 2 fires
		return
	end

	allotment_option_ixs = []
	for ix in eachindex(fires_to_consider)
		fire = fires_to_consider[ix]
		push!(allotment_option_ixs, [i for i in 0:length(fire_allots[fire])])
	end
	cartesian_product = product(allotment_option_ixs...)
	if length(cartesian_product) > 10000
		@info "CGLP too big, returning"
		return
	end

	for ix_per_fire in cartesian_product
		if sum(ix_per_fire) == 0
			continue
		end
		allots_per_fire =
			[
				fire_ix == 0 ? 0 : costs[fires_to_consider[ix]][fire_ix] for
				(ix, fire_ix) in enumerate(ix_per_fire)
			]
		if sum(allots_per_fire) <= available_for_fires
			a = @constraint(
				m,
				sum(
					vars[
						fires_to_consider[ix],
						costs[fires_to_consider[ix]][fire_ix],
					] for (ix, fire_ix) in enumerate(ix_per_fire) if fire_ix > 0
				) <=
				1
			)
		end
	end
	optimize!(m)
	if has_values(m) && (objective_value(m) > 1 + 1e-2)
		num_crews_in_cut = length(inactive_crews)
		crew_coeffs = [j ∈ inactive_crews ? 1.0 : 0.0 for j ∈ 1:num_crews]
		fire_coeffs = Dict{Int64, Vector{Float64}}()
		for (fire, allot) ∈ eachindex(vars)
			cost = value(vars[fire, allot])
			if cost > 1e-20
				if fire ∉ keys(fire_coeffs)
					fire_coeffs[fire] = zeros(Float64, num_crews)
				end
				for j ∈ eachindex(fire_coeffs[fire])
					if j >= allot
						fire_coeffs[fire][j] = max(cost, fire_coeffs[fire][j])
					end
				end
			end
		end
		cut = GUBKnapsackCut(
			time_ix,
			fire_coeffs,
			crew_coeffs,
			1 + num_crews_in_cut,
		)
		@debug "CGLP helped!" cut
		return cut
	end
	@debug "cglp model" length(cartesian_product)
end


function enumerate_minimal_cuts(
	crew_allots::Matrix{Float64},
	fire_allots::Vector{Vector{Tuple}},
	cut_search_limit_per_time::Int64)

	all_cuts = []

	# find the crews that are inactive at time t in all used routes
	# inactive_crews = Int64[]
	# for crew in eachindex(crew_allots)
	# 	if (length(crew_allots[crew]) > 0)
	# 		if (crew_allots[crew][1][1] == true) & (crew_allots[crew][1][2] == 0)
	# 			push!(inactive_crews, crew)
	# 		end
	# 	end
	# end

	crew_usages = vec(mapslices(sum, crew_allots, dims = 2))
	inactive_crews = vec(findall(x -> x < 1e-5, crew_usages))
	# for now, assume that these crews were inactive for a "good" reason and force those to be part of any cover
	available_for_fires = length(crew_usages) - length(inactive_crews)

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
		if loop_count >= cut_search_limit_per_time
			if loop_count > 1
				@info "Broke cover cut enumeration early" loop_count
			end
			break
		end

		possible_allots = [fire_allots[j][i[j]][1] for j in eachindex(i)]

		allot = 0
		cost = 0
		min_allot = 10000
		for j in eachindex(i)
			cur_allot = possible_allots[j]
			allot += cur_allot
			cost += fire_allots[j][i[j]][2]
			if i[j] > 1 # not the dummy (0,0) index
				min_allot = min(min_allot, cur_allot)
			end
		end
		if (cost < 1 - 1e-2) & (allot > available_for_fires) &
		   (allot - min_allot <= available_for_fires)
			push!(
				all_cuts,
				(available_for_fires, cost, possible_allots, inactive_crews),
			)
		end
	end

	return all_cuts
end

function adjust_cut_fire_allotment(
	current_allotment::Vector{Int64},
	incumbent_weighted_average::Vector{Float64},
	budget::Int64,
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

function extract_usages(
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	rmp::RestrictedMasterProblem,
)
	## TODO refactor with copied code below ##

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

	# for each time period
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
	end

	return all_fire_allots, all_crew_allots
end

function find_single_fire_knapsack_cuts(
	fire_allots,
	cdf_crew_assign,
	argsort_crew_assign,
)

	num_crews, num_fires, num_time_periods = size(cdf_crew_assign)

	cuts = []

	# get the potential cuts per fire, time, demand
	cut_violations = copy(cdf_crew_assign)
	cut_violations = 1 .- cdf_crew_assign
	for t ∈ 1:num_time_periods
		for g ∈ 1:num_fires
			for c ∈ 1:num_crews
				for (ix, (allot, survival_f)) ∈ enumerate(fire_allots[g, t])
					if allot > num_crews - c
						lifted_coeff = min(c, allot + c - num_crews)
						if ix == length(fire_allots[g, t])
							fire_demand_weight = 1 - survival_f
						else
							next_survival_f = fire_allots[g, t][ix+1][2]
							fire_demand_weight = next_survival_f - survival_f
						end
						cut_violations[c, g, t] += lifted_coeff * fire_demand_weight
					end
				end
			end
		end
	end

	# choose the max violation per fire, time
	max_per_fire = mapslices(findmax, cut_violations, dims = 1)
	for t ∈ 1:num_time_periods
		for g ∈ 1:num_fires
			if max_per_fire[1, g, t][1] > 1.01
				n_crews_in_cut = max_per_fire[1, g, t][2]
				crews_present = argsort_crew_assign[1:n_crews_in_cut, g, t]
				cut_allotments =
					[
						ix == g ? num_crews - n_crews_in_cut + 1 : 0 for
						ix ∈ 1:num_fires
					]
				push!(
					cuts,
					(
						num_crews - n_crews_in_cut,
						0,
						cut_allotments,
						crews_present,
						t,
					),
				)
			end
		end
	end

	@info "cut_violations" cdf_crew_assign cut_violations cuts


	return cuts
end
function find_knapsack_cuts(
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	rmp::RestrictedMasterProblem,
	cut_search_limit_per_time::Int64,
	general_cuts::Bool,
	single_fire_cuts::Bool,
)

	# get problem dimensions
	num_crews, _, num_fires, num_time_periods = size(crew_routes.fires_fought)

	all_fire_allots = Array{Vector{Tuple}, 2}(undef, num_fires, num_time_periods)

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

	crew_assign = get_crew_incumbent_weighted_average(rmp, crew_routes)

	# intitialize output 
	knapsack_gub_cuts = []

	# for each time period
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
	end

	if single_fire_cuts
		# crew_assign = get_crew_incumbent_weighted_average(rmp, crew_routes)
		# @info "crew_assign" crew_assign
		# cdf_crew_assign =
		# 	mapslices(cumsum, mapslices(sort, crew_assign, dims = 1), dims = 1)
		# argsort_crew_assign = mapslices(sortperm, crew_assign, dims = 1)
		# cuts = find_single_fire_knapsack_cuts(
		# 	all_fire_allots,
		# 	cdf_crew_assign,
		# 	argsort_crew_assign,
		# )
		# for cut in cuts
		# 	orig_fire_allots = cut[3]
		# 	unused_crews = cut[4]
		# 	fire_allot_dict = Dict(
		# 		ix => orig_fire_allots[ix] for
		# 		ix in eachindex(orig_fire_allots) if orig_fire_allots[ix] > 0
		# 	)
		# 	num_crews_in_cut = length(unused_crews)

		# 	gub_cut = GUBCoverCut(
		# 		cut[5],
		# 		fire_allot_dict,
		# 		unused_crews,
		# 		length(keys(fire_allot_dict)) + num_crews_in_cut - 1,
		# 	)
		# 	push!(knapsack_gub_cuts, gub_cut)
		# end

		error("are you sure you want to do this, this did not help in testing")
	end

	if general_cuts
		for t in 1:num_time_periods

			minimal_cuts =
				enumerate_minimal_cuts(
					crew_assign[:, :, t],
					all_fire_allots[:, t],
					cut_search_limit_per_time,
				)

			for cut in minimal_cuts
				orig_fire_allots = cut[3]
				unused_crews = cut[4]
				adjusted_cut =
					adjust_cut_fire_allotment(
						orig_fire_allots,
						fire_alloc[:, t],
						cut[1],
					)
				fire_allot_dict = Dict(
					ix => adjusted_cut[ix] for
					ix in eachindex(orig_fire_allots) if orig_fire_allots[ix] > 0
				)
				num_crews_in_cut = length(unused_crews)

				crew_coeffs = [j ∈ unused_crews ? 1.0 : 0.0 for j ∈ 1:num_crews]
				fire_coeffs = Dict(
					fire => zeros(Float64, num_crews) for
					fire ∈ keys(fire_allot_dict)
				)
				for fire in keys(fire_coeffs)
					for j ∈ 1:num_crews
						if j >= fire_allot_dict[fire]
							fire_coeffs[fire][j] = 1.0

							# lift coefficients for fires in special case of 1 fire involved only
							if length(keys(fire_coeffs)) == 1

								# based on the number of crews involved in the cut, we can lift
								fire_coeffs[fire][j] +=
									min(j - fire_allot_dict[fire], num_crews_in_cut)
							end
						end
					end
				end

				cut = GUBKnapsackCut(
					t,
					fire_coeffs,
					crew_coeffs,
					length(keys(fire_allot_dict)) + num_crews_in_cut - 1,
				)

				push!(knapsack_gub_cuts, cut)
			end
		end
	end


	for t in 1:num_time_periods
		if t ∉ [knapsack_cut.time_ix for knapsack_cut ∈ knapsack_gub_cuts]
			max_viol_cut = cut_generating_LP(GRB_ENV, t, crew_assign[:, :, t], all_fire_allots[:, t])
			if ~isnothing(max_viol_cut)

				# TODO all the domination not needed if we keep this outer if clause
				incorporate_cut = true
				for cut ∈ knapsack_gub_cuts
					if first_cut_dominates(cut, max_viol_cut)
						@debug "cut dominated" cut max_viol_cut
						incorporate_cut = false
						break
					end
				end
				if incorporate_cut
					push!(knapsack_gub_cuts, max_viol_cut)
					@debug "pushing max viol cglp cut!" max_viol_cut
				end
			end
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

	# for each cut ix
	for ix in eachindex(cuts)

		# if the cut ix is not in the cut dict
		if ix ∉ keys(out_dict)

			# grab the cut and initialize a lookup dictionary arc_ix => cut coeff
			cut = cuts[ix]
			costs = Dict{Int64, Float64}()

			# if the fire is involved in the cut 
			if fire in keys(cut.fire_coeffs)

				# for each arc in the subproblem
				for i in 1:length(submodel.arc_costs)

					# if it considers the given time
					if submodel.long_arcs[i, FM.TIME_FROM] == cut.time_ix

						# get cost based on the coefficient
						crews_present = submodel.long_arcs[
							i,
							FM.CREWS_PRESENT,
						]
						if (crews_present > 0) 
							cost = cut.fire_coeffs[fire][crews_present]
							if cost > 1e-20
								costs[i] = cost
							end
						end
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
			costs = Dict{Int64, Float64}()
			if cut.crew_coeffs[crew] > 1e-20
				ixs = findall(submodel.long_arcs[:, CM.CREW_NUMBER] .== crew)
				for i in ixs
					# - coeff for any that do suppress at the time
					if (submodel.long_arcs[i, CM.TIME_TO] == cut.time_ix) &
					   (submodel.long_arcs[i, CM.TO_TYPE] == CM.FIRE_CODE)

						# the route needs to suppress a fire in the cut to avoid the penalty
						fires = [i for i in keys(cut.fire_coeffs)]
						if (submodel.long_arcs[i, CM.LOC_TO] ∈ fires)
							costs[i] = -cut.crew_coeffs[crew]
						end
					end
				end
			end
			out_dict[ix] = costs
		end
	end
end

function update_cut_fire_mp_lookup!(
	cut_fire_mp_lookup::Dict{Tuple{Int64, Int64}, Float64},
	cut::GUBKnapsackCut,
	fire_plans::FirePlanData,
	fire::Int64,
	plan_ix::Int64,
)
	crews_present = fire_plans.crews_present[
		fire,
		plan_ix...,
		cut.time_ix,
	]
	cut_coeff = 0.0
	if crews_present > 0
		cut_coeff = cut.fire_coeffs[fire][crews_present]
	end
	if cut_coeff > 1e-20
		cut_fire_mp_lookup[(fire, plan_ix)] = cut_coeff
		return true
	end
	return false

end

function push_cut_to_rmp!!(
	rmp::RestrictedMasterProblem,
	cut_data::GUBCutData,
	cut::GUBKnapsackCut,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
)

	cut_time = cut.time_ix
	cut_rhs = cut.rhs

	ix = cut_data.cuts_per_time[cut_time] + 1

	cut_fire_mp_lookup = Dict{Tuple{Int64, Int64}, Float64}()
	for f in keys(cut.fire_coeffs)
		for i in eachindex(rmp.plans[f, :])
			update_cut_fire_mp_lookup!(cut_fire_mp_lookup, cut, fire_plans, f, i...)
		end
	end

	cut_crew_mp_lookup = Dict{Tuple{Int64, Int64}, Any}()
	for c in eachindex(cut.crew_coeffs)
		if cut.crew_coeffs[c] > 1e-20
			for i in eachindex(rmp.routes[c, :])

				# crew has to suppress one of the fires in the cut to avoid entering this cut
				fires = [i for i in keys(cut.fire_coeffs)]
				if maximum(crew_routes.fires_fought[c, i..., fires, cut_time]) == 0
					cut_crew_mp_lookup[(c, i...)] = cut.crew_coeffs[c]
				end
			end
		end
	end

	# TODO this method of forming JuMP expressions tosses a warning
	# Warning: The addition operator has been used on JuMP expressions a large number of times. 
	# This warning is safe to ignore but may indicate that model generation 
	# is slower than necessary. For performance reasons, you should not 
	# add expressions in a loop. Instead of x += y, use add_to_expression!(x,y) 
	# to modify x in place. If y is a single variable, you may also use 
	# add_to_expression!(x, coef, y) for x += coef*y.
	# But low priority to fix, not the bottleneck

	lhs = []
	for ix in eachindex(cut_fire_mp_lookup)
		push!(lhs, cut_fire_mp_lookup[ix] * rmp.plans[ix])
	end

	for ix in eachindex(cut_crew_mp_lookup)
		push!(lhs, cut_crew_mp_lookup[ix] * rmp.routes[ix])
	end

	# update data structures
	rmp.gub_cover_cuts[cut_time, ix] = @constraint(rmp.model, -sum(lhs) >= -cut_rhs)
	cut_data.cut_dict[(cut_time, ix)] = cut
	cut_data.cuts_per_time[cut_time] = ix
	cut_data.crew_mp_lookup[(cut_time, ix)] = cut_crew_mp_lookup
	cut_data.fire_mp_lookup[(cut_time, ix)] = cut_fire_mp_lookup
end

function find_and_incorporate_knapsack_gub_cuts!!(
	gub_cut_data::GUBCutData,
	cut_search_limit_per_time::Int64,
	rmp::RestrictedMasterProblem,
	crew_routes::CrewRouteData,
	fire_plans::FirePlanData,
	crew_models,
	fire_models;
	general_cuts = true,
	single_fire_cuts = false,
)

	num_crews, _, num_fires, _ = size(crew_routes.fires_fought)
	num_cuts_found = 0
	keep_iterating = true
	while keep_iterating
		knapsack_cuts = find_knapsack_cuts(
			crew_routes,
			fire_plans,
			rmp,
			cut_search_limit_per_time,
			general_cuts,
			single_fire_cuts,
		)
		num_cuts_found += length(knapsack_cuts)

		## lots of choices for how many times to introduce cuts and re-optimize
		## can look for more after reoptimizing to cut off new solution
		## but maybe don't care if the new incumbent is nowhere near the eventual solution
		## for now let's do one round of minimal cuts but could do all of them or tune based on change in objective

		keep_iterating = false

		# if length(knapsack_cuts) == 0
		#     keep_iterating = false
		# end

		for cut in knapsack_cuts
			push_cut_to_rmp!!(
				rmp,
				gub_cut_data,
				cut,
				crew_routes,
				fire_plans,
			)
		end

		optimize!(rmp.model)
	end

	for fire ∈ 1:num_fires
		incorporate_gub_cover_cuts_into_fire_subproblem!(
			gub_cut_data.fire_sp_lookup[fire],
			gub_cut_data.cut_dict,
			fire,
			fire_models[fire],
		)
	end

	for crew ∈ 1:num_crews
		incorporate_gub_cover_cuts_into_crew_subproblem!(
			gub_cut_data.crew_sp_lookup[crew],
			gub_cut_data.cut_dict,
			crew,
			crew_models[crew],
		)
	end

	return num_cuts_found
end
