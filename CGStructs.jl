struct GlobalData
    
    ff_dist::Matrix{Float64}
    bf_dist::Matrix{Float64}
    ff_tau::Matrix{Int64}
    bf_tau::Matrix{Int64}
    
end

struct CrewStatus
    
    rest_by::Vector{Int64}
    current_fire::Vector{Int64}
    rested_periods::Vector{Int64}
    
end

struct RegionData
    
    crew_regions::Vector{Int64}
    fire_regions::Vector{Int64}
    
end


struct KeyArcIndices
    
    # fire flow data
    f_out::Array{Vector{Int64}}
    f_in::Array{Vector{Int64}}
    
    # base flow data
    b_out::Array{Vector{Int64}}
    b_in::Array{Vector{Int64}}
    
    # total crews suppressing each fire
    supp_fire::Array{Vector{Int64}}
    
    # start constraints
    start::Array{Vector{Int64}}
    
    # assignments out of region
    out_of_region::Array{Vector{Int64}}
    
end 

mutable struct RouteData
    
    routes_per_crew::Vector{Int64} # could add in length
    route_costs::Matrix{Float64}
    fires_fought::BitArray{4}
    out_of_reg::BitArray{3}
    
end

mutable struct SuppressionPlanData
    
    plans_per_fire::Vector{Int64} # could add in length
    plan_costs::Matrix{Float64}
    crews_present::Array{Int16, 3}
    
end

struct FireConstraintIndices
    
    # fire flow data
    flow_out::Array{Vector{Int64}}
    flow_in::Array{Vector{Int64}}
    
    # total crews suppressing each fire
    crews_needed::Array{Vector{Int64}}
    
    # start constraints
    start::Array{Vector{Int64}}
    
end

struct ColumnGeneration
    
    route_sps::Vector{Any}
    plan_sps::Vector{Any}
    routes::RouteData
    suppression_plans::SuppressionPlanData
    
end