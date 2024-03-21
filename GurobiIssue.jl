using JuMP, Gurobi

mutable struct MutGRBsvec
    len::Cint
    ind::Ptr{Cint}
    val::Ptr{Cdouble}
end

mutable struct mySvec
    len::Cint;
    ind::Ptr{Cint};
    val::Ptr{Cdouble};
    end

m=direct_model(Gurobi.Optimizer())
grb = backend(m)
@variable(m,y>=0)
@variable(m,x>=0)
@constraint(m,x<=1)
@constraint(m,y<=2)
@objective(m,Max,x+y)
optimize!(m)

idxs = Vector{Cint}(undef,2) # model has 2 rows
vals = Vector{Cdouble}(undef,2)
ss = mySvec(0,pointer_from_objref(idxs),pointer_from_objref(vals))
res = ccall((:GRBBinvRowi, "gurobi110"), Cint, (Ptr{GRBmodel}, Cint, Ptr{mySvec}), grb, 1, Ref(ss))

println(res)
@show ss.len
@show unsafe_load(ss.ind, 1)
@show unsafe_load(ss.val)

# num_vars_ref = Ref{Cint}()
# Gurobi.GRBgetintattr(backend(rmp.model), "NumVars", num_vars_ref)
# num_vars = num_vars_ref[]
# indices = Vector{Cint}(undef, num_vars)
# values = Vector{Cdouble}(undef, num_vars)
# svec = MutGRBsvec(num_vars, pointer_from_objref(indices), pointer_from_objref(values))

# function _info(
#     model::Gurobi.Optimizer,
#     key::MOI.ConstraintIndex{MOI.ScalarAffineFunction{Float64},<:Any},
# )
#     if haskey(model.affine_constraint_info, key.value)
#         return model.affine_constraint_info[key.value]
#     end
#     return throw(MOI.InvalidIndex(key))
# end

# link_const = Cint(_info(backend(rmp.model), index(rmp.supply_demand_linking[1, 1])).row - 1)
# model_backend = backend(rmp.model)
# GC.@preserve indices values svec model_backend begin
#     res = ccall((:GRBBinvRowi, "gurobi110"), Cint, (Ptr{GRBmodel}, Cint, Ptr{MutGRBsvec}), model_backend, link_const, Ref(svec))	
# end
# println("here")
# println(res)
# println(svec.len)
# @info "tableau" num_vars 