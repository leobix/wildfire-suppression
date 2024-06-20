include("../BranchAndPrice.jl")

using JuMP, Gurobi, JSON, Profile, ArgParse, Logging, IterTools

const GRB_ENV = Gurobi.Env()


error()
# precompile
branch_and_price(3, 10, 14, max_nodes = 3, soft_heuristic_time_limit = 0.0, algo_tracking = false)

io = open("profile_logs.txt", "w")
Profile.init()
@time begin @profile branch_and_price(
	9,
	30,
	14,
	cut_loop_max=10,
	max_nodes=1,
	price_and_cut_file = "profile_timings.json",
	algo_tracking = false,
	soft_heuristic_time_limit = 0.0
)
	end

Profile.print(io, mincount = 20)
if ~Sys.iswindows()
	Profile.print(io, mincount = 20, groupby=:thread)
end

close(io)
