include("../BranchAndPrice.jl")

using JuMP, Gurobi, JSON, Profile, ArgParse, Logging, IterTools

const GRB_ENV = Gurobi.Env()



# precompile
branch_and_price(6, 20, 14, max_nodes = 3, soft_heuristic_time_limit = 0.0, algo_tracking = false)

io = open("profile_logs.txt", "w")
Profile.init()
@profile branch_and_price(
	9,
	30,
	14,
	max_nodes = 1,
	price_and_cut_file = "profile_timings.json",
	algo_tracking = false,
	soft_heuristic_time_limit = 0.0
)
Profile.print(io, mincount = 100)
close(io)
