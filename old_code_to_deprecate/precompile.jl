using PackageCompiler, Libdl

PackageCompiler.create_sysimage(
    ["JuMP", "Gurobi"],
    sysimage_path = "project_image." * Libdl.dlext,
    precompile_execution_file = "Main.jl",
)
