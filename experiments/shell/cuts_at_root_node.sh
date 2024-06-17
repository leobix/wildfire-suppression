#!/bin/bash

#SBATCH --exclusive
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jacobwac@mit.edu

module load julia/1.10.1
module load gurobi/gurobi-1102


julia --threads 32 --project=package_dependencies/julia experiments/cuts_at_root_node.jl