#!/bin/bash
#SBATCH -a 1-2
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=sched_mit_sloan_batch
#SBATCH --time=0-08:00
#SBATCH -o /home/jacobwac/ihc_routing/data/experimental_outputs/20221123/output_\%a.out
#SBATCH -e /home/jacobwac/ihc_routing/data/experimental_outputs/20221123/error_\%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jacobwac@mit.edu

module load julia/1.7.1
module load gurobi/8.1.1


julia cluster_client.jl "20221123" $SLURM_ARRAY_TASK_ID
