#!/bin/bash
#SBATCH --job-name=calmodulin_mcmc
#SBATCH --array=0-14

models=(calmodulin_faas calmodulin_faas_shifman_fig1 calmodulin_full)

julia --threads="${SLURM_CPUS_PER_TASK:-1}" scripts/mcmc.jl "${models[SLURM_ARRAY_TASK_ID / 5]}" $((SLURM_ARRAY_TASK_ID % 5 + 1))
