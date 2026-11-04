#!/bin/bash
#SBATCH --partition=compute
#SBATCH --output=slurm_outputs/slurm-%A-%a.out
#SBATCH --error=slurm_errors/slurm-%A-%a.err
#SBATCH --cpus-per-task=128
#SBATCH --mem=512G
#SBATCH --time=96:00:00
#SBATCH --mail-user=d.linkevicius@oist.jp
#SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS
#SBATCH --job-name=calmodulin_mcmc
#SBATCH --array=0-14

models=(calmodulin_faas calmodulin_faas_shifman_fig1 calmodulin_full)

julia --threads="${SLURM_CPUS_PER_TASK:-1}" scripts/mcmc.jl "${models[SLURM_ARRAY_TASK_ID / 5]}" $((SLURM_ARRAY_TASK_ID % 5 + 1))
