#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL
#SBATCH --output=fastqc-%A.out
#SBATCH --error=fastqc-%A.out

args=("$@")
if [ ! -z "$SLURM_CPUS_PER_TASK" ]
then
  args+=("-t" "$SLURM_CPUS_PER_TASK")
fi

fastqc "${args[@]}"
