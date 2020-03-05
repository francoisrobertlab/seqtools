#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --array=0-0
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL
#SBATCH --output=firstdyadcoverage-%A_%a.out
#SBATCH --error=firstdyadcoverage-%A_%a.out

args+=("$@")
if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
  args+=("-i" "$SLURM_ARRAY_TASK_ID")
fi

mnasetools dyadcov "${args[@]}"
