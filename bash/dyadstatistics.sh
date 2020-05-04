#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --mem=16G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL
#SBATCH --output=dyadcoverage-%A.out
#SBATCH --error=dyadcoverage-%A.out

args+=("$@")
if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
  args+=("-i" "$SLURM_ARRAY_TASK_ID")
fi

mnasetools dyadstatistics "${args[@]}"
