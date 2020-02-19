#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --array=0-0
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL
#SBATCH --output=plot2do-%A_%a.out
#SBATCH --error=plot2do-%A_%a.out

args=("$@")
if [ ! -z "$SLURM_ARRAY_TASK_ID" ]
then
  args+=("-i" "$SLURM_ARRAY_TASK_ID")
fi

# Parameters for Martin.
# -t dyads -r Plus1 -L 400 -m 0.02
plot2do "${args[@]}"
