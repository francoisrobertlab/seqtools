#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL
#SBATCH --output=statistics-%A.out
#SBATCH --error=statistics-%A.out

args=("$@")

statistics "${args[@]}"
