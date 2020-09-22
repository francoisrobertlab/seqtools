#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --mem=16G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL
#SBATCH --output=dyadstatistics-%A.out
#SBATCH --error=dyadstatistics-%A.out

args=("$@")

mnasetools dyadstatistics "${args[@]}"
