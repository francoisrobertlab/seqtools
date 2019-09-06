#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

module load nixpkgs/16.09 python/3.7.4
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0

source bin/venv/bin/activate
python bin/src/main/python/Plot2do.py $@
