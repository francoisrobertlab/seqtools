#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

module load python/3.7.4
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
source ~/projects/def-robertf/mnase-201908/venv/bin/activate
python ~/projects/def-robertf/mnase-201908/src/main/python/Plot2do.py $@
