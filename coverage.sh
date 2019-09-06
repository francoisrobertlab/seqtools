#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

module load python/3.7.4
module load nixpkgs/16.09 gcc/5.4.0 bedtools/2.27.1
# TODO Create a module for kent-tools
module load kent-tools/1.04.00
source bin/venv/bin/activate
python bin/src/main/python/SplitBed.py -T 4
python bin/src/main/python/GenomeCoverage.py -T 4
