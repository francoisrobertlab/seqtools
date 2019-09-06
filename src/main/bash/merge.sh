#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
. "$DIR/loadmodules.sh"

source bin/venv/bin/activate
python bin/src/main/python/MergeSampleBed.py -T 4
python bin/src/main/python/SplitBed.py -T 4 -s merge.txt
python bin/src/main/python/GenomeCoverage.py -T 4 -s merge.txt
