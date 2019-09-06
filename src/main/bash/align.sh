#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
. "$DIR/loadmodules.sh"

source bin/venv/bin/activate
python bin/src/main/python/AlignSample.py -T 4 -t 4
python bin/src/main/python/FilterBam.py -T 4 -t 4
python bin/src/main/python/BamToBed.py -T 4 -t 4
