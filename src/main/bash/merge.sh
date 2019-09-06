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

VENV_DIR="$DIR/../../../venv"
PYTHON_DIR="$DIR/../python"

source $VENV_DIR/bin/activate
python $PYTHON_DIR/MergeSampleBed.py -T 4
python $PYTHON_DIR/SplitBed.py -T 4 -s merge.txt
python $PYTHON_DIR/GenomeCoverage.py -T 4 -s merge.txt
