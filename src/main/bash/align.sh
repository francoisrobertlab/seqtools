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

VENV_DIR="$DIR/../../../venv"
PYTHON_DIR="$DIR/../python"

source $VENV_DIR/venv/bin/activate
python $PYTHON_DIR/AlignSample.py -T 4 -t 4
python $PYTHON_DIR/FilterBam.py -T 4 -t 4
python $PYTHON_DIR/BamToBed.py -T 4 -t 4
