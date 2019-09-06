#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

DIR="${BASH_SOURCE%/*}"
if [[ ! -d "$DIR" ]]; then DIR="$PWD"; fi
. "$DIR/loadmodules.sh"

VENV_DIR="$DIR/../../../venv"
PYTHON_DIR="$DIR/../python"

source $VENV_DIR/venv/bin/activate
python $PYTHON_DIR/Plot2do.py $@
