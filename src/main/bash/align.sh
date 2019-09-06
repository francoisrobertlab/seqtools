#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

. $CHEC_VENV/activate
python $CHEC_PATH/AlignSample.py -T 4 -t 4
python $CHEC_PATH/FilterBam.py -T 4 -t 4
python $CHEC_PATH/BamToBed.py -T 4 -t 4
