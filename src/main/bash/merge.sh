#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=4
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

. $CHEC_VENV/activate
python $CHEC_PATH/MergeSampleBed.py -T 4
python $CHEC_PATH/SplitBed.py -T 4 -s merge.txt
python $CHEC_PATH/GenomeCoverage.py -T 4 -s merge.txt
