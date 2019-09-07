#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --array=0-0
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

. $CHEC_VENV/activate
python $CHEC_PATH/AlignSample.py -t 4 -i $SLURM_ARRAY_TASK_ID
python $CHEC_PATH/FilterBam.py -t 4 -i $SLURM_ARRAY_TASK_ID
python $CHEC_PATH/BamToBed.py -t 4 -i $SLURM_ARRAY_TASK_ID
