#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --array=0-0
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

if [ -z "$SLURM_ARRAY_TASK_ID" ]
then
  SLURM_ARRAY_TASK_ID=0
fi

python $CHEC_PATH/MergeSampleBed.py -i $SLURM_ARRAY_TASK_ID
python $CHEC_PATH/SplitBed.py -s merge.txt -i $SLURM_ARRAY_TASK_ID
python $CHEC_PATH/GenomeCoverage.py -s merge.txt -i $SLURM_ARRAY_TASK_ID
