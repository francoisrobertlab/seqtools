#!/bin/bash
#SBATCH --account=def-robertf
#SBATCH --time=24:00:00
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=christian.poitras@ircm.qc.ca
#SBATCH --mail-type=ALL

module load python/3.7.4
module load nixpkgs/16.09 gcc/5.4.0 bwa/0.7.15
module load nixpkgs/16.09 intel/2018.3 samtools/1.9
module load nixpkgs/16.09 gcc/5.4.0 bedtools/2.27.1
source bin/venv/bin/activate
python bin/src/main/python/AlignSample.py -T 4 -t 4
python bin/src/main/python/FilterBam.py -T 4 -t 4
python bin/src/main/python/BamToBed.py -T 4 -t 4
