#!/bin/bash

echo "loading python module"
module load nixpkgs/16.09 python/3.7.4
echo "loading bwa module"
module load bwa/0.7.17
echo "loading samtools module"
module load nixpkgs/16.09 gcc/7.3.0 samtools/1.9
echo "loading bedtools module"
module load nixpkgs/16.09 gcc/7.3.0 bedtools/2.27.1
echo "loading kentutils (UCSC utilities) module"
module load nixpkgs/16.09 gcc/7.3.0 kentutils/20180716
echo "loading R module"
module load nixpkgs/16.09 gcc/7.3.0 r/3.6.0
