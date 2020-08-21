# ChIP-seq

*All commands assume sacCer3 genome with [2 samples](sbatch.md)*

*[Connecting to Compute Canada server](connect.md)*

#### Steps

* [Upload dataset files to Compute Canada](upload.md)
* [Align FASTQ files](align)
* [Filter reads](filter)
* [Convert BAM to BED](bam2bed)
* [Merge samples into dataset](merge)
* [Genome converage](genomecov)
* [Statistics](statistics)

<a name="align"/>

## Align FASTQ files with genome

### Download the FASTA file of the genome and chromosomes size

1. Go to the [UCSC Genome Browser downloads](http://hgdownload.soe.ucsc.edu/downloads.html)
2. Select organism
3. Select "*Genome sequence files and select annotations*"
4. Download the *2bit* file of the organism like *sacCer3.2bit* for Yeast
4. Download the *chrom.sizes* file of the organism like *sacCer3.chrom.sizes* for Yeast
5. [Upload the files to Compute Canada server in the same folder of the dataset files](upload.md)
6. Convert the *2bit* file into a FASTA file using the following command

```
twoBitToFa sacCer3.2bit sacCer3.fa
```

### Run bowtie2

Run the following commands

```
bowtie2-build sacCer3.fa sacCer3.fa.index
sbatch --array=0-1 bowtie2.sh -x sacCer3.fa.index
```

<a name="filter"/>

## Filter reads to remove poorly map reads and duplicates

```
sbatch --array=0-1 filterbam.sh
```

<a name="bam2bed"/>

## Convert BAM files to fragment BED files

```
sbatch --array=0-1 bam2bed.sh
```

<a name="merge"/>

## Merge dataset samples data

```
sbatch --array=0-1 merge.sh -s dataset.txt
```

<a name="genomecov"/>

## Genome coverage

```
sbatch --array=0-1 genomecov.sh -S sacCer3.chrom.sizes
sbatch --array=0 genomecov.sh -s dataset.txt -S sacCer3.chrom.sizes
```

<a name="statistics"/>

## Statistics

```
sbatch statistics.sh -m dataset.txt
```
