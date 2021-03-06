# ChEC-seq

:information_source: *[Connecting to Compute Canada server](connect.md)*

:bulb: Most `sbatch` commands can be optimized using `--array` argument, see [sbatch](sbatch.md)

#### Steps

* [Upload dataset files to Compute Canada](#upload-dataset-files-to-compute-canada)
* [Align FASTQ files](#align-fastq-files-with-genome)
* [Filter reads](#filter-reads-to-remove-poorly-map-reads-and-duplicates)
* [Quality control check](#quality-control-check)
* [Convert BAM to BED](#convert-bam-files-to-fragment-bed-files)
* [Merge samples into dataset](#merge-dataset-samples-data)
* [Ignore fragment strand](#ignore-fragment-strand)
* [Genome converage](#genome-coverage)
* [Statistics](#statistics)
* [Heatmaps of coverage over genes versus fragment size (Optional)](#heatmaps-of-coverage-over-genes-versus-fragment-size-optional)

## Upload dataset files to Compute Canada

See [Uploading dataset files to Compute Canada server](upload.md)

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
sbatch bowtie2.sh -x sacCer3.fa.index
```

:bulb: To prevent out of memory errors, use `--array` argument for `sbatch`, see [sbatch](sbatch.md)

## Filter reads to remove poorly map reads and duplicates

```
sbatch filterbam.sh
```

:bulb: To prevent out of memory errors, use `--array` argument for `sbatch`, see [sbatch](sbatch.md)

## Quality control check

```
sbatch fastqc.sh *.bam
```

Copy the HTML and ZIP files produced by FastQC on your local computer using an FTP software and check the result. [See the documentation for FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

:bulb: Check HTML files ending with "-dedup.bam" first!

## Convert BAM files to fragment BED files

```
sbatch bam2bed.sh
```

:bulb: To prevent out of memory errors, use `--array` argument for `sbatch`, see [sbatch](sbatch.md)

## Merge dataset samples data

```
sbatch merge.sh
```

## Ignore fragment strand

```
sbatch ignorestrand.sh
sbatch ignorestrand.sh -s dataset.txt
```

:bulb: The previous commands can be called simultaneously

## Genome coverage

```
sbatch genomecov.sh -g sacCer3.chrom.sizes -is -forcov -5
sbatch genomecov.sh -s dataset.txt -g sacCer3.chrom.sizes -is -forcov -5
```

:bulb: The previous commands can be called simultaneously

## Statistics

```
sbatch statistics.sh
```

## Heatmaps of coverage over genes versus fragment size (Optional)

### Upload VAP parameters and related files

Upload your VAP parameters file and related files to Compute Canada, see [VAP](https://bitbucket.org/labjacquespe/vap/src/master/) and [Uploading dataset files to Compute Canada server](upload.md)

An example of such files for yeast is available [here](mnase-chipseq/vap/sacCer3)

For yeast, you can copy the example files using these commands

```
curl https://raw.githubusercontent.com/francoisrobertlab/seqtools/master/computecanada/mnase-chipseq/vap/sacCer3/vap_parameters.txt >> vap_parameters.txt
curl https://raw.githubusercontent.com/francoisrobertlab/seqtools/master/computecanada/mnase-chipseq/vap/sacCer3/Group_AllLongGenes_TxSorted.txt >> Group_AllLongGenes_TxSorted.txt
curl https://raw.githubusercontent.com/francoisrobertlab/seqtools/master/computecanada/mnase-chipseq/vap/sacCer3/sgdGeneAndOthers_UTR_TIF-seq_sacCer3_july_2018.txt >> sgdGeneAndOthers_UTR_TIF-seq_sacCer3_july_2018.txt
```

### Run the analysis

```
sbatch split.sh -s dataset.txt --binLength 10 --binMinLength 50 --binMaxLength 500
```

:bulb: To prevent out of memory errors when running `split.sh`, use `--array` argument for `sbatch`, see [sbatch](sbatch.md)

```
sbatch ignorestrand.sh -s dataset.txt
sbatch genomecov.sh -s dataset.txt -g sacCer3.chrom.sizes -is -forcov -5
sbatch vap.sh -s dataset.txt -p vap_parameters.txt
remove-bins.sh
```
