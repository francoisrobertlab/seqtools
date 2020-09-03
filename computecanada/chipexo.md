# ChIP-exo

:information_source: *[Connecting to Compute Canada server](connect.md)*

:bulb: Most `sbatch` commands can be optimized using `--array` argument, see [sbatch](sbatch.md)

#### Steps

* [Upload dataset files to Compute Canada](upload)
* [Align FASTQ files](align)
* [Filter reads](filter)
* [Remove second mate](removesecondmate)
* [Convert BAM to BED](bam2bed)
* [Merge samples into dataset](merge)
* [Move annotations](moveannotations)
* [Genome converage](genomecov)
* [Statistics](statistics)

<a name="upload"/>

## Upload dataset files to Compute Canada

See [Uploading dataset files to Compute Canada server](upload.md)

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

### Run bwa

Run the following commands

```
bwa index sacCer3.fa
sbatch bwa.sh --fasta sacCer3.fa
```

:bulb: To prevent out of memory errors, use `--array` argument for `sbatch`, see [sbatch](sbatch.md)

<a name="filter"/>

## Filter reads to remove poorly map reads and duplicates

```
sbatch filterbam.sh
```

<a name="removesecondmate"/>

## Remove second mate

```
sbatch removesecondmate.sh
```

<a name="bam2bed"/>

## Convert BAM files to fragment BED files

```
sbatch bam2bed.sh --unpaired --suffix -mate1
```

<a name="merge"/>

## Merge dataset samples data

```
sbatch merge.sh -m dataset.txt
```

<a name="moveannotations"/>

## Move annotations

```
sbatch moveannotations.sh -d 6 -rn -dn
sbatch moveannotations.sh -s dataset.txt -d 6 -rn -dn
```

:bulb: These commands can be called simultaneously

<a name="genomecov"/>

## Genome coverage

```
sbatch genomecov.sh -S sacCer3.chrom.sizes -5
sbatch genomecov.sh -S sacCer3.chrom.sizes -5 --strand +
sbatch genomecov.sh -S sacCer3.chrom.sizes -5 --strand -
sbatch genomecov.sh -s dataset.txt -S sacCer3.chrom.sizes -5
sbatch genomecov.sh -s dataset.txt -S sacCer3.chrom.sizes -5 --strand +
sbatch genomecov.sh -s dataset.txt -S sacCer3.chrom.sizes -5 --strand -
```

:bulb: These commands can be called simultaneously

<a name="statistics"/>

## Statistics

```
sbatch statistics.sh -m dataset.txt
```
