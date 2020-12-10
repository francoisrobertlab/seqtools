:warning: seqtools was renamed to [robtools](https://github.com/francoisrobertlab/robtools), see https://github.com/francoisrobertlab/robtools

# seqtools : Tools to analyze NGS data

### Install
Install requirements:
* [python version 3.7.4 or newer](https://www.python.org)
* [Git](https://git-scm.com)

Install seqtools.

```
pip install git+https://git@github.com/francoisrobertlab/seqtools.git
```

### Usage

The following executables are installed:
* seqtools - tools used by most ChIP-seq analysis.
* mnasetools - tools specific to MNase-ChIP-seq analysis.
* exotools - tools specific to ChIP-exo analysis.
* chectools - tools specific to ChEC-seq analysis.

You can see the parameters details by using the `-h` parameter:

```
seqtools -h
```

### Requirements

The following are required depending on the commands used:
* [perl](https://www.perl.org)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [bwa](http://bio-bwa.sourceforge.net)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools](http://www.htslib.org)
* [bedtools](https://bedtools.readthedocs.io/en/latest/)
* [utilities from the UCSC Genome Browser](http://genome.ucsc.edu)
* [sra-toolkit from GeoArchive](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
* [vap](https://bitbucket.org/labjacquespe/vap/src/master/)
* [plot2do](https://github.com/rchereji/plot2DO)
