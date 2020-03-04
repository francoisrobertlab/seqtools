import click

from seqtools import Bam2Bed, Bowtie2, Bwa, Download, FilterBam, GenomeCoverage, Intersect, Merge, MergeBigwigs, Plot2do, SlowSplit, Split, Statistics, Vap


@click.group()
def seqtools():
    pass


seqtools.add_command(Bam2Bed.bam2bed)
seqtools.add_command(Bowtie2.bowtie2)
seqtools.add_command(Bwa.bwa)
seqtools.add_command(Download.download)
seqtools.add_command(FilterBam.filterbam)
seqtools.add_command(GenomeCoverage.genomecov)
seqtools.add_command(Intersect.intersect)
seqtools.add_command(Merge.merge)
seqtools.add_command(MergeBigwigs.mergebw)
seqtools.add_command(Plot2do.plot2do)
seqtools.add_command(SlowSplit.slowsplit)
seqtools.add_command(Split.split)
seqtools.add_command(Statistics.statistics)
seqtools.add_command(Vap.vap)

if __name__ == '__main__':
   seqtools()
