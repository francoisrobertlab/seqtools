import click

from seqtools import Bam2Bed, Bowtie2, Bwa, DownloadSample, FilterBam, GenomeCoverage, IntersectBed, MergeBed, MergeBigwigs, Plot2do, SlowSplitBed, SplitBed, StatisticsFile, Vap


@click.group()
def seqtools():
    pass


seqtools.add_command(Bam2Bed.bam2bed)
seqtools.add_command(Bowtie2.bowtie2)
seqtools.add_command(Bwa.bwa)
seqtools.add_command(DownloadSample.download)
seqtools.add_command(FilterBam.filterbam)
seqtools.add_command(GenomeCoverage.genomecov)
seqtools.add_command(IntersectBed.intersect)
seqtools.add_command(MergeBed.merge)
seqtools.add_command(MergeBigwigs.mergebw)
seqtools.add_command(Plot2do.plot2do)
seqtools.add_command(SlowSplitBed.slowsplit)
seqtools.add_command(SplitBed.split)
seqtools.add_command(StatisticsFile.statistics)
seqtools.add_command(Vap.vap)

if __name__ == '__main__':
   seqtools()
