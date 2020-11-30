import click

from seqtools import Bam2Bed, Bowtie2, Bwa, CenterAnnotations, ChipexoQual, Download, FilterBam, Fixmd5, GenomeCoverage, IgnoreStrand, Intersect, IntersectAnnotations, Merge, MergeBam, MergeBigwigs, Plot2do, RemoveSecondMate, Rename, ShiftAnnotations, SlowSplit, Split, Statistics, Vap


@click.group()
def seqtools():
    pass


seqtools.add_command(Bam2Bed.bam2bed)
seqtools.add_command(Bowtie2.bowtie2)
seqtools.add_command(Bwa.bwa)
seqtools.add_command(CenterAnnotations.centerannotations)
seqtools.add_command(ChipexoQual.chipexoqual)
seqtools.add_command(Download.download)
seqtools.add_command(FilterBam.filterbam)
seqtools.add_command(Fixmd5.fixmd5)
seqtools.add_command(GenomeCoverage.genomecov)
seqtools.add_command(IgnoreStrand.ignorestrand)
seqtools.add_command(Intersect.intersect)
seqtools.add_command(IntersectAnnotations.intersectannotations)
seqtools.add_command(Merge.merge)
seqtools.add_command(MergeBam.mergebam)
seqtools.add_command(MergeBigwigs.mergebw)
seqtools.add_command(Plot2do.plot2do)
seqtools.add_command(RemoveSecondMate.removesecondmate)
seqtools.add_command(Rename.rename)
seqtools.add_command(ShiftAnnotations.shiftannotations)
seqtools.add_command(SlowSplit.slowsplit)
seqtools.add_command(Split.split)
seqtools.add_command(Statistics.statistics)
seqtools.add_command(Vap.vap)

if __name__ == '__main__':
   seqtools()
