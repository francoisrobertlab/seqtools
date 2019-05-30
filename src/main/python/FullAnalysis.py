import logging

import AlignSample
import BamToBed
import DownloadSample
import FilterBam
import GenomeCoverage
import MergeSampleBed
import SplitGenomeCoverage
import click


@click.command()
@click.option('--samples', '-s', type=click.File('r'), default='samples.txt',
              help='Sample names listed one sample name by line.'
              ' The first line is ignored.'
              ' An SRR id can be provided (tab-separated) to download the FASTQ file automatically, otherwise it must be provided.')
@click.option('--merge', '-m', type=click.Path, default=None,
              help='Merge name if first columns and sample names to merge on following columns - tab delimited.')
@click.option('--fasta', '-f', type=click.Path(exists=True), default='sacCer3.fa',
              help='FASTA file used for alignment.')
@click.option('--sizes', type=click.Path(exists=True), default='sacCer3.chrom.sizes',
              help='Size of chromosomes.')
@click.option('--threads', '-t', default=1,
              help='Number of threads used to process data.')
@click.option('--splitLength', type=int, default=None,
              help='Split reads in bins by their length.')
@click.option('--splitMinLength', default=100,
              help='Split reads minimum length.')
@click.option('--splitMaxLength', default=400,
              help='Split reads maximum length.')
def main(samples, merge, fasta, sizes, threads, splitlength, splitminlength, splitmaxlength):
    '''Analyse Martin et al. data from November 2018 in Genetics.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    AlignSample.bwa_index(fasta)
    next(samples)
    samples_lines = samples.read().splitlines()
    for sample_line in samples_lines:
        if sample_line.startswith('#'):
            continue
        sample_info = sample_line.split('\t');
        sample = sample_info[0]
        srr = sample_info[1] if len(sample_info) > 1 else None
        analyse(sample, srr, fasta, sizes, threads, splitlength, splitminlength, splitmaxlength)
    if merge is not None:
        with open(merge) as infile:
            for merge_line in infile:
                if merge_line.startswith('#'):
                    continue
                merge_info = merge_line.split('\t');
                sample = merge_info[0]
                MergeSampleBed.merge_samples(sample, merge_info[1:])
                coverage(sample, sizes, splitlength, splitminlength, splitmaxlength)


def analyse(sample, srr, fasta, sizes, threads, splitlength, splitminlength, splitmaxlength):
    '''Analyse a single sample.'''
    print ('Analyse sample {}'.format(sample))
    DownloadSample.download(sample, srr)
    AlignSample.align(sample, fasta, threads)
    FilterBam.filter_bam(sample, fasta, threads)
    BamToBed.bam_to_bed(sample, threads)
    coverage(sample, sizes, splitlength, splitminlength, splitmaxlength)


def coverage(sample, sizes, splitlength, splitminlength, splitmaxlength):
    '''Compute coverage of a single sample.'''
    print ('Analyse sample {}'.format(sample))
    GenomeCoverage.genome_coverage(sample, sizes)
    if splitlength is not None:
        SplitGenomeCoverage.split_genome_coverage(sample, sizes, splitlength, splitminlength, splitmaxlength)


if __name__ == '__main__':
    main()
