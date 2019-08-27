import logging
import os

import AlignSample
import BamToBed
import DownloadSample
import FilterBam
import GenomeCoverage
import MergeSampleBed
import SplitGenomeCoverage
import click


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.'
              ' An SRR id can be provided (tab-separated) to download the FASTQ file automatically, otherwise FASTQ file must be provided.')
@click.option('--merge', '-m', type=click.Path(), default='merge.txt',
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
@click.option('--splitMaxLength', default=500,
              help='Split reads maximum length.')
def main(samples, merge, fasta, sizes, threads, splitlength, splitminlength, splitmaxlength):
    '''Analyse Martin et al. data from November 2018 in Genetics.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    AlignSample.bwa_index(fasta)
    samples_columns = all_columns(samples)
    for sample_columns in samples_columns:
        sample = sample_columns[0]
        fastq = sample_columns[1] if len(sample_columns) > 1 else None
        srr = sample_columns[2] if len(sample_columns) > 2 else None
        analyse(sample, fastq, srr, fasta, sizes, threads, splitlength, splitminlength, splitmaxlength)
    merges_columns = all_columns(merge) if os.path.isfile(merge) else []
    for merge_columns in merges_columns:
        sample = merge_columns[0]
        samples_to_merge = merge_columns[1:] if len(merge_columns) > 1 else None
        MergeSampleBed.merge_samples(sample, samples_to_merge)
        coverage(sample, sizes, splitlength, splitminlength, splitmaxlength)


def all_columns(file):
    all = []
    with open(file, 'r') as lines:
        for line in lines:
            if line.startswith('#'):
                continue
            columns = line.rstrip('\n\r').split('\t');
            all.extend([columns])
    return all
    
    
def first_column(file):
    all = all_columns(file)
    return [columns[0] for columns in all]


def analyse(sample, fastq, srr, fasta, sizes, threads, splitlength, splitminlength, splitmaxlength):
    '''Analyse a single sample.'''
    print ('Analyse sample {}'.format(sample))
    DownloadSample.download(sample, fastq, srr)
    AlignSample.align(sample, fastq, fasta, threads)
    FilterBam.filter_bam(sample, threads)
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
