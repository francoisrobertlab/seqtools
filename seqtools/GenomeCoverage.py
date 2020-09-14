from distutils.command.check import check
import logging
import os
import subprocess
import tempfile

import click
from seqtools.bed import Bed
from seqtools.txt import Parser

from . import Split

BASE_SCALE = 1000000


@click.command(context_settings=dict(
    ignore_unknown_options=True,
))
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--genome', '-g', type=click.Path(exists=True), default='sacCer3.chrom.sizes', show_default=True,
              help='Size of chromosomes.')
@click.option('-scale', type=float, default=None,
              help='Scale for genome coverage. Defaults to 1000000 / number of reads.')
@click.option('-strand', type=click.Choice(['+', '-']), default=None, show_default=True,
              help='Calculate coverage of intervals from a specific strand.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
@click.argument('genomecov_args', nargs=-1, type=click.UNPROCESSED)
def genomecov(samples, genome, scale, strand, index, genomecov_args):
    '''Compute genome coverage on samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    genome_coverage_samples(samples, genome, scale, strand, index, genomecov_args)


def genome_coverage_samples(samples='samples.txt', genome='sacCer3.chrom.sizes', scale=None, strand=None, index=None, genomecov_args=()):
    '''Compute genome coverage on samples.'''
    sample_names = Parser.first(samples)
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        sample_splits_genome_coverage(sample, genome, scale, strand, genomecov_args)


def sample_splits_genome_coverage(sample, genome, scale=None, strand=None, genomecov_args=()):
    '''Compute genome coverage on a single sample.'''
    print ('Computing genome coverage on sample {}'.format(sample))
    genome_coverage(sample, genome, scale, strand, genomecov_args)
    splits = Split.splits(sample)
    for split in splits:
        genome_coverage(split, genome, scale, strand, genomecov_args)


def genome_coverage(sample, genome, scale=None, strand=None, genomecov_args=()):
    bed_source = sample + '-forcov.bed'
    if not os.path.exists(bed_source):
        logging.info('File {} does not exists, using {} for coverage'.format(bed_source, sample + '.bed'))
        bed_source = sample + '.bed'
    print ('Computing genome coverage on BED {}'.format(bed_source))
    if not scale:
        count = Bed.count_bed(bed_source)
        scale = BASE_SCALE / max(count, 1)
    bed = sample + '-cov.bed'
    bigwig = sample + '-cov.bw'
    if strand:
        bed = sample + '-cov' + ('-neg' if strand == '-' else '-pos') + '.bed'
        bigwig = sample + '-cov' + ('-neg' if strand == '-' else '-pos') + '.bw'
    coverage(bed_source, bed, genome, sample, scale, strand, genomecov_args)
    Bed.bedgraph_to_bigwig(bed, bigwig, genome)


def coverage(bed_input, bed_output, genome, sample, scale=None, strand=None, genomecov_args=()):
    '''Compute genome coverage.'''
    coverage_output_o, coverage_output = tempfile.mkstemp(suffix='.bed')
    cmd = ['bedtools', 'genomecov', '-bg', '-i', bed_input, '-g', genome] + list(genomecov_args)
    if scale:
        cmd.extend(['-scale', str(scale)]) 
    if strand:
        cmd.extend(['-strand', strand]) 
    logging.debug('Running {}'.format(cmd))
    with open(coverage_output_o, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)
    sort_output_o, sort_output = tempfile.mkstemp(suffix='.bed')
    Bed.sort(coverage_output, sort_output)
    os.remove(coverage_output)
    track = 'track type=bedGraph name="' + sample
    if strand:
        track += ' Minus' if strand == '-' else ' Plus'
    track += '"'
    with open(sort_output_o, 'r') as infile, open(bed_output, 'w') as outfile:
        outfile.write(track + '\n')
        outfile.writelines(infile)
    os.remove(sort_output)


if __name__ == '__main__':
    genomecov()
