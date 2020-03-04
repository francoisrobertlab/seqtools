from distutils.command.check import check
import logging
import os
import subprocess

import click
import pandas as pd
from seqtools.bed import Bed

from . import Split

BASE_SCALE = 1000000


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--sizes', '-S', type=click.Path(exists=True), default='sacCer3.chrom.sizes', show_default=True,
              help='Size of chromosomes.')
@click.option('--five', is_flag=True, default=False,
              help='Calculate coverage of 5’ positions (instead of entire interval).')
@click.option('--three', is_flag=True, default=False,
              help='Calculate coverage of 3’ positions (instead of entire interval).')
@click.option('--scale', '-C', type=float, default=None,
              help='Scale for genome coverage. Defaults to 1000000 / number of reads.')
@click.option('--strand', '-T', type=click.Choice(['+', '-']), default=None, show_default=True,
              help='Calculate coverage of intervals from a specific strand.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def genomecov(samples, sizes, five, three, scale, strand, index):
    '''Compute genome coverage on samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        sample_splits_genome_coverage(sample, sizes, five, three, scale, strand)


def sample_splits_genome_coverage(sample, sizes, five=False, three=False, scale=None, strand=None):
    '''Compute genome coverage on a single sample.'''
    print ('Computing genome coverage on sample {}'.format(sample))
    genome_coverage(sample, sizes, five, three, scale, strand)
    splits = Split.splits(sample)
    for split in splits:
        genome_coverage(split, sizes, five, three, scale, strand)


def genome_coverage(sample, sizes, five=False, three=False, scale=None, strand=None):
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
    coverage(bed_source, bed, sizes, sample, five, three, scale, strand)
    Bed.bedgraph_to_bigwig(bed, bigwig, sizes)


def coverage(bed_input, bed_output, sizes, sample, five=False, three=False, scale=None, strand=None):
    '''Compute genome coverage.'''
    coverage_output = bed_input + '.cov'
    cmd = ['bedtools', 'genomecov', '-bg', '-i', bed_input, '-g', sizes]
    if five:
        cmd.extend(['-5']) 
    if three:
        cmd.extend(['-3']) 
    if scale:
        cmd.extend(['-scale', str(scale)]) 
    if strand:
        cmd.extend(['-strand', strand]) 
    logging.debug('Running {}'.format(cmd))
    with open(coverage_output, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)
    sort_output = bed_input + '.sort'
    Bed.sort(coverage_output, sort_output)
    os.remove(coverage_output)
    track = 'track type=bedGraph name="' + sample
    if strand:
        track += ' Minus' if strand == '-' else ' Plus'
    track += '"'
    with open(sort_output, 'r') as infile, open(bed_output, 'w') as outfile:
        outfile.write(track + '\n')
        outfile.writelines(infile)
    os.remove(sort_output)


if __name__ == '__main__':
    genomecov()
