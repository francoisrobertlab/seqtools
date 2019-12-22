from distutils.command.check import check
import logging
import os
import subprocess

import click
import pandas as pd
from seqtools.bed import Bed

from . import SplitBed

BASE_SCALE = 1000000


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--sizes', '-S', type=click.Path(exists=True), default='sacCer3.chrom.sizes',
              help='Size of chromosomes.')
@click.option('--scale', '-C', type=float, default=None,
              help='Scale for genome coverage. Defaults to 1000000 / number of reads.')
@click.option('--strand', '-T', type=click.Choice(['+', '-']), default=None,
              help='Calculate coverage of intervals from a specific strand.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def main(samples, sizes, scale, strand, index):
    '''Compute genome coverage on samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        genome_coverage(sample, sizes, scale, strand)


def genome_coverage(sample, sizes, scale=None, strand=None):
    '''Compute genome coverage on a single sample.'''
    print ('Computing genome coverage on sample {}'.format(sample))
    do_genome_coverage(sample, sizes, scale, strand)
    splits = SplitBed.splits(sample)
    for split in splits:
        do_genome_coverage(split, sizes, scale, strand)


def do_genome_coverage(sample, sizes, scale=None, strand=None):
    bed_source = sample + '-forcov.bed'
    if not os.path.exists(bed_source):
        logging.info('File {} does not exists, using {} for coverage'.format(bed_source, sample + '.bed'))
        bed_source = sample + '.bed'
    print ('Computing genome coverage on BED {}'.format(bed_source))
    count = Bed.count_bed(bed_source)
    if not scale:
        scale = BASE_SCALE / max(count, 1)
    bed = sample + '-cov.bed'
    bigwig = sample + '-cov.bw'
    if strand:
        bed = sample + '-cov' + ('-neg' if strand == '-' else '-pos') + '.bed'
        bigwig = sample + '-cov' + ('-neg' if strand == '-' else '-pos') + '.bw'
    coverage(bed_source, bed, sizes, sample, scale, strand)
    bedgraph_to_bigwig(bed, bigwig, sizes)


def coverage(bed_input, bed_output, sizes, sample, scale=None, strand=None):
    '''Compute genome coverage.'''
    coverage_output = bed_input + '.cov'
    cmd = ['bedtools', 'genomecov', '-bg', '-5', '-i', bed_input, '-g', sizes]
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


def bedgraph_to_bigwig(bed, bigwig, sizes):
    '''Converts bedgraph file to bigwig.'''
    cmd = ['bedGraphToBigWig', bed, sizes, bigwig]
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
