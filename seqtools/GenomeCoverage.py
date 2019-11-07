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
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def main(samples, sizes, index):
    '''Compute genome coverage on samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        genome_coverage(sample, sizes)


def genome_coverage(sample, sizes):
    '''Compute genome coverage on a single sample.'''
    print ('Computing genome coverage on sample {}'.format(sample))
    do_genome_coverage(sample, sizes)
    splits = SplitBed.splits(sample)
    for split in splits:
        do_genome_coverage(split, sizes)


def do_genome_coverage(sample, sizes):
    bed_source = sample + '-forcov.bed'
    if not os.path.exists(bed_source):
        logging.info('File {} does not exists, using {} for coverage'.format(bed_source, sample + '.bed'))
        bed_source = sample + '.bed'
    print ('Computing genome coverage on BED {}'.format(bed_source))
    count = Bed.count_bed(bed_source)
    scale = BASE_SCALE / max(count, 1)
    bed = sample + '-cov.bed'
    bigwig = sample + '-cov.bw'
    coverage(bed_source, bed, sizes, sample, scale)
    bedgraph_to_bigwig(bed, bigwig, sizes)


def coverage(bed_input, bed_output, sizes, sample, scale=None, strand=None):
    '''Compute genome coverage.'''
    coverage_output = bed_input + '.cov'
    cmd = ['bedtools', 'genomecov', '-bg', '-5', '-i', bed_input, '-g', sizes]
    if not scale is None:
        cmd.extend(['-scale', str(scale)]) 
    if not strand is None:
        cmd.extend(['-strand', strand]) 
    logging.debug('Running {}'.format(cmd))
    with open(coverage_output, 'w') as outfile:
        subprocess.call(cmd, stdout=outfile)
    if not os.path.isfile(coverage_output):
        raise AssertionError('Error when computing genome coverage on ' + bed_input)
    sort_output = bed_input + '.sort'
    Bed.sort(coverage_output, sort_output)
    os.remove(coverage_output)
    track = 'track type=bedGraph name="' + sample + '"'
    if not strand is None:
        track += ' Minus' if strand == '-' else ' Plus'
    with open(sort_output, 'r') as infile, open(bed_output, 'w') as outfile:
        outfile.write(track + '\n')
        outfile.writelines(infile)
    os.remove(sort_output)


def bedgraph_to_bigwig(bed, bigwig, sizes):
    '''Converts bedgraph file to bigwig.'''
    cmd = ['bedGraphToBigWig', bed, sizes, bigwig]
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bigwig):
        raise AssertionError('Error when converting BED to BIGWIG ' + bed)


if __name__ == '__main__':
    main()
