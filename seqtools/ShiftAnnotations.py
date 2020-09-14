import logging
import os
import subprocess

import click

import seqtools.Split as sb
from seqtools.txt import Parser


@click.command(context_settings=dict(
    ignore_unknown_options=True,
))
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--genome', '-g', type=click.Path(exists=True), default='sacCer3.chrom.sizes', show_default=True,
              help='Size of chromosomes.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
@click.argument('bedtools_args', nargs=-1, type=click.UNPROCESSED)
def shiftannotations(samples, genome, index, bedtools_args):
    '''Moves annotations contained in BED files.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    shift_annotations_samples(samples, genome, index, bedtools_args)


def shift_annotations_samples(samples='samples.txt', genome='sacCer3.chrom.sizes', index=None, bedtools_args=()):
    '''Moves annotations contained in BED files.'''
    sample_names = Parser.first(samples)
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        shift_annotations_sample(sample, genome, bedtools_args)


def shift_annotations_sample(sample, genome, bedtools_args=()):
    '''Moves annotations contained in BED file of sample.'''
    print ('Moves annotations contained in BED file of sample {}'.format(sample))
    bed = sample + '.bed'
    moved = sample + '-forcov.bed'
    cmd = ['bedtools', 'shift', '-i', bed, '-g', genome] + list(bedtools_args)
    logging.debug('Running {}'.format(cmd))
    with open(moved, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)


if __name__ == '__main__':
    shiftannotations()
