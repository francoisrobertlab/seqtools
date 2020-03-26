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
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
@click.argument('ma_args', nargs=-1, type=click.UNPROCESSED)
def moveannotation(samples, index, ma_args):
    '''Moves annotations contained in BED files.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    moveannotation_samples(samples, index, ma_args)


def moveannotation_samples(samples='samples.txt', index=None, ma_args=()):
    '''Moves annotations contained in BED files.'''
    sample_names = Parser.first(samples)
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        moveannotation_sample(sample, ma_args)


def moveannotation_sample(sample, ma_args=()):
    '''Moves annotations contained in BED file of sample.'''
    print ('Moves annotations contained in BED file of sample {}'.format(sample))
    bed = sample + '.bed'
    moved = sample + '-forcov.bed'
    cmd = ['java', '-jar', 'bed-tools-j-2.1.jar', 'moveannotations'] + list(ma_args)
    cmd.extend(['-i', bed, '-o', moved])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    removesecondmate()
