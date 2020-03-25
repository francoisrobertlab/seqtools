import logging
import os
import subprocess

import click
import seqtools.Split as sb
from seqtools.txt import Parser


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--threads', '-t', default=1, show_default=True,
              help='Number of threads used to process data per sample.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def removesecondmate(samples, threads, index):
    '''Removes second mate from BAM.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    prep_genomecov(samples, threads, index)


def removesecondmate_samples(samples='samples.txt', threads=None, index=None):
    '''Removes second mate from BAM.'''
    sample_names = Parser.first(samples)
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        sample_splits_prepgenomecov(sample, threads)


def removesecondmate_sample(sample, threads=None):
    '''Removes second mate from BAM for a single sample.'''
    print ('Removes second mate from BAM for sample {}'.format(sample))
    bam = sample + '.bam'
    mate1_bam = sample + '-mate1.bam'
    cmd = ['samtools', 'view']
    if not threads is None and threads > 1:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-f', '64', '-b', '-o', mate1_bam, bam])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    removesecondmate()
