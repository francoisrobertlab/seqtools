from distutils.command.check import check
import logging
import os
from pathlib import Path
import subprocess

import click

from seqtools.txt import Parser


@click.command(context_settings=dict(
    ignore_unknown_options=True,
))
@click.option('--file', '-f', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
@click.argument('plot2do_args', nargs=-1, type=click.UNPROCESSED)
def plot2do(file, index, plot2do_args):
    '''Run plot2DO on samples.'''
    logging.basicConfig(filename='seqtools.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    plot2do_samples(file, index, plot2do_args)


def plot2do_samples(file, index=None, plot2do_args=()):
    '''Run plot2DO on samples.'''
    file_parent = Path(file).parent
    sample_names = Parser.first(file)
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        plot2do_sample(file_parent / sample, plot2do_args)


def plot2do_sample(sample, plot2do_args=()):
    '''Run plot2DO on a single sample.'''
    print ('Running plot2DO on sample {}'.format(sample))
    cmd = ['Rscript', 'plot2DO.R'] + list(plot2do_args)
    raw_bed = '{}.bed'.format(sample)
    if os.path.isfile(raw_bed):
        raw_cmd = cmd
        raw_cmd.extend(['-f', raw_bed])
        logging.debug('Running {}'.format(raw_cmd))
        subprocess.run(raw_cmd, check=True)


if __name__ == '__main__':
    plot2do()
