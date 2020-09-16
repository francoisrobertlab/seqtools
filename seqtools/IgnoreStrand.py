import logging
import os
import subprocess

import click

import seqtools.Split as sb
from seqtools.txt import Parser


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def ignorestrand(samples, index):
    '''Prepare BED file used for genome coverage on samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ignore_strand_samples(samples, index)


def ignore_strand_samples(samples='samples.txt', index=None):
    '''Prepare BED file used for genome coverage on samples.'''
    sample_names = Parser.first(samples)
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        ignore_strand_sample_splits(sample)


def ignore_strand_sample_splits(sample):
    '''Prepare BED file used for genome coverage on a single sample.'''
    print ('Prepare BED file used for genome coverage on sample {}'.format(sample))
    ignore_strand_sample(sample)
    splits = sb.splits(sample)
    for split in splits:
        ignore_strand_sample(split)


def ignore_strand_sample(sample):
    bed = sample + '.bed'
    bed_forcoverage = sample + '-forcov.bed'
    ignore_strand(bed, bed_forcoverage)


def ignore_strand(bed, output):
    '''Duplicate all annotations with opposed strand.'''
    with open(bed, 'r') as infile, open(output, 'w') as outfile:
        for line in infile:
            if line.startswith('track') or line.startswith('browser') or line.startswith('#'):
                outfile.write(line)
                continue
            columns = line.rstrip('\r\n').split('\t')
            if len(columns) >= 6:
                outfile.write('\t'.join(columns))
                outfile.write("\n")
                columns[5] = '+' if columns[5] == '-' else '-'
                outfile.write('\t'.join(columns))
                outfile.write("\n")


if __name__ == '__main__':
    ignorestrand()