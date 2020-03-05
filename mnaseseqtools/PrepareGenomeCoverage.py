import logging
import os
import subprocess

import click
import pandas as pd
import seqtools.Split as sb


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def prepgenomecov(samples, index):
    '''Prepare BED file used for genome coverage on samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        sample_splits_prepgenomecov(sample)


def sample_splits_prepgenomecov(sample):
    '''Prepare BED file used for genome coverage on a single sample.'''
    print ('Compute genome coverage on sample {}'.format(sample))
    prepgenomecov_sample(sample)
    splits = sb.splits(sample)
    for split in splits:
        prepgenomecov_sample(split)


def prepgenomecov_sample(sample):
    bed_raw = sample + '.bed'
    bed_forcoverage = sample + '-forcov.bed'
    center_annotations(bed_raw, bed_forcoverage)


def center_annotations(bed, output):
    '''Resize annotations to 1 positioned at the center.'''
    with open(bed, 'r') as infile:
        with open(output, 'w') as outfile:
            for line in infile:
                if line.startswith('track') or line.startswith('browser') or line.startswith('#'):
                    outfile.write(line)
                    continue
                columns = line.rstrip('\r\n').split('\t')
                if len(columns) >= 3:
                    start = int(columns[1])
                    end = int(columns[2])
                    length = end - start
                    start = start + int(length / 2)
                    end = start + 1
                    outfile.write(columns[0])
                    outfile.write('\t')
                    outfile.write(str(start))
                    outfile.write('\t')
                    outfile.write(str(end))
                    for i in range(3, len(columns)):
                        outfile.write('\t')
                        outfile.write(columns[i])
                    outfile.write('\n')


if __name__ == '__main__':
    prepgenomecov()
