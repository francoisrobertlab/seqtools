import logging
import os
import re

import click
import pandas as pd


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
@click.option('--binLength', '-b', type=int, default=10,
              help='Split reads in bins by their length.')
@click.option('--binMinLength', '-l', type=int, default=100,
              help='First bin minimum length.')
@click.option('--binMaxLength', '-L', type=int, default=500,
              help='Last bin maximum length.')
def main(samples, index, binlength, binminlength, binmaxlength):
    '''Split BED files from samples based on lenght of annotations.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        split_bed(sample, binlength, binminlength, binmaxlength)


def split_bed(sample, binlength, binminlength, binmaxlength):
    '''Split BED file from a single sample based on lenght of annotations.'''
    print ('Split BED file of sample {}'.format(sample))
    if binlength is not None:
        bed_raw = sample + '.bed'
        for bin_start in range(binminlength, binmaxlength, binlength):
            bin_end = min(bin_start + binlength, binmaxlength)
            sample_bin = '{}-{}-{}'.format(sample, bin_start, bin_end)
            bed_bin_raw = sample_bin + '.bed'
            print ('Splitting BED {} to BIN {}'.format(bed_raw, bed_bin_raw))
            filter_bed_by_length(bed_raw, bed_bin_raw, bin_start, bin_end)


def filter_bed_by_length(bed, output, minLength, maxLength):
    '''Filter BED file and keep only annotations that have specified size. Minimum size is included but max size is excluded.'''
    with open(bed, "r") as infile, open(output, "w") as outfile:
        for line in infile:
            if line.startswith('track') or line.startswith('browser') or line.startswith('#'):
                outfile.write(line)
                continue
            columns = line.rstrip('\r\n').split('\t')
            if len(columns) >= 3:
                start = int(columns[1])
                end = int(columns[2])
                length = end - start
                if length >= minLength and length < maxLength:
                    outfile.write(line)


if __name__ == '__main__':
    main()
