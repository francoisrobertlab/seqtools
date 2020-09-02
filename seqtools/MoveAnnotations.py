import logging
import os
import subprocess

import click
from seqtools.txt import Parser

import seqtools.Split as sb


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--distance', '-d', type=int,
              help='Distance - can be negative to move annotations backwards.')
@click.option('--discard-negatives/--keep-negatives', '-dn/-kn', default=True,
              help='Discard annotations that would have a negative coordinate if moved.', show_default=True)
@click.option('--reverse-negative/--forward-negative', '-rn/-fn', default=True,
              help='Remove distance for negative strand instead of adding.', show_default=True)
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def moveannotations(samples, distance, discard_negatives, reverse_negative, index):
    '''Moves annotations contained in BED files.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    moveannotations_samples(samples, distance, discard_negatives, reverse_negative, index)


def moveannotations_samples(samples='samples.txt', distance=0, discard_negatives=True, reverse_negative=True, index=None):
    '''Moves annotations contained in BED files.'''
    sample_names = Parser.first(samples)
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        moveannotations_sample(sample, distance, discard_negatives, reverse_negative)


def moveannotations_sample(sample, distance=0, discard_negatives=True, reverse_negative=True):
    '''Moves annotations contained in BED file of sample.'''
    print ('Moves annotations contained in BED file of sample {}'.format(sample))
    bed = sample + '.bed'
    moved = sample + '-forcov.bed'
    with open(bed, 'r') as infile, open(moved, 'w') as outfile:
        for line in infile:
            if line.startswith('browser') or line.startswith('track'):
                outfile.write(line)
                continue
            columns = line.rstrip('\r\n').split('\t')
            if reverse_negative and len(columns) >= 5 and columns[5] == '-':
                columns[1] = str(int(columns[1]) - distance)
                columns[2] = str(int(columns[2]) - distance)
            else:
                columns[1] = str(int(columns[1]) + distance)
                columns[2] = str(int(columns[2]) + distance)
            if discard_negatives and (columns[1].startswith("-") or columns[2].startswith("-")):
                # Discard annotation.
                logging.info("Discarding annotation {}".format(line))
                continue
            outfile.write('\t'.join(columns))
            outfile.write('\n')


if __name__ == '__main__':
    moveannotations()
