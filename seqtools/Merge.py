import logging
from multiprocessing import Pool
import os
import subprocess
import tempfile

import click
from seqtools.bed import Bed
from seqtools.txt import Parser


@click.command()
@click.option('--merge', '-m', type=click.Path(exists=True), default='merge.txt', show_default=True,
              help='Merge name if first columns and sample names to merge on following columns - tab delimited.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def merge(merge, index):
    '''Merge BED files related to samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    merge_samples(merge, index)


def merge_samples(merge='merge.txt', index=None):
    '''Merge BED files related to samples.'''
    merge_columns = Parser.columns(merge)
    if index != None:
        merge_columns = [merge_columns[index]]
    for columns in merge_columns:
        name = columns[0]
        samples = [sample for sample in columns[1:]]
        merge_sample(name, samples)


def merge_sample(name, samples):
    '''Merge BED files related to samples.'''
    print ('Merging samples {} into a single sample {}'.format(samples, name))
    merge_temp_o, merge_temp = tempfile.mkstemp(suffix='.bed')
    with open(merge_temp_o, 'w') as outfile:
        for sample in samples:
            sample_bed = sample + '.bed'
            with open(sample_bed, 'r') as infile:
                for line in infile:
                    if line.startswith('browser') or line.startswith('track'):
                        continue
                    outfile.write(line)
    merged_bed = name + '.bed'
    Bed.sort(merge_temp, merged_bed)
    os.remove(merge_temp)


if __name__ == '__main__':
    merge()
