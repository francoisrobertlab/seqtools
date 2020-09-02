from distutils.command.check import check
import logging
import os
import subprocess

import click
from seqtools.bed import Bed

import pandas as pd


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples-filter.txt', show_default=True,
              help='Tab-separated file containing '
              'the sample names with tag in the first column and '
              'the original sample names in the second column.')
@click.option('--annotations', '-a', type=click.Path(exists=True), default='annotations.bed', show_default=True,
              help='Keep reads for which their center is located on specified annotations.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def intersect(samples, annotations, index):
    '''Keep only reads that intersects specified annotations.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    intersect_samples(samples, annotations, index)


def intersect_samples(samples='samples-filter.txt', annotations='annotations.bed', index=None):
    '''Keep only reads that intersects specified annotations.'''
    annot_length = annotations_length(annotations)
    sample_columns = pd.read_csv(samples, header=None, sep='\t', comment='#')
    if index != None:
        sample_columns = sample_columns.iloc[index:index + 1]
    for index, columns in sample_columns.iterrows():
        tag = columns[0]
        sample = columns[1] if len(columns) > 1 else None
        intersect_sample(sample, tag, annotations, annot_length)


def annotations_length(annotations):
    '''Find the length of an annotation'''
    with open(annotations, 'r') as infile:
        for line in infile:
            if not line.startswith('#'):
                return len(line.split('\t'))


def intersect_sample(sample, tag, annotations, annot_length):
    '''Keep only reads that intersects specified annotations for a single sample.'''
    print ('Keep only reads that intersects specified annotations for sample {}'.format(sample))
    bed = sample + '.bed'
    bed_tag = tag + '.bed'
    bed_intersect_tmp = tag + '-tosort.bed'
    bed_sort_tmp = tag + '-tmp.bed'
    cmd = ['bedtools', 'intersect', '-a', annotations, '-b', bed, '-wb']
    logging.debug('Running {}'.format(cmd))
    with open(bed_intersect_tmp, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile, check=True)
    with open(bed_intersect_tmp, 'r') as infile, open(bed_sort_tmp, 'w') as outfile:
        for line in infile:
            if line.startswith('#'):
                outfile.write(line)
            else:
                columns = line.rstrip('\n\r').split('\t')
                for i in range(annot_length, len(columns)):
                    if i != annot_length:
                        outfile.write('\t')
                    outfile.write(columns[i])
                outfile.write('\n')
    os.remove(bed_intersect_tmp)
    Bed.sort(bed_sort_tmp, bed_tag)
    os.remove(bed_sort_tmp)

    
if __name__ == '__main__':
    intersect()
