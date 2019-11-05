import datetime
import logging
import os
import re
import subprocess

import click
import pandas as pd


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
@click.option('--splitLength', type=int, default=10,
              help='Split reads in bins by their length.')
@click.option('--splitMinLength', default=100,
              help='Split reads minimum length.')
@click.option('--splitMaxLength', default=500,
              help='Split reads maximum length.')
def main(samples, index, splitlength, splitminlength, splitmaxlength):
    '''Split BED files from samples based on lenght of annotations.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        split_bed(sample, splitlength, splitminlength, splitmaxlength)


def split_bed(sample, splitlength, splitminlength, splitmaxlength):
    '''Split BED file from a single sample based on lenght of annotations.'''
    print ('Split BED file of sample {}'.format(sample))
    if splitlength is not None:
        bed_raw = sample + '-raw.bed'
        bed_sort = sample + '-raw-sort.bed'
        sort_bed_by_size(bed_raw, bed_sort)
        with open(bed_sort, 'r') as infile:
            line = infile.readline()
            length = annotation_length(line)
            for bin_start in range(splitminlength, splitmaxlength, splitlength):
                bin_end = bin_start + splitlength
                bin_file_tmp = '{}-{}-{}-tmp.bed'.format(sample, bin_start, bin_end)
                bin_file = '{}-{}-{}-raw.bed'.format(sample, bin_start, bin_end)
                with open(bin_file_tmp, 'w') as outfile:
                    while length < bin_end:
                        if length >= bin_start:
                            outfile.write(line)
                        line = infile.readline()
                        length = annotation_length(line)
                sort_bed(bin_file_tmp, bin_file)
                os.remove(bin_file_tmp)


def sort_bed_by_size(bed, output):
    '''Sort BED file by size'''
    cmd = ['bedtools', 'sort', '-sizeA', '-i', bed]
    logging.debug('Running {}'.format(cmd))
    with open(output, 'w') as outfile:
        subprocess.call(cmd, stdout=outfile)
    if not os.path.isfile(output):
        raise AssertionError('Error when sorting BED ' + bed)


def sort_bed(bed, output):
    '''Sort BED file by '''
    cmd = ['bedtools', 'sort', '-i', bed]
    logging.debug('Running {}'.format(cmd))
    with open(output, 'w') as outfile:
        subprocess.call(cmd, stdout=outfile)
    if not os.path.isfile(output):
        raise AssertionError('Error when sorting BED ' + bed)


def annotation_length(line):
    columns = line.rstrip('\r\n').split('\t')
    length = -1
    if len(columns) >= 3:
        length = int(columns[2]) - int(columns[1])
    return length


def splits(sample):
    '''Returns all splits for sample, sorted.'''
    regex = re.compile(sample + '-(\d+)-\d+-raw\.bed')
    files = os.listdir()
    beds = filter(regex.match, files)
    sample_splits = [bed[:-8] for bed in beds]
    sample_splits.sort(key=splitkey)
    return sample_splits


def splitkey(split):
    return int(re.search('\d+', split)[0])


if __name__ == '__main__':
    main()
