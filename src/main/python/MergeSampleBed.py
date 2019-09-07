import logging
from multiprocessing import Pool
import os
import subprocess

import FullAnalysis
import click


@click.command()
@click.option('--merge', '-s', type=click.File('r'), default='merge.txt',
              help='Merge name if first columns and sample names to merge on following columns - tab delimited.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def main(merge, index):
    '''Merge BED files related to samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    merges_columns = FullAnalysis.all_columns(merge)
    if index != None:
        merges_columns = [merges_columns[index]]
    for merge_columns in merges_columns:
        name = merge_info[0]
        samples = merge_info[1:]
        merge_samples(name, samples)


def merge_samples(name, samples):
    '''Merge BED files related to samples.'''
    print ('Merging samples {} into a single sample {}'.format(samples, name))
    merged_bed_tmp = name + '-tmp.bed'
    with open(merged_bed_tmp, "w") as outfile:
        for sample in samples:
            sample_bed = sample + '-raw.bed'
            with open(sample_bed, "r") as infile:
                for line in infile:
                    outfile.write(line)
    merged_bed = name + '-raw.bed'
    cmd = ['bedtools', 'sort', '-i', merged_bed_tmp]
    logging.debug('Running {}'.format(cmd))
    with open(merged_bed, "w") as outfile:
        subprocess.call(cmd, stdout=outfile)
    if not os.path.isfile(merged_bed):
        raise AssertionError('Error when sorting BED ' + merged_bed_tmp)
    os.remove(merged_bed_tmp)


if __name__ == '__main__':
    main()
