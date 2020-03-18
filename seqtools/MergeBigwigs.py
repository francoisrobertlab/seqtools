import logging
import math
import os

import click

import pyBigWig as pbw
from seqtools.bed import Bed
from seqtools.txt import Parser


@click.command()
@click.option('--merge', '-m', type=click.Path(exists=True), default='merge.txt', show_default=True,
              help='Merge name if first columns and sample names to merge on following columns - tab delimited.')
@click.option('--sizes', '-S', type=click.Path(exists=True), default='sacCer3.chrom.sizes', show_default=True,
              help='Size of chromosomes.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def mergebw(merge, sizes, index):
    '''Merge bigWig files related to samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    merge_samples(merge, sizes, index)


def merge_samples(merge='merge.txt', sizes='sacCer3.chrom.sizes', index=None):
    '''Merge bigWig files related to samples.'''
    merge_columns = Parser.columns(merge)
    if index != None:
        merge_columns = [merge_columns[index]]
    for columns in merge_columns:
        name = columns[0]
        samples = [sample for sample in columns[1:]]
        merge_sample(name, samples, sizes)

    
def merge_sample(name, samples, sizes):
    '''Merge bigWig files related to samples.'''
    print ('Merging samples {} into a single sample {}'.format(samples, name))
    sizes_columns = Parser.columns(sizes)
    bws = [pbw.open(sample + '.bw') for sample in samples]
    merged_bed_tmp = name + '-tmp.bed'
    with open(merged_bed_tmp, 'w') as output:
        output.write('track type=bedGraph name="' + name + '"\n')
        for size_columns in sizes_columns:
            chromosome = size_columns[0]
            size = size_columns[1]
            sums = [0] * size
            for bw in bws:
                bw_size = bw.chroms(chromosome) if bw.chroms(chromosome) else 0
                if bw_size == 0:
                    continue
                values = bw.values(chromosome, 0, min(size, bw_size))
                sums = [sums[i] + (values[i] if not math.isnan(values[i]) else 0) for i in range(0, min(size, bw_size))]
            for i in range(0, len(sums)):
                output.write(chromosome)
                output.write('\t')
                output.write(str(i))
                output.write('\t')
                output.write(str(i + 1))
                output.write('\t')
                output.write(str(sums[i]))
                output.write('\n')
    merged_bed_sort_tmp = name + '-tmp-sort.bed'
    Bed.sort(merged_bed_tmp, merged_bed_sort_tmp)
    merged_bw = name + '.bw'
    Bed.bedgraph_to_bigwig(merged_bed_sort_tmp, merged_bw, sizes)
    os.remove(merged_bed_sort_tmp)
    os.remove(merged_bed_tmp)


if __name__ == '__main__':
    mergebw()
