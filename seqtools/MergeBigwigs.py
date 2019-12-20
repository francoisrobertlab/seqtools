import argparse
from distutils.command.check import check
import logging
import os
import subprocess

import pyBigWig as pbw


def main():
    '''Merge bigWig files.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser(description='Merge bigWig files.')
    parser.add_argument('bigwigs', type=lambda x: is_valid_file(parser, x), nargs='+')
    parser.add_argument('--output', '-o', default='merge.bw')
    args = parser.parse_args()
    
    size = len(args.bigwigs)
    bws = [pbw.open(bw) for bw in args.bigwigs]
    chromosomes = {}
    for bw in bws:
        bw_chromosomes = bw.chroms()
        for chromosome in bw_chromosomes:
            if not chromosome in chromosomes or bw_chromosomes[chromosome] > chromosomes[chromosome]:
                chromosomes[chromosome] = bw_chromosomes[chromosome]
    chomosomes_sizes = 'chomosomes.size.tmp'
    with open('chomosomes.size.tmp', 'w') as output:
        for chromosome in chromosomes:
            output.write(chromosome)
            output.write('\t')
            output.write(str(chromosomes[chromosome]))
            output.write('\n')
    bed = 'bed.tmp'
    with open('bed.tmp', 'w') as output:
        output.write('track type=wiggle_0\n')
        for chromosome in chromosomes:
            size = chromosomes[chromosome]
            sums = [0] * size
            for bw in bws:
                values = bw.values(chromosome, 0, size)
                sums = [sums[i] + values[i] for i in range(0, size)]
            for i in range(0, len(sums)):
                output.write(chromosome)
                output.write('\t')
                output.write(str(i))
                output.write('\t')
                output.write(str(i + 1))
                output.write('\t')
                output.write(str(sums[i]))
                output.write('\n')
    bedgraph_to_bigwig(bed, args.output, chomosomes_sizes)
    os.remove(chomosomes_sizes)
    os.remove(bed)


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist" % arg)
    else:
        return arg


def bedgraph_to_bigwig(bed, bigwig, sizes):
    '''Converts bedgraph file to bigwig.'''
    cmd = ['bedGraphToBigWig', bed, sizes, bigwig]
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
