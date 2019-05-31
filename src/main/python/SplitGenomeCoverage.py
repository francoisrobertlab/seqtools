import logging
import os

import GenomeCoverage
import click


@click.command()
@click.option('--samples', '-s', type=click.File('r'), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--sizes', type=click.Path(exists=True), default='sacCer3.chrom.sizes',
              help='Size of chromosomes.')
@click.option('--splitLength', type=int, default=25,
              help='Split reads in bins by their length.')
@click.option('--splitMinLength', default=100,
              help='Split reads minimum length.')
@click.option('--splitMaxLength', default=400,
              help='Split reads maximum length.')
def main(samples, sizes, splitlength, splitminlength, splitmaxlength):
    '''Split BED files from samples based on lenght of annotations.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    samples_lines = samples.read().splitlines()
    for sample_line in samples_lines:
        if sample_line.startswith('#'):
            continue
        sample_info = sample_line.rstrip("\n\r").split('\t');
        sample = sample_info[0]
        split_genome_coverage(sample, sizes, splitlength, splitminlength, splitmaxlength)


def split_genome_coverage(sample, sizes, splitlength, splitminlength, splitmaxlength):
    '''Split BED file from a single sample based on lenght of annotations.'''
    print ('Split BED file of sample {}'.format(sample))
    if splitlength is not None:
        bed_raw = sample + '-raw.bed'
        for bin_start in range(splitminlength, splitmaxlength, splitlength):
            bin_end = bin_start + splitlength
            bed_bin_raw = '{}-raw-{}-{}.bed'.format(sample, bin_start, bin_end)
            filter_annotations_by_length(bed_raw, bed_bin_raw, bin_start, bin_end)
            print ('Compute genome coverage on file {} for sample {}'.format(bed_bin_raw, sample))
            bed_bin_center = '{}-{}-{}-center.bed'.format(sample, bin_start, bin_end)
            GenomeCoverage.center_annotations(bed_bin_raw, bed_bin_center)
            bed = '{}-{}-{}.bed'.format(sample, bin_start, bin_end)
            bigwig = '{}-{}-{}.bw'.format(sample, bin_start, bin_end)
            count = GenomeCoverage.count_bed(bed_bin_center)
            if count == 0:
                GenomeCoverage.empty_bed(bed, sample)
                GenomeCoverage.bedgraph_to_bigwig(bed, bigwig, sizes)
            else :
                scale = GenomeCoverage.BASE_SCALE / count
                GenomeCoverage.coverage(bed_bin_center, bed, sizes, sample, scale)
                GenomeCoverage.bedgraph_to_bigwig(bed, bigwig, sizes)
            os.remove(bed_bin_center)


def filter_annotations_by_length(bed, output, minLength, maxLength):
    '''Filter BED file and keep only annotations that have specified size. Minimum size is included but max size is excluded.'''
    count = 0
    with open(bed, "r") as infile:
        with open(output, "w") as outfile:
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
