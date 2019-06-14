import logging
import FullAnalysis
import click


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--annotations', '-a', type=click.Path(exists=True), default='annotations.bed',
              help='Keep reads for which their center is located on specified annotations.')
def main(samples, annotations):
    '''Filter raw BED files to keep reads for which their center overlap specified annotations.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    samples_names = FullAnalysis.first_column(samples)
    parsed_annotations = FullAnalysis.all_columns(annotations)
    annotations_chrom = {}
    for annotation in parsed_annotations:
        chromosome = annotation[0]
        strand = annotation[5]
        if not chromosome in annotations_chrom:
            annotations_chrom[chromosome] = {}
        if not strand in annotations_chrom[chromosome]:
            annotations_chrom[chromosome][strand] = []
        annotations_chrom[chromosome][strand].append(annotation)
    for sample in samples_names:
        filter_bed(sample, annotations_chrom)


def filter_bed(sample, annotations_chrom):
    '''Filter raw BED files to keep reads for which their center overlap specified annotations.'''
    print ('Filter raw BED files with annotations for sample {}'.format(sample))
    bed_raw = sample + '-raw.bed'
    bed_filtered = sample + '-annotations.bed'
    with open(bed_raw) as infile:
        with open(bed_filtered, 'w') as outfile:
            for line in infile:
                if line.startswith('track') or line.startswith('browser') or line.startswith('#'):
                    continue
                read = line.rstrip('\n\r').split('\t')
                chromosome = read[0]
                strand = read[5]
                if chromosome in annotations_chrom and strand in annotations_chrom[chromosome]:
                    if overlap_any_annotation(read, annotations_chrom[read[0]][read[5]]):
                        outfile.write('\t'.join(read))
                        outfile.write('\n')


def overlap_any_annotation(read, annotations):
    '''Returns true if read overlaps any annotations, otherwise, returns false.'''
    center = (int(read[2]) - int(read[1])) / 2 + int(read[1])
    ret = False
    for annotation in annotations:
        ret |= read[0] == annotation[0] and read[5] == annotation[5] and center >= int(annotation[1]) and center < int(annotation[2])
    return ret


if __name__ == '__main__':
    main()
