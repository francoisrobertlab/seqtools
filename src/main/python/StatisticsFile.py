import logging
import os
import re
import subprocess
import GenomeCoverage

import SplitGenomeCoverage
import click


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--merge', '-m', type=click.Path(), default='merge.txt',
              help='Merge name if first columns and sample names to merge on following columns - tab delimited.')
@click.option('--output', '-o', type=click.Path(), default='statistics.txt',
              help='Output file were statistics are written.')
def main(samples, merge, output):
    '''Creates statistics file for samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = first_column(samples)
    merge_names = first_column(merge) if os.path.isfile(merge) else None
    all_statistics(sample_names, merge_names, output)


def first_column(file):
    fcs = []
    with open(file, 'r') as lines:
        for line in lines:
            if line.startswith('#'):
                continue
            fc = line.rstrip("\n\r").split('\t')[0];
            fcs.extend([fc])
    return fcs


def all_statistics(samples, merges, output):
    stats_headers = None
    samples_stats = []
    for sample in samples:
        if stats_headers is None:
            stats_headers = headers(sample)
        sample_stats = sample_statistics(sample)
        samples_stats.extend([sample_stats])
    if merges:
        for merge in merges:
            sample_stats = sample_statistics(merge)
            samples_stats.extend([sample_stats])
    with open(output, 'w') as out:
        out.write('\t'.join(stats_headers))
        out.write('\n')
        for sample_stats in samples_stats:
            out.write('\t'.join([str(value) for value in sample_stats]))
            out.write('\n')


def headers(sample):
    '''Statistics headers'''
    headers = ['Sample', 'Total reads', 'Mapped reads', 'Deduplicated reads']
    splits = SplitGenomeCoverage.splits(sample)
    splits_headers = [split[len(sample) + 1:] for split in splits]
    headers.extend(splits_headers)
    return headers
    
    
def sample_statistics(sample):
    '''Statistics of a single sample.'''
    print ('Computing statistics for sample {}'.format(sample))
    sample_stats = [sample]
    bam_raw = sample + '-raw.bam'
    sample_stats.extend([flagstat_total(bam_raw) if os.path.isfile(bam_raw) else ''])
    bam_filtered = sample + '-filtered.bam'
    sample_stats.extend([flagstat_total(bam_filtered) if os.path.isfile(bam_filtered) else ''])
    bed_raw = sample + '-raw.bed'
    sample_stats.extend([GenomeCoverage.count_bed(bed_raw) * 2 if os.path.isfile(bed_raw) else ''])
    splits = SplitGenomeCoverage.splits(sample)
    if splits:
        beds = [split + '-raw.bed' for split in splits]
        counts = [GenomeCoverage.count_bed(bed) for bed in beds]
        sample_stats.extend(counts)
    return sample_stats


def flagstat_total(bam):
    cmd = ['samtools', 'flagstat', bam]
    logging.debug('Running {}'.format(cmd))
    output = subprocess.run(cmd, capture_output=True, check=True)
    return re.search("^\d+", output.stdout.decode("utf-8")).group()


if __name__ == '__main__':
    main()
