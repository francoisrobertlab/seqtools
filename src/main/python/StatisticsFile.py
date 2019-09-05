import logging
import os
import re
import subprocess

import FilterBed
import FullAnalysis
import GenomeCoverage
import SplitBed
import click


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--merge', '-m', type=click.Path(), default='merge.txt',
              help='Merge name if first columns and sample names to merge on following columns - tab delimited.')
@click.option('--annotations', '-a', type=click.Path(), default='annotations.bed',
              help='File used for FilterBed script, if any.')
@click.option('--output', '-o', type=click.Path(), default='statistics.txt',
              help='Output file were statistics are written.')
def main(samples, merge, annotations, output):
    '''Creates statistics file for samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = FullAnalysis.first_column(samples)
    merge_names = FullAnalysis.first_column(merge) if os.path.isfile(merge) else None
    all_statistics(sample_names, merge_names, annotations, output)


def all_statistics(samples, merges, annotations, output):
    stats_headers = None
    samples_stats = []
    for sample in samples:
        if stats_headers is None:
            stats_headers = headers(sample, annotations)
        sample_stats = sample_statistics(sample, annotations)
        samples_stats.append(sample_stats)
    if merges:
        for merge in merges:
            sample_stats = sample_statistics(merge, annotations)
            samples_stats.append(sample_stats)
    with open(output, 'w') as out:
        out.write('\t'.join(stats_headers))
        out.write('\n')
        for sample_stats in samples_stats:
            out.write('\t'.join([str(value) for value in sample_stats]))
            out.write('\n')


def headers(sample, annotations):
    '''Statistics headers'''
    headers = ['Sample', 'Total reads', 'Mapped reads', 'Deduplicated reads']
    splits = SplitBed.splits(sample)
    splits_headers = [split[len(sample) + 1:] for split in splits]
    headers.extend(splits_headers)
    filtereds = FilterBed.filtered(sample, annotations)
    filtereds_headers = [filtered[len(sample) + 1:] for filtered in filtereds]
    headers.extend(filtereds_headers)
    return headers
    
    
def sample_statistics(sample, annotations):
    '''Statistics of a single sample.'''
    print ('Computing statistics for sample {}'.format(sample))
    sample_stats = [sample]
    bam_raw = sample + '-raw.bam'
    sample_stats.extend([flagstat_total(bam_raw) if os.path.isfile(bam_raw) else ''])
    bam_filtered = sample + '-filtered.bam'
    sample_stats.extend([flagstat_total(bam_filtered) if os.path.isfile(bam_filtered) else ''])
    bed_raw = sample + '-raw.bed'
    sample_stats.extend([GenomeCoverage.count_bed(bed_raw) * 2 if os.path.isfile(bed_raw) else ''])
    splits = SplitBed.splits(sample)
    if splits:
        beds = [split + '-raw.bed' for split in splits]
        counts = [GenomeCoverage.count_bed(bed) for bed in beds]
        sample_stats.extend(counts)
    filtereds = FilterBed.filtered(sample, annotations)
    if filtereds:
        beds = [filtered + '.bed' for filtered in filtereds]
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
