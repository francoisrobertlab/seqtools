import logging
import subprocess

import click
from seqtools.txt import Parser


@click.command()
@click.option('--merge', '-m', type=click.Path(exists=True), default='merge.txt', show_default=True,
              help='Merge name if first columns and sample names to merge on following columns - tab delimited.')
@click.option('--suffix', '-s', default='', show_default=True,
              help='Suffix added to merge name and sample name in BAM filename.')
@click.option('--threads', '-t', default=1, show_default=True,
              help='Number of threads used to process data per sample.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def mergebam(merge, suffix, threads, index):
    '''Merge BAM files related to samples.'''
    logging.basicConfig(filename='seqtools.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    merge_samples(merge, suffix, threads, index)


def merge_samples(merge='merge.txt', suffix='', threads=1, index=None):
    '''Merge BAM files related to samples.'''
    merge_columns = Parser.columns(merge)
    if index != None:
        merge_columns = [merge_columns[index]]
    for columns in merge_columns:
        name = columns[0]
        samples = [sample for sample in columns[1:]]
        merge_sample(name, samples, suffix, threads)


def merge_sample(name, samples, suffix='', threads=1):
    '''Merge BAM files related to samples.'''
    print ('Merging samples {} into a single sample {}'.format(samples, name))
    for sample in samples:
        bam_input = sample + suffix + '.bam'
        cmd = ['samtools', 'index']
        if not threads is None and threads > 1:
            cmd.extend(['-@', str(threads - 1)])
        cmd.extend([bam_input])
        logging.debug('Running {}'.format(cmd))
        subprocess.run(cmd, check=True)
    bams_input = [sample + suffix + '.bam' for sample in samples]
    bam_output = name + suffix + '.bam'
    cmd = ['samtools', 'merge']
    if not threads is None and threads > 1:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend([bam_output])
    cmd.extend(bams_input)
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    merge()
