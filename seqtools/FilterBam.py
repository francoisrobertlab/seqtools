from distutils.command.check import check
import logging
import os
import subprocess

import click
import pandas as pd


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--paired/--unpaired', '-p/-u', default=True, show_default=True,
              help='Sample reads are paired')
@click.option('--threads', '-t', default=1, show_default=True,
              help='Number of threads used to process data per sample.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def main(samples, paired, threads, index):
    '''Filter BAM file to keep only properly paired reads and remove supplementary alignments and duplicates.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        filter_bam(sample, paired, threads)


def filter_bam(sample, paired, threads=None):
    '''Filter BAM file to keep only properly paired reads and remove supplementary alignments and duplicates.'''
    print ('Filtering BAM for sample {}'.format(sample))
    bam_raw = sample + '-raw.bam'
    bam_filtered = sample + '-filtered.bam'
    filter_mapped(bam_raw, bam_filtered, paired, threads)
    bam_sort_input = bam_filtered
    if paired:
        bam_dedup = sample + '-dedup.bam'
        remove_duplicates(bam_filtered, bam_dedup, threads)
        bam_sort_input = bam_dedup
    bam = sample + '.bam'
    sort(bam_sort_input, bam, threads)


def filter_mapped(bam_input, bam_output, paired, threads=None):
    '''Filter BAM file to remove poorly mapped sequences.'''
    print ('Filtering BAM {} to remove poorly mapped sequences'.format(bam_input))
    cmd = ['samtools', 'view', '-b', '-F', '2048']
    if bool(paired):
        cmd.extend(['-f', '2'])
    else:
        cmd.extend(['-F', '4'])
    if not threads is None and threads > 1:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', bam_output, bam_input])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)


def remove_duplicates(bam_input, bam_output, threads=None):
    '''Remove duplicated sequences from BAM file.'''
    print ('Removing duplicated sequences from BAM {}'.format(bam_input))
    fixmate_output = bam_input + '.fix'
    cmd = ['samtools', 'fixmate', '-m']
    if not threads is None and threads > 1:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend([bam_input, fixmate_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)
    sort_output = bam_input + '.sort'
    cmd = ['samtools', 'sort']
    if not threads is None and threads > 1:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', sort_output, fixmate_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)
    os.remove(fixmate_output)
    cmd = ['samtools', 'markdup', '-r']
    if not threads is None and threads > 1:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend([sort_output, bam_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)
    os.remove(sort_output)


def sort(bam_input, bam_output, threads=None):
    '''Sort BAM file.'''
    print ('Sorting BAM {}'.format(bam_input))
    cmd = ['samtools', 'sort']
    if not threads is None and threads > 1:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', bam_output, bam_input])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)


if __name__ == '__main__':
    main()
