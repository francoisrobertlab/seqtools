import logging
import os
import subprocess

import click


@click.command()
@click.option('--samples', '-s', type=click.File('r'), default='samples.txt',
              help='Sample names listed one sample name by line.'
              ' The first line is ignored.')
@click.option('--threads', '-t', default=1,
              help='Number of threads used to process data.')
def main(samples, threads):
    '''Filter BAM files to remove unmapped sequences and duplicates.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    next(samples)
    samples_lines = samples.read().splitlines()
    for sample_line in samples_lines:
        if sample_line.startswith('#'):
            continue
        sample_info = sample_line.split('\t');
        sample = sample_info[0]
        filter_bam(sample, threads)


def filter_bam(sample, threads):
    '''Filter BAM files to remove unmapped sequences and duplicates on a single sample.'''
    print ('Filtering BAM and removing duplicates on sample {}'.format(sample))
    bam_raw = sample + '-raw.bam'
    bam_filtered = sample + '-filtered.bam'
    filter_mapped(bam_raw, bam_filtered, threads)
    bam_dedup = sample + '-dedup.bam'
    remove_duplicates(bam_filtered, bam_dedup, threads)
    bam = sample + '.bam'
    sort(bam_dedup, bam, threads)


def filter_mapped(bam_input, bam_output, threads=None):
    '''Filter BAM file to remove poorly mapped sequences.'''
    cmd = ['samtools', 'view', '-f', '2', '-F', '2048', '-b']
    if not threads is None:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', bam_output, bam_input])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when filtering BAM ' + bam_input)


def remove_duplicates(bam_input, bam_output, threads=None):
    '''Remove duplicated sequences from BAM file.'''
    fixmate_output = bam_input + '.fix'
    cmd = ['samtools', 'fixmate', '-m']
    if not threads is None:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend([bam_input, fixmate_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(fixmate_output):
        raise AssertionError('Error when fixing duplicates in BAM ' + bam_input)
    sort_output = bam_input + '.sort'
    cmd = ['samtools', 'sort']
    if not threads is None:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', sort_output, fixmate_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(sort_output):
        raise AssertionError('Error when sorting BAM ' + fixmate_output)
    os.remove(fixmate_output)
    cmd = ['samtools', 'markdup', '-r']
    if not threads is None:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend([sort_output, bam_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when removing duplicates from BAM ' + sort_output)
    os.remove(sort_output)


def sort(bam_input, bam_output, threads=None):
    '''Sort BAM file.'''
    cmd = ['samtools', 'sort']
    if not threads is None:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', bam_output, bam_input])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when sorting BAM ' + bam_input)


if __name__ == '__main__':
    main()
