from distutils.command.check import check
import logging
import os
import re
import subprocess

import click
import pandas as pd


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--fasta', '-f', type=click.Path(exists=True), default='sacCer3.fa',
              help='FASTA file used for alignment.')
@click.option('--threads', '-t', default=1,
              help='Number of threads used to process data per sample.')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def main(samples, fasta, threads, index):
    '''Align samples using bwa program.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    if index == None:
        bwa_index(fasta)
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        bwa(sample, fasta, threads)


def bwa(sample, fasta, threads=None):
    '''Align one sample using bwa program.
    bwa output will be saved to a file called {sample}-raw.bam.
    FASTQ files are expected to follow the regular expression ``{sample}_R?1\.fastq(\.gz)?``.
    bwa output will be filtered to keep only properly paired reads and supplementary alignments will be removed.

    :param sample: sample name
    :param fasta: fasta file
    :param threads: number of threads that BWA will use - optional'''
    print ('Running BWA on sample {}'.format(sample))
    fastq1 = fastq(sample, 1)
    if fastq1 is None:
        raise AssertionError('Cannot find FASTQ files for sample ' + sample)
    fastq2 = fastq(sample, 2)
    paired = fastq2 is not None and os.path.isfile(fastq2)
    bam_raw = sample + '-raw.bam'
    run_bwa(fastq1, fastq2, fasta, bam_raw, threads)
    bam_filtered = sample + '-filtered.bam'
    filter_mapped(bam_raw, bam_filtered, paired, threads)
    bam_sort_input = bam_filtered
    if paired:
        bam_dedup = sample + '-dedup.bam'
        remove_duplicates(bam_filtered, bam_dedup, threads)
        bam_sort_input = bam_dedup
    bam = sample + '.bam'
    sort(bam_sort_input, bam, threads)


def fastq(sample, read=1):
    '''Returns existing FASTQ file for sample - read 1 or read 2, defaults to read 1.'''
    files = [f for f in os.listdir('.') if re.match('^' + re.escape(sample) + r'_R?' + str(read) + r'\.fastq(\.gz)?$', f)]
    if len(files) == 0:
        return None
    else:
        return files[0]


def bwa_index(fasta):
    '''Run BWA on FASTQ files.'''
    bwa_index_cmd = ['bwa', 'index', fasta]
    fasta_indexed = fasta + '.bwt';
    logging.debug('Running {}'.format(bwa_index_cmd))
    subprocess.run(bwa_index_cmd, check=True)
    if not os.path.isfile(fasta_indexed):
        raise AssertionError('Error when indexing FASTA ' + fasta)


def run_bwa(fastq1, fastq2, fasta, bam_output, threads=None):
    '''Run BWA on FASTQ files.'''
    sam_output = bam_output + '.sam'
    cmd = ['bwa', 'mem']
    if not threads is None:
        cmd.extend(['-t', str(threads)])
    cmd.extend(['-o', sam_output, fasta, fastq1])
    if fastq2 is not None and os.path.isfile(fastq2):
        cmd.append(fastq2)
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)
    if not os.path.isfile(sam_output):
        raise AssertionError('Error when running BWA with command ' + bwa_cmd)
    cmd = ['samtools', 'view', '-b']
    if not threads is None:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', bam_output, sam_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when converting SAM ' + sam_output + ' to BAM ' + bam_output)
    os.remove(sam_output)


def filter_mapped(bam_input, bam_output, paired, threads=None):
    '''Filter BAM file to remove poorly mapped sequences.'''
    cmd = ['samtools', 'view', '-b', '-F', '2048']
    if bool(paired):
        cmd.extend(['-f', '2'])
    else:
        cmd.extend(['-F', '4'])
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
