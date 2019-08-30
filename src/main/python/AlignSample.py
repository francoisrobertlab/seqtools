from distutils.command.check import check
import logging
from multiprocessing import Pool
import os
import re
import subprocess

import FullAnalysis
import click


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--fasta', '-f', type=click.Path(exists=True), default='sacCer3.fa',
              help='FASTA file used for alignment.')
@click.option('--threads', '-t', default=1,
              help='Number of threads used to process data per sample.')
@click.option('--poolThreads', '-T', default=2,
              help='Samples to process in parallel.')
def main(samples, fasta, threads, poolthreads):
    '''Align samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    bwa_index(fasta)
    samples_columns = FullAnalysis.all_columns(samples)
    samples_pool_args = []
    for sample_columns in samples_columns:
        sample = sample_columns[0]
        fastq = sample_columns[1] if len(sample_columns) > 1 else None
        samples_pool_args.append((sample, fastq, fasta, threads))
    with Pool(poolthreads) as pool:
        pool.starmap(align, samples_pool_args)


def align(sample, sample_fastq, fasta, threads=None):
    '''Align a single sample.'''
    print ('Running BWA on sample {}'.format(sample))
    if not sample_fastq:
        sample_fastq = sample
    fastq1 = fastq(sample_fastq, 1)
    fastq2 = fastq(sample_fastq, 2)
    bam_raw = sample + '-raw.bam'
    bwa(fastq1, fastq2, fasta, bam_raw, threads)


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


def bwa(fastq1, fastq2, fasta, bam_output, threads=None):
    '''Run BWA on FASTQ files.'''
    sam_output = bam_output + '.sam'
    cmd = ['bwa', 'mem']
    if not threads is None:
        cmd.extend(['-t', str(threads)])
    cmd.extend(['-o', sam_output, fasta, fastq1])
    if os.path.isfile(fastq2):
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


if __name__ == '__main__':
    main()
