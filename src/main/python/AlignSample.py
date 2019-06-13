from distutils.command.check import check
import logging
import os
import subprocess

import FullAnalysis
import click


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--fasta', '-f', type=click.Path(exists=True), default='sacCer3.fa',
              help='FASTA file used for alignment.')
@click.option('--threads', '-t', default=1,
              help='Number of threads used to process data.')
def main(samples, fasta, threads):
    '''Align samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    bwa_index(fasta)
    sample_names = FullAnalysis.first_column(samples)
    for sample in sample_names:
        align(sample, fasta, threads)


def align(sample, fasta, threads):
    '''Align a single sample.'''
    print ('Running BWA on sample {}'.format(sample))
    fastq1 = fastq(sample, 1)
    fastq2 = fastq(sample, 2)
    bam_raw = sample + '-raw.bam'
    bwa(fastq1, fastq2, fasta, bam_raw, threads)


def fastq(sample, read=1):
    '''Returns existing FASTQ file for sample - read 1 or read 2, defaults to read 1.'''
    fastq = sample + '_' + str(read) + '.fastq'
    if not os.path.isfile(fastq):
        fastq = sample + '_' + str(read) + '.fastq.gz'
    if not os.path.isfile(fastq):
        return None
    return fastq
    

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
