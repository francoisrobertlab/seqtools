from distutils.command.check import check
import logging
import os
import subprocess

import click
import pandas as pd
from seqtools.seq import Fastq


@click.command(context_settings=dict(
    ignore_unknown_options=True,
))
@click.option('--samples', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--fasta', type=click.Path(exists=True), default='sacCer3.fa', show_default=True,
              help='FASTA file used for alignment.')
@click.option('--threads', '-t', default=1, show_default=True,
              help='Number of threads used to process data per sample.')
@click.option('--index', type=int, default=None,
              help='Index of sample to process in samples file.')
@click.argument('bwa_args', nargs=-1, type=click.UNPROCESSED)
def main(samples, fasta, threads, index, bwa_args):
    '''Align samples using bwa program.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    if index == None:
        bwa_index(fasta)
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        bwa(sample, fasta, threads, bwa_args)


def bwa(sample, fasta, threads=None, bwa_args=None):
    '''Align one sample using bwa program.'''
    print ('Running BWA on sample {}'.format(sample))
    fastq1 = Fastq.fastq(sample, 1)
    if fastq1 is None:
        raise AssertionError('Cannot find FASTQ files for sample ' + sample)
    fastq2 = Fastq.fastq(sample, 2)
    paired = fastq2 is not None and os.path.isfile(fastq2)
    bam_raw = sample + '-raw.bam'
    run_bwa(fastq1, fastq2, fasta, bam_raw, threads, bwa_args)


def bwa_index(fasta):
    '''Run BWA on FASTQ files.'''
    bwa_index_cmd = ['bwa', 'index', fasta]
    fasta_indexed = fasta + '.bwt';
    logging.debug('Running {}'.format(bwa_index_cmd))
    subprocess.run(bwa_index_cmd, check=True)


def run_bwa(fastq1, fastq2, fasta, bam_output, threads=None, bwa_args=None):
    '''Run BWA on FASTQ files.'''
    sam_output = bam_output + '.sam'
    cmd = ['bwa', 'mem'] + list(bwa_args)
    if not threads is None:
        cmd.extend(['-t', str(threads)])
    cmd.extend(['-o', sam_output, fasta, fastq1])
    if fastq2 is not None and os.path.isfile(fastq2):
        cmd.append(fastq2)
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)
    cmd = ['samtools', 'view', '-b']
    if not threads is None:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', bam_output, sam_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.run(cmd, check=True)
    os.remove(sam_output)


if __name__ == '__main__':
    main()
