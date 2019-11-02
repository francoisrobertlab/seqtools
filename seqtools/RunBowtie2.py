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
@click.option('--samples', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--threads', '-p', default=1,
              help='Number of threads used to process data per sample.')
@click.option('--index', type=int, default=None,
              help='Index of sample to process in samples file.')
@click.argument('bowtie_args', nargs=-1, type=click.UNPROCESSED)
def main(samples, threads, index, bowtie_args):
    '''Align samples using bowtie2 program.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        bowtie(sample, threads, bowtie_args)


def bowtie(sample, threads=None, bowtie_args=None):
    '''Align one sample using bowtie2 program.'''
    print ('Running bowtie2 on sample {}'.format(sample))
    fastq1 = Fastq.fastq(sample, 1)
    if fastq1 is None:
        raise AssertionError('Cannot find FASTQ files for sample ' + sample)
    fastq2 = Fastq.fastq(sample, 2)
    paired = fastq2 is not None and os.path.isfile(fastq2)
    bam_raw = sample + '-raw.bam'
    run_bowtie(fastq1, fastq2, bam_raw, threads, bowtie_args)


def run_bowtie(fastq1, fastq2, bam_output, threads=None, bowtie_args=None):
    '''Run bowtie2 on FASTQ files.'''
    sam_output = bam_output + '.sam'
    cmd = ['bowtie2'] + list(bowtie_args)
    if not threads is None:
        cmd.extend(['-p', str(threads)])
    cmd.extend(['-S', sam_output])
    if fastq2 is not None and os.path.isfile(fastq2):
        cmd.extend(['-1', fastq1, '-2', fastq2])
    else:
        cmd.extend(['-U', fastq1])
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
