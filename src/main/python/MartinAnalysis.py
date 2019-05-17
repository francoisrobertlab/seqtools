import logging
import os
import subprocess

import click


@click.command()
@click.option('--samples', '-s', type=click.File('r'), default='samples.txt',
              help='Sample names listed one sample name by line.'
              ' The first line is ignored.'
              ' An SRR id can be provided (tab-separated) to download the FASTQ file automatically, otherwise it must be provided.')
@click.option('--merge', '-m', type=click.File('r'), default='merge.txt',
              help='Merge samples listed on the same line (tab-separated).'
              ' The first line is ignored.')
@click.option('--fasta', '-f', type=click.Path(exists=True), default='sacCer3.fa',
              help='FASTA file used for alignment.')
@click.option('--sizes', type=click.Path(exists=True), default='sacCer3.chrom.sizes',
              help='Size of chromosomes.')
@click.option('--threads', '-t', default=1,
              help='Number of threads used to process data.')
def main(samples, merge, fasta, sizes, threads):
    '''Analyse Martin et al. data from November 2018 in Genetics'''
    # logging.basicConfig(level=logging.WARN)
    logging.basicConfig(filename='debug.log', level=logging.DEBUG)
    next(samples)
    samples_lines = samples.read().splitlines()
    for sample_line in samples_lines:
        if sample_line.startswith('#'):
            continue
        sample_info = sample_line.split('\t');
        sample = sample_info[0]
        srr = sample_info[1] if len(sample_info) > 1 else None
        analyse(sample, srr, fasta)


def analyse(sample, srr, fasta):
    '''Analyse sample'''
    print ('Analyse sample {}'.format(sample))
    fastq = sample + '_1.fastq'
    fastq2 = sample + '_2.fastq'
    if not os.path.isfile(fastq):
        fastq = sample + '_1.fastq.gz'
        fastq2 = sample + '_2.fastq.gz'
    if not os.path.isfile(fastq):
        fastq = sample + '.fastq'
    if not os.path.isfile(fastq):
        fastq = sample + '.fastq.gz'
    if not os.path.isfile(fastq) and srr:
        fastq = sample + '_1.fastq.gz'
        fastq2 = sample + '_2.fastq.gz'
        print ('Downloading FASTQ for {} with SRR {}'.format(sample, srr))
        download(sample, srr)
    bam_raw = sample + '-raw.bam'
    print ('Running BWA on sample {}'.format(sample))
    bwa(fastq, fastq2, fasta, bam_raw)


def download(sample, srr):
    '''Download reads of a sample'''
    fastq = sample + '_1.fastq.gz'
    fastq2 = sample + '_2.fastq.gz'
    srr_output = srr + '_1.fastq.gz'
    srr_output2 = srr + '_2.fastq.gz'
    fastqdump_cmd = ['fastq-dump', '--split-files', '--gzip', srr];
    logging.debug('Running {}'.format(fastqdump_cmd))
    subprocess.call(fastqdump_cmd)
    os.rename(srr_output, fastq)
    if os.path.isfile(srr_output2):
        os.rename(srr_output2, fastq2)


def bwa(fastq, fastq2, fasta, bam_output):
    '''Run BWA on FASTQ files'''
    bwa_index_cmd = ['bwa', 'index', fasta]
    fasta_indexed = fasta + '.amb';
    if not os.path.isfile(fasta_indexed):
        logging.debug('Running {}'.format(bwa_index_cmd))
        subprocess.call(bwa_index_cmd)
    sam_output = bam_output + '.sam'
    bwa_cmd = ['bwa', 'mem', '-o', sam_output, fasta, fastq]
    if os.path.isfile(fastq2):
        bwa_cmd.append(fastq2)
    logging.debug('Running {}'.format(bwa_cmd))
    subprocess.call(bwa_cmd)
    sam_to_bam_cmd = ['samtools', 'view', '-b', '-o', bam_output, sam_output]
    logging.debug('Running {}'.format(sam_to_bam_cmd))
    subprocess.call(sam_to_bam_cmd)
    os.remove(sam_output)


if __name__ == '__main__':
    main()
