import logging
import os
import subprocess

import click

BASE_SCALE = 1000000


@click.command()
@click.option('--samples', '-s', type=click.File('r'), default='samples.txt',
              help='Sample names listed one sample name by line.'
              ' The first line is ignored.')
@click.option('--fasta', '-f', type=click.Path(exists=True), default='sacCer3.fa',
              help='FASTA file used for alignment.')
@click.option('--threads', '-t', default=1,
              help='Number of threads used to process data.')
def main(samples, fasta, threads):
    '''Align samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    bwa_index(fasta)
    next(samples)
    samples_lines = samples.read().splitlines()
    for sample_line in samples_lines:
        if sample_line.startswith('#'):
            continue
        sample_info = sample_line.split('\t');
        sample = sample_info[0]
        align(sample, fasta, threads)


def align(sample, fasta, threads):
    '''Align a single sample.'''
    print ('Analyse sample {}'.format(sample))
    fastq1 = fastq(sample, 1)
    fastq2 = fastq(sample, 2)
    bam_raw = sample + '-raw.bam'
    print ('Running BWA on sample {}'.format(sample))
    bwa(fastq1, fastq2, fasta, bam_raw, threads)
    print ('Filtering BAM and removing duplicates on sample {}'.format(sample))
    bam_filtered = sample + '-filtered.bam'
    filter_mapped(bam_raw, bam_filtered, threads)
    bam_dedup = sample + '-dedup.bam'
    remove_duplicates(bam_filtered, bam_dedup, threads)
    bam = sample + '.bam'
    sort(bam_dedup, bam, threads)


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
    subprocess.call(bwa_index_cmd)
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
    subprocess.call(cmd)
    if not os.path.isfile(sam_output):
        raise AssertionError('Error when running BWA with command ' + bwa_cmd)
    cmd = ['samtools', 'view', '-b']
    if not threads is None:
        cmd.extend(['--threads', str(threads - 1)])
    cmd.extend(['-o', bam_output, sam_output])
    logging.debug('Running {}'.format(cmd))
    subprocess.call(cmd)
    if not os.path.isfile(bam_output):
        raise AssertionError('Error when converting SAM ' + sam_output + ' to BAM ' + bam_output)
    os.remove(sam_output)


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
