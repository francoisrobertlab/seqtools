import logging
import os
import subprocess

import AlignSample
import click


@click.command()
@click.option('--samples', '-s', type=click.File('r'), default='samples.txt',
              help='Sample names listed one sample name by line.'
              ' An SRR id can be provided (tab-separated) to download the FASTQ file automatically, otherwise it must be provided.')
def main(samples):
    '''Download reads of all samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    samples_lines = samples.read().splitlines()
    for sample_line in samples_lines:
        if sample_line.startswith('#'):
            continue
        sample_info = sample_line.split('\t');
        sample = sample_info[0]
        srr = sample_info[1] if len(sample_info) > 1 else None
        download(sample, srr)


def download(sample, srr):
    '''Download reads of a sample.'''
    fastq = AlignSample.fastq(sample, 1)
    if fastq is None and srr:
        print ('Downloading FASTQ for sample {} with SRR {}'.format(sample, srr))
        fastq = sample + '_1.fastq.gz'
        fastq2 = sample + '_2.fastq.gz'
        srr_output = srr + '_1.fastq.gz'
        srr_output2 = srr + '_2.fastq.gz'
        cmd = ['fastq-dump', '--split-files', '--gzip', srr];
        logging.debug('Running {}'.format(cmd))
        subprocess.call(cmd)
        os.rename(srr_output, fastq)
        if os.path.isfile(srr_output2):
            os.rename(srr_output2, fastq2)


if __name__ == '__main__':
    main()
