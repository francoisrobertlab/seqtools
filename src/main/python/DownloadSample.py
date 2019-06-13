import logging
import os
import subprocess

import AlignSample
import FullAnalysis
import click


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.'
              ' An SRR id can be provided (tab-separated) to download the FASTQ file automatically, otherwise FASTQ file must be provided.')
def main(samples):
    '''Download reads of all samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    samples_columns = FullAnalysis.all_columns(samples)
    for sample_columns in samples_columns:
        sample = sample_columns[0]
        srr = sample_columns[1] if len(sample_columns) > 1 else None
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
