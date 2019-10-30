import logging
import os
import subprocess

import click
import pandas as pd
from pathlib import Path


@click.command()
@click.option('--file', '-f', type=click.Path(exists=True), default='samples.txt',
              help='Sample names listed one sample name by line.')
@click.option('--type', '-t', default=None,
              help='Type of distribution to plot.'
              ' [options: occ, dyads, fivePrime_ends, threePrime_ends; default = occ]')
@click.option('--genome', '-g', default=None,
              help='Genome version.'
              ' [options: sacCer3 (default) (S. cerevisiae); EF2 (S. pombe); dm3, dm6 (D. melanogaster);'
              ' ce10, ce11 (C. elegans); mm9, mm10 (M. musculus);'
              ' hg18, hg19, hg38 (H. sapiens); tair10 (A. thaliana)]')
@click.option('--reference', '-r', default=None,
              help='Reference points to be aligned.'
              ' [options: TSS (default), TTS, Plus1]')
@click.option('--sites', '-s', default=None,
              help='User-provided sites to be aligned (BED file).')
@click.option('--align', '-a', default=None,
              help='Points of the provided intervals to be aligned?'
              ' [options: center (default), fivePrime, threePrime]')
@click.option('--siteLabel', default=None,
              help='Label for the aligned sites.'
              ' [default = Sites]')
@click.option('--minLength', '-l', default=None,
              help='The smallest DNA fragment to be considered.'
              ' [default = 50]')
@click.option('--maxLength', '-L', default=None,
              help='The largest DNA fragment to be considered.'
              ' [default = 200]')
@click.option('--upstream', '-u', default=None,
              help='Length of the upstream region to be plotted.'
              ' [default = 1000]')
@click.option('--downstream', '-d', default=None,
              help='Length of the downstream region to be plotted.'
              ' [default = 1000]')
@click.option('--colorScaleMax', '-m', default=None,
              help='Maximum value on the color scale (e.g. 0.02).')
@click.option('--simplifyPlot', default=None,
              help='Simplify the plot (show only the 2D heat map).'
              ' [options: on, off (default)]')
@click.option('--squeezePlot', default=None,
              help='Simplify the plot and squeeze the heat map.'
              ' [options: on, off (default)]')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def main(file, type, genome, reference, sites, align, sitelabel, minlength, maxlength, upstream, downstream, colorscalemax, simplifyplot, squeezeplot, index):
    '''Run plot2DO on samples.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    file_parent = Path(file).parent
    sample_names = pd.read_csv(file, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        plot2do(file_parent / sample, type, genome, reference, sites, align, sitelabel, minlength, maxlength, upstream, downstream, colorscalemax, simplifyplot, squeezeplot)


def plot2do(sample, type, genome, reference, sites, align, sitelabel, minlength, maxlength, upstream, downstream, colorscalemax, simplifyplot, squeezeplot):
    '''Run plot2DO on a single sample.'''
    print ('Running plot2DO on sample {}'.format(sample))
    cmd = ['Rscript', 'plot2DO.R']
    if type:
        cmd.extend(['--type', type])
    if genome:
        cmd.extend(['--genome', genome])
    if reference:
        cmd.extend(['--reference', reference])
    if sites:
        cmd.extend(['--sites', sites])
    if align:
        cmd.extend(['--align', align])
    if sitelabel:
        cmd.extend(['--siteLabel', sitelabel])
    if minlength:
        cmd.extend(['--minLength', minlength])
    if maxlength:
        cmd.extend(['--maxLength', maxlength])
    if upstream:
        cmd.extend(['--upstream', upstream])
    if downstream:
        cmd.extend(['--downstream', downstream])
    if colorscalemax:
        cmd.extend(['--colorScaleMax', colorscalemax])
    if simplifyplot:
        cmd.extend(['--simplifyPlot', simplifyplot])
    if squeezeplot:
        cmd.extend(['--squeezePlot', squeezeplot])
    raw_bed = '{}-raw.bed'.format(sample)
    if os.path.isfile(raw_bed):
        raw_cmd = cmd
        raw_cmd.extend(['-f', raw_bed])
        logging.debug('Running {}'.format(raw_cmd))
        subprocess.call(raw_cmd)


if __name__ == '__main__':
    main()
