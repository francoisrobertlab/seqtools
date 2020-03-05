import logging

import click
from lmfit.models import GaussianModel, ConstantModel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seqtools.Split as sb


@click.command()
@click.option('--samples', '-s', type=click.Path(exists=True), default='samples.txt', show_default=True,
              help='Sample names listed one sample name by line.')
@click.option('--components', '-c', is_flag=True,
              help='Shows fit components and initial fit in plot.')
@click.option('--svg', is_flag=True,
              help='Save vectorial plot.')
@click.option('--verbose', '-v', is_flag=True,
              help='Shows fit report.')
@click.option('--center', type=float, default=None,
              help='Center of first gaussian. Defaults to 0')
@click.option('--cmin', '-cm', type=float, default=None,
              help='Minimum value for center of first gaussian. Defaults to unbounded')
@click.option('--cmax', '-cM', type=float, default=None,
              help='Maximum value for center of first gaussian. Defaults to unbounded')
@click.option('--amp', type=float, default=None,
              help='Amplitude of first gaussian. Defaults to maximum relative frequency, approximately')
@click.option('--amin', '-am', type=float, default=None,
              help='Minimum amplitude of first gaussian. Defaults to unbounded')
@click.option('--sigma', type=float, default=None,
              help='Width (sigma) of first gaussian. Defaults to half of maximum index')
@click.option('--smin', '-sm', type=float, default=None,
              help='Minimum width (sigma) of first gaussian. Defaults unbounded')
@click.option('--index', '-i', type=int, default=None,
              help='Index of sample to process in samples file.')
def fitgaussian(samples, components, svg, verbose, center, cmin, cmax, amp, amin, sigma, smin, index):
    '''Fits gaussian curve to dyad coverage.'''
    logging.basicConfig(filename='debug.log', level=logging.DEBUG, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sample_names = pd.read_csv(samples, header=None, sep='\t', comment='#')[0]
    if index != None:
        sample_names = [sample_names[index]]
    for sample in sample_names:
        fitgaussian_sample(sample, components, svg, verbose, center, cmin, cmax, amp, amin, sigma, smin)
        splits = sb.splits(sample)
        for split in splits:
            fitgaussian_sample(split, components, svg, verbose, center, cmin, cmax, amp, amin, sigma, smin)


def fitgaussian_sample(sample, components, svg, verbose, center, cmin, cmax, amp, amin, sigma, smin):
    '''Fits gaussian curve to dyad coverage for a single sample.'''
    print ('Fits gaussian curve to dyad coverage of sample {}'.format(sample))
    input = sample + '-dyad.txt'
    dyads = pd.read_csv(input, sep='\t', index_col=0, comment='#')
    x = dyads.index.values
    y = dyads['Relative Frequency'].values
    if not amp:
        amp = dyads['Relative Frequency'].max() * 100
    if not center:
        center = 0.0
    if not sigma:
        sigma = dyads.index.max() / 2
    plt.figure()
    plt.title(sample)
    plt.xlabel('Position relative to dyad (bp)')
    plt.ylabel('Relative Frequency')
    plt.xlim(x[0], x[len(x) - 1])
    plt.xticks(list(range(x[0], x[len(x) - 1] + 1, 25)))
    plt.plot(dyads.index.values, dyads['Relative Frequency'].values, color='red')
    plot_output = sample + '-dyad-gaussian.png'
    try:
        constant = ConstantModel(prefix='c_')
        pars = constant.make_params()
        pars['c_c'].set(value=dyads['Relative Frequency'].min(), min=0.0, max=dyads['Relative Frequency'].max())
        gauss = GaussianModel(prefix='g_')
        pars.update(gauss.make_params())
        pars['g_center'].set(value=center, min=cmin, max=cmax)
        pars['g_sigma'].set(value=sigma, min=smin)
        pars['g_amplitude'].set(value=amp, min=amin)
        mod = constant + gauss
        init = mod.eval(pars, x=x)
        out = mod.fit(y, pars, x=x)
        if components:
            plt.plot(x, init, 'b--', label='Initial fit')
        if verbose:
            print(out.fit_report(min_correl=0.5))
        plt.plot(x, out.best_fit, 'b-', label='Best fit')
        if components:
            comps = out.eval_components(x=x)
            plt.plot(x, np.repeat(comps['c_'], len(x)), 'g--', label='Constant component')
            plt.plot(x, comps['g_'], 'm--', label='Gaussian component')
    except Exception as e:
        logging.warning('could not fit gaussian curve to sample {}'.format(sample), e)
    if components:
        plt.legend(loc='lower right')
    plt.savefig(plot_output)
    if svg:
        plot_svg_output = sample + '-dyad-gaussian.svg'
        plt.savefig(plot_svg_output, transparent=True)
    plt.close()


if __name__ == '__main__':
    fitgaussian()
