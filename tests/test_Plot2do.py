import logging
import math
import os
from pathlib import Path
from shutil import copyfile
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import Plot2do as p


@pytest.fixture
def mock_testclass():
    plot2do_sample = p.plot2do_sample
    run = subprocess.run
    yield
    p.plot2do_sample = plot2do_sample
    subprocess.run = run
    

def test_plot2do(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    samples_parent = samples.parent
    p.plot2do_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(p.plot2do, ['-f', samples])
    assert result.exit_code == 0
    p.plot2do_sample.assert_any_call(samples_parent / 'POLR2A', None, None, None, None, None, None, None, None, None, None, None, None, None)
    p.plot2do_sample.assert_any_call(samples_parent / 'ASDURF', None, None, None, None, None, None, None, None, None, None, None, None, None)
    p.plot2do_sample.assert_any_call(samples_parent / 'POLR1C', None, None, None, None, None, None, None, None, None, None, None, None, None)


def test_plot2do_second(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    samples_parent = samples.parent
    p.plot2do_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(p.plot2do, ['-f', samples, '-i', '1'])
    assert result.exit_code == 0
    p.plot2do_sample.assert_any_call(samples_parent / 'ASDURF', None, None, None, None, None, None, None, None, None, None, None, None, None)


def test_plot2do_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    samples_parent = samples.parent
    type = 'dyads'
    genome = 'mm9'
    reference = 'Plus1'
    sites = 'sites.txt'
    align = 'fivePrime'
    sitelabel = 'site-test'
    minlength = '100'
    maxlength = '300'
    upstream = '500'
    downstream = '600'
    colorscalemax = '0.05'
    simplifyplot = 'on'
    squeezeplot = 'on'
    p.plot2do_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(p.plot2do, ['-f', samples, '--type', type, '--genome', genome, '--reference', reference, '--sites', sites, '--align', align, '--siteLabel', sitelabel, '--minLength', minlength, '--maxLength', maxlength, '--upstream', upstream, '--downstream', downstream, '--colorScaleMax', colorscalemax, '--simplifyPlot', simplifyplot, '--squeezePlot', squeezeplot])
    assert result.exit_code == 0
    p.plot2do_sample.assert_any_call(samples_parent / 'POLR2A', type, genome, reference, sites, align, sitelabel, minlength, maxlength, upstream, downstream, colorscalemax, simplifyplot, squeezeplot)
    p.plot2do_sample.assert_any_call(samples_parent / 'ASDURF', type, genome, reference, sites, align, sitelabel, minlength, maxlength, upstream, downstream, colorscalemax, simplifyplot, squeezeplot)
    p.plot2do_sample.assert_any_call(samples_parent / 'POLR1C', type, genome, reference, sites, align, sitelabel, minlength, maxlength, upstream, downstream, colorscalemax, simplifyplot, squeezeplot)


def test_plot2do_filenotexists(testdir, mock_testclass):
    samples = 'samples.txt'
    p.plot2do_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(p.plot2do, ['-f', samples])
    assert result.exit_code != 0
    p.plot2do_sample.assert_not_called()


def test_plot2do_sample(testdir, mock_testclass):
    sample = 'POLR2A'
    bed = sample + '.bed'
    copyfile(Path(__file__).parent.joinpath('sample.bed'), bed)
    subprocess.run = MagicMock()
    p.plot2do_sample(sample)
    subprocess.run.assert_any_call(['Rscript', 'plot2DO.R', '-f', bed], check=True)


def test_plot2do_sample_parameters(testdir, mock_testclass):
    sample = 'POLR2A'
    bed = sample + '.bed'
    copyfile(Path(__file__).parent.joinpath('sample.bed'), bed)
    type = 'dyads'
    genome = 'mm9'
    reference = 'Plus1'
    sites = 'sites.txt'
    align = 'fivePrime'
    sitelabel = 'site-test'
    minlength = '100'
    maxlength = '300'
    upstream = '500'
    downstream = '600'
    colorscalemax = '0.05'
    simplifyplot = 'on'
    squeezeplot = 'on'
    subprocess.run = MagicMock()
    p.plot2do_sample(sample, type, genome, reference, sites, align, sitelabel, minlength, maxlength, upstream, downstream, colorscalemax, simplifyplot, squeezeplot)
    subprocess.run.assert_any_call(['Rscript', 'plot2DO.R', '--type', type, '--genome', genome, '--reference', reference, '--sites', sites, '--align', align, '--siteLabel', sitelabel, '--minLength', minlength, '--maxLength', maxlength, '--upstream', upstream, '--downstream', downstream, '--colorScaleMax', colorscalemax, '--simplifyPlot', simplifyplot, '--squeezePlot', squeezeplot, '-f', bed], check=True)


def test_plot2do_sample_bednotexists(testdir, mock_testclass):
    sample = 'POLR2A'
    bed = sample + '.bed'
    subprocess.run = MagicMock()
    p.plot2do_sample(sample)
    subprocess.run.assert_not_called()
