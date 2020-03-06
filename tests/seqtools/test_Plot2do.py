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
    plot2do_samples = p.plot2do_samples
    plot2do_sample = p.plot2do_sample
    run = subprocess.run
    yield
    p.plot2do_samples = plot2do_samples
    p.plot2do_sample = plot2do_sample
    subprocess.run = run
    

def test_plot2do(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    p.plot2do_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(p.plot2do, ['-f', samples])
    assert result.exit_code == 0
    p.plot2do_samples.assert_called_once_with(samples, None, ())


def test_plot2do_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    type = 'dyads'
    genome = 'mm9'
    index = 1
    p.plot2do_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(p.plot2do, ['-f', samples, '--type', type, '--genome', genome, '--index', index])
    assert result.exit_code == 0
    p.plot2do_samples.assert_called_once_with(samples, index, ('--type', type, '--genome', genome,))


def test_plot2do_filenotexists(testdir, mock_testclass):
    samples = 'samples.txt'
    p.plot2do_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(p.plot2do, ['-f', samples])
    assert result.exit_code != 0
    p.plot2do_samples.assert_not_called()


def test_plot2do_samples(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    samples_parent = samples.parent
    p.plot2do_sample = MagicMock()
    p.plot2do_samples(samples)
    p.plot2do_sample.assert_any_call(samples_parent / 'POLR2A', ())
    p.plot2do_sample.assert_any_call(samples_parent / 'ASDURF', ())
    p.plot2do_sample.assert_any_call(samples_parent / 'POLR1C', ())


def test_plot2do_samples_second(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    samples_parent = samples.parent
    p.plot2do_sample = MagicMock()
    p.plot2do_samples(samples, index=1)
    p.plot2do_sample.assert_any_call(samples_parent / 'ASDURF', ())


def test_plot2do_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    samples_parent = samples.parent
    plot2do_args = ('--type', 'dyads', '--genome', 'mm9',)
    p.plot2do_sample = MagicMock()
    p.plot2do_samples(samples, plot2do_args=plot2do_args)
    p.plot2do_sample.assert_any_call(samples_parent / 'POLR2A', plot2do_args)
    p.plot2do_sample.assert_any_call(samples_parent / 'ASDURF', plot2do_args)
    p.plot2do_sample.assert_any_call(samples_parent / 'POLR1C', plot2do_args)


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
    plot2do_args = ('--type', 'dyads', '--genome', 'mm9',)
    subprocess.run = MagicMock()
    p.plot2do_sample(sample, plot2do_args)
    subprocess.run.assert_any_call(['Rscript', 'plot2DO.R', '--type', 'dyads', '--genome', 'mm9', '-f', bed], check=True)


def test_plot2do_sample_bednotexists(testdir, mock_testclass):
    sample = 'POLR2A'
    bed = sample + '.bed'
    subprocess.run = MagicMock()
    p.plot2do_sample(sample)
    subprocess.run.assert_not_called()
