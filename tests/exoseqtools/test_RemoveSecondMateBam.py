import logging
import os
from pathlib import Path
import pytest
from shutil import copyfile
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
from exoseqtools import RemoveSecondMateBam as rsm
from seqtools.txt import Parser


@pytest.fixture
def mock_testclass():
    removesecondmate_samples = rsm.removesecondmate_samples
    removesecondmate_sample = rsm.removesecondmate_sample
    first = Parser.first
    yield
    rsm.removesecondmate_samples = removesecondmate_samples
    rsm.removesecondmate_sample = removesecondmate_sample
    Parser.first = first
    
    
def test_removesecondmate(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    rsm.removesecondmate_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(rsm.removesecondmate, ['-s', samples])
    assert result.exit_code == 0
    rsm.removesecondmate_samples.assert_called_once_with(samples, 1, None)


def test_removesecondmate_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    rsm.removesecondmate_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(rsm.removesecondmate, ['-s', samples, '-t', threads])
    assert result.exit_code == 0
    rsm.removesecondmate_samples.assert_called_once_with(samples, threads, None)


def test_removesecondmate_second(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    index = 1
    rsm.removesecondmate_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(rsm.removesecondmate, ['-s', samples, '-i', index])
    assert result.exit_code == 0
    rsm.removesecondmate_samples.assert_called_once_with(samples, 1, index)


def test_removesecondmate_samplesnotexists(testdir, mock_testclass):
    samples = 'samples.txt'
    rsm.removesecondmate_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(rsm.removesecondmate, ['-s', samples])
    assert result.exit_code != 0
    rsm.removesecondmate_samples.assert_not_called()


def test_removesecondmate_samples(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    Parser.first = MagicMock(return_value=samples)
    rsm.removesecondmate_sample = MagicMock()
    rsm.removesecondmate_samples(samples_file)
    Parser.first.assert_called_once_with(samples_file)
    for sample in samples:
        rsm.removesecondmate_sample.assert_any_call(sample, None)


def test_removesecondmate_samples_second(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    Parser.first = MagicMock(return_value=samples)
    rsm.removesecondmate_sample = MagicMock()
    rsm.removesecondmate_samples(samples_file, index=1)
    Parser.first.assert_called_once_with(samples_file)
    rsm.removesecondmate_sample.assert_any_call(samples[1], None)

    
def test_removesecondmate_sample(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    output = sample + '-mate1.bam'
    subprocess.run = MagicMock()
    rsm.removesecondmate_sample(sample)
    subprocess.run.assert_any_call(['samtools', 'view', '-f', '64', '-b', '-o', output, bam], check=True)


def test_removesecondmate_sample_threads(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    output = sample + '-mate1.bam'
    subprocess.run = MagicMock()
    rsm.removesecondmate_sample(sample, 3)
    subprocess.run.assert_any_call(['samtools', 'view', '--threads', '2', '-f', '64', '-b', '-o', output, bam], check=True)

