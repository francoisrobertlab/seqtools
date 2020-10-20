import logging
import os
from pathlib import Path
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import MergeBam as mb


@pytest.fixture
def mock_testclass():
    merge_samples = mb.merge_samples
    merge_sample = mb.merge_sample
    run = subprocess.run
    yield
    mb.merge_samples = merge_samples
    mb.merge_sample = merge_sample
    subprocess.run = run
    

def create_file(*args, **kwargs):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_mergebam(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.mergebam, ['-m', merge])
    assert result.exit_code == 0
    mb.merge_samples.assert_called_once_with(merge, '', 1, None)


def test_mergebam_parameters(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    suffix = '-dedup'
    threads = 3
    index = 1
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.mergebam, ['-m', merge, '--suffix', suffix, '--threads', threads, '--index', index])
    assert result.exit_code == 0
    mb.merge_samples.assert_called_once_with(merge, suffix, threads, index)


def test_mergebam_mergenotexists(testdir, mock_testclass):
    merge = 'merge.txt'
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.mergebam, ['-m', merge])
    assert result.exit_code != 0
    mb.merge_samples.assert_not_called()


def test_merge_samples(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    mb.merge_sample = MagicMock()
    mb.merge_samples(merge)
    mb.merge_sample.assert_any_call('POLR2A', ['POLR2A_1', 'POLR2A_2'], '', 1)
    mb.merge_sample.assert_any_call('ASDURF', ['ASDURF_1', 'ASDURF_2'], '', 1)
    mb.merge_sample.assert_any_call('POLR1C', ['POLR1C_1', 'POLR1C_2'], '', 1)


def test_merge_samples_parameters(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    suffix = '-dedup'
    threads = 3
    mb.merge_sample = MagicMock()
    mb.merge_samples(merge, suffix, threads)
    mb.merge_sample.assert_any_call('POLR2A', ['POLR2A_1', 'POLR2A_2'], suffix, threads)
    mb.merge_sample.assert_any_call('ASDURF', ['ASDURF_1', 'ASDURF_2'], suffix, threads)
    mb.merge_sample.assert_any_call('POLR1C', ['POLR1C_1', 'POLR1C_2'], suffix, threads)


def test_merge_samples_second(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    mb.merge_sample = MagicMock()
    mb.merge_samples(merge, index=1)
    mb.merge_sample.assert_called_once_with('ASDURF', ['ASDURF_1', 'ASDURF_2'], '', 1)


def test_merge_sample(testdir, mock_testclass):
    merge = 'POLR2A'
    merge_bam = merge + '.bam'
    sample1 = merge + '_1'
    sample1_bam = sample1 + '.bam'
    sample2 = merge + '_2'
    sample2_bam = sample2 + '.bam'
    samples = [sample1, sample2]
    subprocess.run = MagicMock()
    mb.merge_sample(merge, samples)
    subprocess.run.assert_any_call(['samtools', 'index', sample1_bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'index', sample2_bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'merge', merge_bam, sample1_bam, sample2_bam], check=True)


def test_merge_sample_parameter(testdir, mock_testclass):
    suffix = '-dedup'
    merge = 'POLR2A'
    merge_bam = merge + suffix + '.bam'
    sample1 = merge + '_1'
    sample1_bam = sample1 + suffix + '.bam'
    sample2 = merge + '_2'
    sample2_bam = sample2 + suffix + '.bam'
    samples = [sample1, sample2]
    threads = 3
    subprocess.run = MagicMock()
    mb.merge_sample(merge, samples, suffix, threads)
    subprocess.run.assert_any_call(['samtools', 'index', '-@', str(threads - 1), sample1_bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'index', '-@', str(threads - 1), sample2_bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'merge', '--threads', str(threads - 1), merge_bam, sample1_bam, sample2_bam], check=True)
