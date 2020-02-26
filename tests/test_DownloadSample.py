import logging
import os
from pathlib import Path
import subprocess
from unittest.mock import MagicMock

import click
from click.testing import CliRunner
import pytest

from seqtools import DownloadSample as ds


def create_file(*args):
    for name in args:
        with open(name, 'w') as outfile:
            outfile.write('test')


def test_main(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    mem = '100MB'
    threads = 6
    subprocess.run = MagicMock(side_effect=[create_file('SRR8518913_1.fastq'), create_file('SRX5322424_1.fastq', 'SRX5322424_2.fastq'), create_file('SRR8518915_1.fastq', 'SRR8518915_2.fastq')])
    runner = CliRunner()
    result = runner.invoke(ds.main, ['-s', samples ])
    assert result.exit_code == 0
    subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRR8518913'], check=True)
    subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRX5322424'], check=True)
    subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRR8518915'], check=True)
    assert os.path.exists('POLR2A_1.fastq')
    assert not os.path.exists('POLR2A_2.fastq')
    assert os.path.exists('ASDURF_1.fastq')
    assert os.path.exists('ASDURF_2.fastq')
    assert os.path.exists('POLR1C_1.fastq')
    assert os.path.exists('POLR1C_2.fastq')


def test_main_allfast(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    mem = '200MB'
    threads = 2
    subprocess.run = MagicMock(side_effect=[create_file('SRR8518913_1.fastq'), create_file('SRX5322424_1.fastq', 'SRX5322424_2.fastq'), create_file('SRR8518915_1.fastq', 'SRR8518915_2.fastq')])
    runner = CliRunner()
    result = runner.invoke(ds.main, ['-s', samples, '--fast', '-e', str(threads), '-m', mem])
    assert result.exit_code == 0
    subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRR8518913'], check=True)
    subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRX5322424'], check=True)
    subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRR8518915'], check=True)
    assert os.path.exists('POLR2A_1.fastq')
    assert not os.path.exists('POLR2A_2.fastq')
    assert os.path.exists('ASDURF_1.fastq')
    assert os.path.exists('ASDURF_2.fastq')
    assert os.path.exists('POLR1C_1.fastq')
    assert os.path.exists('POLR1C_2.fastq')


def test_main_allslow(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    mem = '200MB'
    threads = 2
    subprocess.run = MagicMock(side_effect=[create_file('SRR8518913_1.fastq'), create_file('SRX5322424_1.fastq', 'SRX5322424_2.fastq'), create_file('SRR8518915_1.fastq', 'SRR8518915_2.fastq')])
    runner = CliRunner()
    result = runner.invoke(ds.main, ['-s', samples, '--slow', '-e', str(threads), '-m', mem])
    assert result.exit_code == 0
    subprocess.run.assert_any_call(['fastq-dump', '--split-files', 'SRR8518913'], check=True)
    subprocess.run.assert_any_call(['fastq-dump', '--split-files', 'SRX5322424'], check=True)
    subprocess.run.assert_any_call(['fastq-dump', '--split-files', 'SRR8518915'], check=True)
    assert os.path.exists('POLR2A_1.fastq')
    assert not os.path.exists('POLR2A_2.fastq')
    assert os.path.exists('ASDURF_1.fastq')
    assert os.path.exists('ASDURF_2.fastq')
    assert os.path.exists('POLR1C_1.fastq')
    assert os.path.exists('POLR1C_2.fastq')


def test_main_second(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    mem = '100MB'
    threads = 6
    subprocess.run = MagicMock(side_effect=create_file('SRX5322424_1.fastq', 'SRX5322424_2.fastq'))
    runner = CliRunner()
    result = runner.invoke(ds.main, ['-s', samples, '-i', 1])
    assert result.exit_code == 0
    subprocess.run.assert_called_once_with(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRX5322424'], check=True)
    assert os.path.exists('ASDURF_1.fastq')
    assert os.path.exists('ASDURF_2.fastq')


def test_download_singleend(testdir):
    srr = 'SRR8518913'
    mem = '100MB'
    threads = 1
    subprocess.run = MagicMock(side_effect=create_file(srr + '_1.fastq'))
    ds.download('POLR2A', srr, True, threads, mem)
    subprocess.run.assert_called_with(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, srr], check=True)
    assert os.path.exists('POLR2A_1.fastq')


def test_download_pairedend(testdir):
    srr = 'SRR8518913'
    mem = '100MB'
    threads = 1
    subprocess.run = MagicMock(side_effect=create_file(srr + '_1.fastq', srr + '_2.fastq'))
    ds.download('POLR2A', srr, True, threads, mem)
    subprocess.run.assert_called_with(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, srr], check=True)
    assert os.path.exists('POLR2A_1.fastq')
    assert os.path.exists('POLR2A_2.fastq')


def test_download_slow_singleend(testdir):
    srr = 'SRR8518913'
    mem = '100MB'
    threads = 1
    subprocess.run = MagicMock(side_effect=create_file(srr + '_1.fastq'))
    ds.download('POLR2A', srr, False, threads, mem)
    subprocess.run.assert_called_with(['fastq-dump', '--split-files', srr], check=True)
    assert os.path.exists('POLR2A_1.fastq')


def test_download_slow_pairedend(testdir):
    srr = 'SRR8518913'
    mem = '100MB'
    threads = 1
    subprocess.run = MagicMock(side_effect=create_file(srr + '_1.fastq', srr + '_2.fastq'))
    ds.download('POLR2A', srr, False, threads, mem)
    subprocess.run.assert_called_with(['fastq-dump', '--split-files', srr], check=True)
    assert os.path.exists('POLR2A_1.fastq')
    assert os.path.exists('POLR2A_2.fastq')


def test_download_failed(testdir):
    srr = 'SRR8518913'
    mem = '100MB'
    threads = 1
    subprocess.run = MagicMock(side_effect=subprocess.CalledProcessError('Could not download file', ['test']))
    with pytest.raises(subprocess.CalledProcessError):
        ds.download('POLR2A', srr, True, threads, mem)
    subprocess.run.assert_called_with(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, srr], check=True)
    assert not os.path.exists('POLR2A_1.fastq')
    assert not os.path.exists('POLR2A_2.fastq')
