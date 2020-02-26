import logging
import os
from pathlib import Path
import subprocess
from unittest.mock import MagicMock

import click
from click.testing import CliRunner
import pytest

from seqtools import DownloadSample as ds


@pytest.fixture
def mock_testclass():
    download = ds.download
    yield download
    ds.download = download
    
    
def create_file(*args):
    for name in args:
        with open(name, 'w') as outfile:
            outfile.write('test')


def test_main(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    mem = '100MB'
    threads = 6
    ds.download = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ds.main, ['-s', samples ])
    assert result.exit_code == 0
    ds.download.assert_any_call('POLR2A', 'SRR8518913', True, threads, mem)
    ds.download.assert_any_call('ASDURF', 'SRX5322424', True, threads, mem)
    ds.download.assert_any_call('POLR1C', 'SRR8518915', True, threads, mem)


def test_main_allfast(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    mem = '200MB'
    threads = 2
    ds.download = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ds.main, ['-s', samples, '--fast', '-e', str(threads), '-m', mem])
    assert result.exit_code == 0
    ds.download.assert_any_call('POLR2A', 'SRR8518913', True, threads, mem)
    ds.download.assert_any_call('ASDURF', 'SRX5322424', True, threads, mem)
    ds.download.assert_any_call('POLR1C', 'SRR8518915', True, threads, mem)


def test_main_allslow(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    mem = '200MB'
    threads = 2
    ds.download = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ds.main, ['-s', samples, '--slow', '-e', str(threads), '-m', mem])
    assert result.exit_code == 0
    ds.download.assert_any_call('POLR2A', 'SRR8518913', False, threads, mem)
    ds.download.assert_any_call('ASDURF', 'SRX5322424', False, threads, mem)
    ds.download.assert_any_call('POLR1C', 'SRR8518915', False, threads, mem)


def test_main_second(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    mem = '100MB'
    threads = 6
    ds.download = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ds.main, ['-s', samples, '-i', 1])
    assert result.exit_code == 0
    ds.download.assert_called_once_with('ASDURF', 'SRX5322424', True, threads, mem)


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
