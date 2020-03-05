import logging
import os
from pathlib import Path
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import FilterBam as fb


@pytest.fixture
def mock_testclass():
    filterbam_sample = fb.filterbam_sample
    filter_mapped = fb.filter_mapped
    remove_duplicates = fb.remove_duplicates
    sort = fb.sort
    run = subprocess.run
    yield
    fb.filterbam_sample = filterbam_sample
    fb.filter_mapped = filter_mapped
    fb.remove_duplicates = remove_duplicates
    fb.sort = sort
    subprocess.run = run
    
    
def create_file(*args, **kwargs):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_filterbam(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 1
    fb.filterbam_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(fb.filterbam, ['-s', samples ])
    assert result.exit_code == 0
    fb.filterbam_sample.assert_any_call('POLR2A', True, threads)
    fb.filterbam_sample.assert_any_call('ASDURF', True, threads)
    fb.filterbam_sample.assert_any_call('POLR1C', True, threads)


def test_filterbam_all_single_threads(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    fb.filterbam_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(fb.filterbam, ['-s', samples, '--unpaired', '--threads', threads ])
    assert result.exit_code == 0
    fb.filterbam_sample.assert_any_call('POLR2A', False, threads)
    fb.filterbam_sample.assert_any_call('ASDURF', False, threads)
    fb.filterbam_sample.assert_any_call('POLR1C', False, threads)


def test_filterbam_second_single_threads(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    fb.filterbam_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(fb.filterbam, ['-s', samples, '--unpaired', '--threads', threads, '-i', 1 ])
    assert result.exit_code == 0
    fb.filterbam_sample.assert_called_once_with('ASDURF', False, threads)


def test_filterbam_all_paired_threads(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    fb.filterbam_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(fb.filterbam, ['-s', samples, '--paired', '--threads', threads ])
    assert result.exit_code == 0
    fb.filterbam_sample.assert_any_call('POLR2A', True, threads)
    fb.filterbam_sample.assert_any_call('ASDURF', True, threads)
    fb.filterbam_sample.assert_any_call('POLR1C', True, threads)


def test_filterbam_second_paired_threads(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    fb.filterbam_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(fb.filterbam, ['-s', samples, '--paired', '--threads', threads, '-i', 1 ])
    assert result.exit_code == 0
    fb.filterbam_sample.assert_called_once_with('ASDURF', True, threads)


def test_filterbam_sample_single(testdir, mock_testclass):
    sample = 'POLR2A'
    bam_raw = sample + '-raw.bam'
    bam_filtered = sample + '-filtered.bam'
    bam = sample + '.bam'
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock()
    fb.sort = MagicMock(create_file(['-o', bam]))
    fb.filterbam_sample(sample, False)
    fb.filter_mapped.assert_called_with(bam_raw, bam_filtered, False, None)
    fb.remove_duplicates.assert_not_called()
    fb.sort.assert_called_with(bam_filtered, bam, None)
    assert os.path.exists(bam_filtered)
    assert os.path.exists(bam)


def test_filterbam_sample_single_threads(testdir, mock_testclass):
    sample = 'POLR2A'
    bam_raw = sample + '-raw.bam'
    bam_filtered = sample + '-filtered.bam'
    bam = sample + '.bam'
    threads = 3
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock()
    fb.sort = MagicMock(create_file(['-o', bam]))
    fb.filterbam_sample(sample, False, threads)
    fb.filter_mapped.assert_called_with(bam_raw, bam_filtered, False, threads)
    fb.remove_duplicates.assert_not_called()
    fb.sort.assert_called_with(bam_filtered, bam, threads)
    assert os.path.exists(bam_filtered)
    assert os.path.exists(bam)


def test_filterbam_sample_paired(testdir, mock_testclass):
    sample = 'POLR2A'
    bam_raw = sample + '-raw.bam'
    bam_filtered = sample + '-filtered.bam'
    bam_dedup = sample + '-dedup.bam'
    bam = sample + '.bam'
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock(create_file(['-o', bam_dedup]))
    fb.sort = MagicMock(create_file(['-o', bam]))
    fb.filterbam_sample(sample, True)
    fb.filter_mapped.assert_called_with(bam_raw, bam_filtered, True, None)
    fb.remove_duplicates.assert_called_with(bam_filtered, bam_dedup, None)
    fb.sort.assert_called_with(bam_dedup, bam, None)
    assert os.path.exists(bam_filtered)
    assert os.path.exists(bam_dedup)
    assert os.path.exists(bam)


def test_filterbam_sample_paired_threads(testdir, mock_testclass):
    sample = 'POLR2A'
    bam_raw = sample + '-raw.bam'
    bam_filtered = sample + '-filtered.bam'
    bam_dedup = sample + '-dedup.bam'
    bam = sample + '.bam'
    threads = 3
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock(create_file(['-o', bam_dedup]))
    fb.sort = MagicMock(create_file(['-o', bam]))
    fb.filterbam_sample(sample, True, threads)
    fb.filter_mapped.assert_called_with(bam_raw, bam_filtered, True, threads)
    fb.remove_duplicates.assert_called_with(bam_filtered, bam_dedup, threads)
    fb.sort.assert_called_with(bam_dedup, bam, threads)
    assert os.path.exists(bam_filtered)
    assert os.path.exists(bam_dedup)
    assert os.path.exists(bam)


def test_filter_mapped_single(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    subprocess.run = MagicMock()
    fb.filter_mapped(bam, output, False)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-F', '4', '-o', output, bam], check=True)


def test_filter_mapped_single_threads(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    threads = 3
    subprocess.run = MagicMock()
    fb.filter_mapped(bam, output, False, threads)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-F', '4', '--threads', str(threads - 1), '-o', output, bam], check=True)


def test_filter_mapped_single_singlethread(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    threads = 1
    subprocess.run = MagicMock()
    fb.filter_mapped(bam, output, False, threads)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-F', '4', '-o', output, bam], check=True)


def test_filter_mapped_paired(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    subprocess.run = MagicMock()
    fb.filter_mapped(bam, output, True)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-f', '2', '-o', output, bam], check=True)


def test_filter_mapped_paired_threads(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    threads = 3
    subprocess.run = MagicMock()
    fb.filter_mapped(bam, output, True, threads)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-f', '2', '--threads', str(threads - 1), '-o', output, bam], check=True)


def test_filter_mapped_paired_singlethread(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    threads = 1
    subprocess.run = MagicMock()
    fb.filter_mapped(bam, output, True, threads)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-f', '2', '-o', output, bam], check=True)


def test_remove_duplicates(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    fix = bam + '.fix'
    sort = bam + '.sort'
    output = 'POLR2A-out.bam'
    subprocess.run = MagicMock(side_effect=[create_file(['-o', fix]), create_file(['-o', sort]), create_file(['-o', output])])
    fb.remove_duplicates(bam, output)
    subprocess.run.assert_any_call(['samtools', 'fixmate', '-m', bam, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', sort, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'markdup', '-r', sort, output], check=True)
    assert not os.path.exists(fix)
    assert not os.path.exists(sort)


def test_remove_duplicates_threads(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    fix = bam + '.fix'
    sort = bam + '.sort'
    output = 'POLR2A-out.bam'
    threads = 3
    subprocess.run = MagicMock(side_effect=[create_file(['-o', fix]), create_file(['-o', sort]), create_file(['-o', output])])
    fb.remove_duplicates(bam, output, threads)
    subprocess.run.assert_any_call(['samtools', 'fixmate', '-m', '--threads', str(threads - 1), bam, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '--threads', str(threads - 1), '-o', sort, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'markdup', '-r', '--threads', str(threads - 1), sort, output], check=True)
    assert not os.path.exists(fix)
    assert not os.path.exists(sort)


def test_remove_duplicates_singlethread(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    fix = bam + '.fix'
    sort = bam + '.sort'
    output = 'POLR2A-out.bam'
    threads = 1
    subprocess.run = MagicMock(side_effect=[create_file(['-o', fix]), create_file(['-o', sort]), create_file(['-o', output])])
    fb.remove_duplicates(bam, output, threads)
    subprocess.run.assert_any_call(['samtools', 'fixmate', '-m', bam, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', sort, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'markdup', '-r', sort, output], check=True)
    assert not os.path.exists(fix)
    assert not os.path.exists(sort)


def test_sort(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    subprocess.run = MagicMock()
    fb.sort(bam, output)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', output, bam], check=True)


def test_sort_threads(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    threads = 3
    subprocess.run = MagicMock()
    fb.sort(bam, output, threads)
    subprocess.run.assert_any_call(['samtools', 'sort', '--threads', str(threads - 1), '-o', output, bam], check=True)


def test_sort_singlethread(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    output = 'POLR2A-out.bam'
    threads = 1
    subprocess.run = MagicMock()
    fb.sort(bam, output, threads)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', output, bam], check=True)
