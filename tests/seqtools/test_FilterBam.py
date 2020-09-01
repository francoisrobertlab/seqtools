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
    filter_bam = fb.filter_bam
    filter_bam_sample = fb.filter_bam_sample
    filter_mapped = fb.filter_mapped
    remove_duplicates = fb.remove_duplicates
    sort = fb.sort
    run = subprocess.run
    yield
    fb.filter_bam = filter_bam
    fb.filter_bam_sample = filter_bam_sample
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
    fb.filter_bam = MagicMock()
    runner = CliRunner()
    result = runner.invoke(fb.filterbam, ['-s', samples])
    assert result.exit_code == 0
    fb.filter_bam.assert_called_once_with(samples, True, True, threads, None)


def test_filterbam_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    index = 1
    fb.filter_bam = MagicMock()
    runner = CliRunner()
    result = runner.invoke(fb.filterbam, ['-s', samples, '--unpaired', '--no-dedup', '--threads', threads, '--index', index])
    assert result.exit_code == 0
    fb.filter_bam.assert_called_once_with(samples, False, False, threads, index)


def test_filter_bam(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    fb.filter_bam_sample = MagicMock()
    fb.filter_bam(samples)
    fb.filter_bam_sample.assert_any_call('POLR2A', True, True, None)
    fb.filter_bam_sample.assert_any_call('ASDURF', True, True, None)
    fb.filter_bam_sample.assert_any_call('POLR1C', True, True, None)


def test_filter_bam_all_single_threads(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    fb.filter_bam_sample = MagicMock()
    fb.filter_bam(samples, False, True, threads)
    fb.filter_bam_sample.assert_any_call('POLR2A', False, True, threads)
    fb.filter_bam_sample.assert_any_call('ASDURF', False, True, threads)
    fb.filter_bam_sample.assert_any_call('POLR1C', False, True, threads)


def test_filter_bam_second_single_threads(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    fb.filter_bam_sample = MagicMock()
    fb.filter_bam(samples, False, True, threads, 1)
    fb.filter_bam_sample.assert_called_once_with('ASDURF', False, True, threads)


def test_filter_bam_all_paired_threads(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    fb.filter_bam_sample = MagicMock()
    fb.filter_bam(samples, True, True, threads)
    fb.filter_bam_sample.assert_any_call('POLR2A', True, True, threads)
    fb.filter_bam_sample.assert_any_call('ASDURF', True, True, threads)
    fb.filter_bam_sample.assert_any_call('POLR1C', True, True, threads)


def test_filter_bam_second_paired_threads(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 2
    fb.filter_bam_sample = MagicMock()
    fb.filter_bam(samples, True, True, threads, 1)
    fb.filter_bam_sample.assert_called_once_with('ASDURF', True, True, threads)


def test_filter_bam_sample_single(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    bam_filtered = sample + '-filtered.bam'
    bam_dedup = sample + '-dedup.bam'
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock(create_file(['-o', bam_dedup]))
    fb.filter_bam_sample(sample, False, True)
    fb.filter_mapped.assert_called_with(bam, bam_filtered, False, None)
    fb.remove_duplicates.assert_called_with(bam_filtered, bam_dedup, None)
    assert os.path.exists(bam_filtered)
    assert os.path.exists(bam_dedup)


def test_filter_bam_sample_single_nodedup(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    bam_filtered = sample + '-filtered.bam'
    bam_dedup = sample + '-dedup.bam'
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock()
    fb.filter_bam_sample(sample, False, False)
    fb.filter_mapped.assert_called_with(bam, bam_filtered, False, None)
    fb.remove_duplicates.assert_not_called()
    assert os.path.exists(bam_filtered)
    assert not os.path.exists(bam_dedup)


def test_filter_bam_sample_single_threads(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    bam_filtered = sample + '-filtered.bam'
    bam_dedup = sample + '-dedup.bam'
    threads = 3
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock()
    fb.filter_bam_sample(sample, False, False, threads)
    fb.filter_mapped.assert_called_with(bam, bam_filtered, False, threads)
    fb.remove_duplicates.assert_not_called()
    assert os.path.exists(bam_filtered)
    assert not os.path.exists(bam_dedup)


def test_filter_bam_sample_paired(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    bam_filtered = sample + '-filtered.bam'
    bam_dedup = sample + '-dedup.bam'
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock(create_file(['-o', bam_dedup]))
    fb.filter_bam_sample(sample, True, True)
    fb.filter_mapped.assert_called_with(bam, bam_filtered, True, None)
    fb.remove_duplicates.assert_called_with(bam_filtered, bam_dedup, None)
    assert os.path.exists(bam_filtered)
    assert os.path.exists(bam_dedup)


def test_filter_bam_sample_paired_nodedup(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    bam_filtered = sample + '-filtered.bam'
    bam_dedup = sample + '-dedup.bam'
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock()
    fb.filter_bam_sample(sample, True, False)
    fb.filter_mapped.assert_called_with(bam, bam_filtered, True, None)
    fb.remove_duplicates.assert_not_called()
    assert os.path.exists(bam_filtered)
    assert not os.path.exists(bam_dedup)


def test_filter_bam_sample_paired_threads(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    bam_filtered = sample + '-filtered.bam'
    bam_dedup = sample + '-dedup.bam'
    threads = 3
    fb.filter_mapped = MagicMock(create_file(['-o', bam_filtered]))
    fb.remove_duplicates = MagicMock(create_file(['-o', bam_dedup]))
    fb.filter_bam_sample(sample, True, True, threads)
    fb.filter_mapped.assert_called_with(bam, bam_filtered, True, threads)
    fb.remove_duplicates.assert_called_with(bam_filtered, bam_dedup, threads)
    assert os.path.exists(bam_filtered)
    assert os.path.exists(bam_dedup)


def test_filter_mapped_single(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    temp = bam + '.tmp'
    output = 'POLR2A-out.bam'
    subprocess.run = MagicMock(side_effect=create_file)
    fb.filter_mapped(bam, output, False)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-F', '4', '-o', temp, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', output, temp], check=True)


def test_filter_mapped_single_threads(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    temp = bam + '.tmp'
    output = 'POLR2A-out.bam'
    threads = 3
    subprocess.run = MagicMock(side_effect=create_file)
    fb.filter_mapped(bam, output, False, threads)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-F', '4', '--threads', str(threads - 1), '-o', temp, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '--threads', str(threads - 1), '-o', output, temp], check=True)


def test_filter_mapped_single_singlethread(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    temp = bam + '.tmp'
    output = 'POLR2A-out.bam'
    threads = 1
    subprocess.run = MagicMock(side_effect=create_file)
    fb.filter_mapped(bam, output, False, threads)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-F', '4', '-o', temp, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', output, temp], check=True)


def test_filter_mapped_paired(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    temp = bam + '.tmp'
    output = 'POLR2A-out.bam'
    subprocess.run = MagicMock(side_effect=create_file)
    fb.filter_mapped(bam, output, True)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-f', '2', '-o', temp, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', output, temp], check=True)


def test_filter_mapped_paired_threads(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    temp = bam + '.tmp'
    output = 'POLR2A-out.bam'
    threads = 3
    subprocess.run = MagicMock(side_effect=create_file)
    fb.filter_mapped(bam, output, True, threads)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-f', '2', '--threads', str(threads - 1), '-o', temp, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '--threads', str(threads - 1), '-o', output, temp], check=True)


def test_filter_mapped_paired_singlethread(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    temp = bam + '.tmp'
    output = 'POLR2A-out.bam'
    threads = 1
    subprocess.run = MagicMock(side_effect=create_file)
    fb.filter_mapped(bam, output, True, threads)
    subprocess.run.assert_any_call(['samtools', 'view', '-b', '-F', '2048', '-f', '2', '-o', temp, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', output, temp], check=True)


def test_remove_duplicates(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    sort_bam = bam + '.fixsort'
    fix = bam + '.fix'
    dedup_tmp = bam + '.dedup.tmp'
    sort = bam + '.sort'
    output = 'POLR2A-out.bam'
    subprocess.run = MagicMock(side_effect=[create_file(['-o', sort_bam]), create_file(['-o', fix]), create_file(['-o', sort]), create_file(['-o', dedup_tmp]), create_file(['-o', output])])
    fb.remove_duplicates(bam, output)
    subprocess.run.assert_any_call(['samtools', 'sort', '-n', '-o', sort_bam, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'fixmate', '-m', sort_bam, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', sort, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'markdup', '-r', sort, dedup_tmp], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', output, dedup_tmp], check=True)
    assert not os.path.exists(sort_bam)
    assert not os.path.exists(fix)
    assert not os.path.exists(dedup_tmp)
    assert not os.path.exists(sort)


def test_remove_duplicates_threads(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    sort_bam = bam + '.fixsort'
    fix = bam + '.fix'
    dedup_tmp = bam + '.dedup.tmp'
    sort = bam + '.sort'
    output = 'POLR2A-out.bam'
    threads = 3
    subprocess.run = MagicMock(side_effect=[create_file(['-o', sort_bam]), create_file(['-o', fix]), create_file(['-o', sort]), create_file(['-o', dedup_tmp]), create_file(['-o', output])])
    fb.remove_duplicates(bam, output, threads)
    subprocess.run.assert_any_call(['samtools', 'sort', '-n', '--threads', str(threads - 1), '-o', sort_bam, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'fixmate', '-m', '--threads', str(threads - 1), sort_bam, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '--threads', str(threads - 1), '-o', sort, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'markdup', '-r', '--threads', str(threads - 1), sort, dedup_tmp], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '--threads', str(threads - 1), '-o', output, dedup_tmp], check=True)
    assert not os.path.exists(sort_bam)
    assert not os.path.exists(fix)
    assert not os.path.exists(dedup_tmp)
    assert not os.path.exists(sort)


def test_remove_duplicates_singlethread(testdir, mock_testclass):
    bam = 'POLR2A.bam'
    sort_bam = bam + '.fixsort'
    fix = bam + '.fix'
    dedup_tmp = bam + '.dedup.tmp'
    sort = bam + '.sort'
    output = 'POLR2A-out.bam'
    threads = 1
    subprocess.run = MagicMock(side_effect=[create_file(['-o', sort_bam]), create_file(['-o', fix]), create_file(['-o', sort]), create_file(['-o', dedup_tmp]), create_file(['-o', output])])
    fb.remove_duplicates(bam, output, threads)
    subprocess.run.assert_any_call(['samtools', 'sort', '-n', '-o', sort_bam, bam], check=True)
    subprocess.run.assert_any_call(['samtools', 'fixmate', '-m', sort_bam, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', sort, fix], check=True)
    subprocess.run.assert_any_call(['samtools', 'markdup', '-r', sort, dedup_tmp], check=True)
    subprocess.run.assert_any_call(['samtools', 'sort', '-o', output, dedup_tmp], check=True)
    assert not os.path.exists(sort_bam)
    assert not os.path.exists(fix)
    assert not os.path.exists(dedup_tmp)
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
