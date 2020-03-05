import logging
import os
from pathlib import Path
from shutil import copyfile
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import Merge as mb
from seqtools.bed import Bed


@pytest.fixture
def mock_testclass():
    merge_samples = mb.merge_samples
    yield merge_samples
    mb.merge_samples = merge_samples
    

def create_file(*args, **kwargs):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_merge(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.merge, ['-m', merge])
    assert result.exit_code == 0
    mb.merge_samples.assert_any_call('POLR2A', ['POLR2A_1', 'POLR2A_2'])
    mb.merge_samples.assert_any_call('ASDURF', ['ASDURF_1', 'ASDURF_2'])
    mb.merge_samples.assert_any_call('POLR1C', ['POLR1C_1', 'POLR1C_2'])


def test_merge_second(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.merge, ['-m', merge, '-i', '1'])
    assert result.exit_code == 0
    mb.merge_samples.assert_called_once_with('ASDURF', ['ASDURF_1', 'ASDURF_2'])


def test_merge_mergenotexists(testdir, mock_testclass):
    merge = 'merge.txt'
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.merge, ['-m', merge])
    assert result.exit_code != 0
    mb.merge_samples.assert_not_called()


def test_merge_samples(testdir):
    merge = 'POLR2A'
    merge_bed_tmp = merge + '-tmp.bed'
    merge_bed = merge + '.bed'
    sample1 = merge + '_1'
    sample1_bed = sample1 + '.bed'
    sample2 = merge + '_2'
    sample2_bed = sample2 + '.bed'
    copyfile(Path(__file__).parent.joinpath('sample.bed'), sample1_bed)
    copyfile(Path(__file__).parent.joinpath('sample2.bed'), sample2_bed)
    Bed.sort = MagicMock(side_effect=create_file(['-o', merge_bed]))
    os.remove = MagicMock()
    mb.merge_samples(merge, [sample1, sample2])
    Bed.sort.assert_called_once_with(merge_bed_tmp, merge_bed)
    os.remove.assert_called_once_with(merge_bed_tmp)
    assert os.path.exists(merge_bed)
    with open(merge_bed_tmp, 'r') as infile:
        assert infile.readline() == 'chr1\t100\t150\ttest1\t1\t+\n'
        assert infile.readline() == 'chr2\t400\t450\ttest2\t2\t+\n'
        assert infile.readline() == 'chr3\t500\t650\ttest3\t3\t+\n'
        assert infile.readline() == 'chr4\t800\t750\ttest4\t4\t+\n'
        assert infile.readline() == 'chr5\t100\t150\ttest5\t1\t-\n'
        assert infile.readline() == 'chr6\t400\t450\ttest6\t2\t-\n'
        assert infile.readline() == 'chr7\t500\t650\ttest7\t3\t-\n'
        assert infile.readline() == 'chr8\t800\t750\ttest8\t4\t-\n'
        assert infile.readline() == 'chr1\t200\t250\ttest1\t1\t+\n'
        assert infile.readline() == 'chr2\t300\t350\ttest2\t2\t+\n'
        assert infile.readline() == 'chr3\t600\t550\ttest3\t3\t+\n'
        assert infile.readline() == 'chr4\t700\t850\ttest4\t4\t+\n'
        assert infile.readline() == 'chr5\t200\t250\ttest5\t1\t-\n'
        assert infile.readline() == 'chr6\t300\t350\ttest6\t2\t-\n'
        assert infile.readline() == 'chr7\t600\t550\ttest7\t3\t-\n'
        assert infile.readline() == 'chr8\t700\t850\ttest8\t4\t-\n'
        assert infile.readline() == ''
