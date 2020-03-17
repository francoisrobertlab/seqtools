import logging
import os
from pathlib import Path
from shutil import copyfile
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import Split as s
from seqtools.bed import Bed


@pytest.fixture
def mock_testclass():
    split_samples = s.split_samples
    split_sample = s.split_sample
    sort_bysize = Bed.sort_bysize
    sort = Bed.sort
    remove = os.remove
    yield
    s.split_samples = split_samples
    s.split_sample = split_sample
    Bed.sort_bysize = sort_bysize
    Bed.sort = sort
    os.remove = remove
   

def write_file(file, annotations):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_split(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    s.split_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(s.split, ['-s', samples])
    assert result.exit_code == 0
    s.split_samples.assert_called_once_with(samples, None, 10, 100, 500)


def test_split_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    s.split_samples = MagicMock()
    binlength = 20
    binminlength = 200
    binmaxlength = 400
    index = 1
    s.split_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(s.split, ['-s', samples, '--binLength', binlength, '--binMinLength', binminlength, '--binMaxLength', binmaxlength, '--index', index])
    assert result.exit_code == 0
    s.split_samples.assert_called_once_with(samples, index, binlength, binminlength, binmaxlength)


def test_split_samplesnotexists(testdir, mock_testclass):
    samples = 'samples.txt'
    s.split_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(s.split, ['-s', samples])
    assert result.exit_code != 0
    s.split_samples.assert_not_called()


def test_split_samples(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    s.split_sample = MagicMock()
    s.split_samples(samples)
    s.split_sample.assert_any_call('POLR2A', 10, 100, 500)
    s.split_sample.assert_any_call('ASDURF', 10, 100, 500)
    s.split_sample.assert_any_call('POLR1C', 10, 100, 500)


def test_split_samples_second(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    s.split_sample = MagicMock()
    s.split_samples(samples, 1)
    s.split_sample.assert_called_once_with('ASDURF', 10, 100, 500)


def test_split_samples_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    s.split_sample = MagicMock()
    binlength = 20
    binminlength = 200
    binmaxlength = 400
    s.split_samples(samples, binlength=binlength, binminlength=binminlength, binmaxlength=binmaxlength)
    s.split_sample.assert_any_call('POLR2A', binlength, binminlength, binmaxlength)
    s.split_sample.assert_any_call('ASDURF', binlength, binminlength, binmaxlength)
    s.split_sample.assert_any_call('POLR1C', binlength, binminlength, binmaxlength)


def test_split_samples_second_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    s.split_sample = MagicMock()
    binlength = 20
    binminlength = 200
    binmaxlength = 400
    s.split_samples(samples, 1, binlength, binminlength, binmaxlength)
    s.split_sample.assert_called_once_with('ASDURF', binlength, binminlength, binmaxlength)


def test_split_sample(testdir, mock_testclass):
    sample = 'POLR2A'
    bed = sample + '.bed'
    sort = sample + '-sort.bed'
    sort_copy = Path(__file__).parent.joinpath('sample-split.bed')
    copyfile(sort_copy, sort)
    Bed.sort_bysize = MagicMock()
    Bed.sort = MagicMock()
    os.remove = MagicMock()
    binlength = 10
    binminlength = 100
    binmaxlength = 130
    s.split_sample(sample, binlength, binminlength, binmaxlength)
    Bed.sort_bysize.assert_called_once_with(bed, sort)
    Bed.sort.assert_any_call(sample + '-100-110-tmp.bed', sample + '-100-110.bed')
    Bed.sort.assert_any_call(sample + '-110-120-tmp.bed', sample + '-110-120.bed')
    Bed.sort.assert_any_call(sample + '-120-130-tmp.bed', sample + '-120-130.bed')
    os.remove.assert_any_call(sort)
    os.remove.assert_any_call(sample + '-100-110-tmp.bed')
    os.remove.assert_any_call(sample + '-110-120-tmp.bed')
    os.remove.assert_any_call(sample + '-120-130-tmp.bed')
    with open(sample + '-100-110-tmp.bed', 'r') as infile:
        assert infile.readline() == 'chr4\t800\t900\ttest4\t4\t+\n'
        assert infile.readline() == ''
    with open(sample + '-110-120-tmp.bed', 'r') as infile:
        assert infile.readline() == 'chr8\t800\t910\ttest8\t4\t-\n'
        assert infile.readline() == ''
    with open(sample + '-120-130-tmp.bed', 'r') as infile:
        assert infile.readline() == 'chr5\t100\t220\ttest5\t1\t-\n'
        assert infile.readline() == 'chr1\t100\t229\ttest1\t1\t+\n'
        assert infile.readline() == ''
