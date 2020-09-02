import logging
import math
import os
from pathlib import Path
from shutil import copyfile
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import MergeBigwigs as mb
from seqtools.bed import Bed


@pytest.fixture
def mock_testclass():
    merge_samples = mb.merge_samples
    merge_sample = mb.merge_sample
    sort = Bed.sort
    bedgraph_to_bigwig = Bed.bedgraph_to_bigwig
    remove = os.remove
    yield merge_samples
    mb.merge_samples = merge_samples
    mb.merge_sample = merge_sample
    Bed.sort = sort
    Bed.bedgraph_to_bigwig = bedgraph_to_bigwig
    os.remove = remove
    

def create_file(*args, **kwargs):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_mergebw(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.mergebw, ['-m', merge, '--sizes', sizes])
    assert result.exit_code == 0
    mb.merge_samples.assert_called_once_with(merge, sizes, None)


def test_mergebw_parameters(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    index = 1
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.mergebw, ['-m', merge, '--sizes', sizes, '--index', index])
    assert result.exit_code == 0
    mb.merge_samples.assert_called_once_with(merge, sizes, index)


def test_mergebw_mergenotexists(testdir, mock_testclass):
    merge = 'merge.txt'
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.mergebw, ['-m', merge, '--sizes', sizes])
    assert result.exit_code != 0
    mb.merge_samples.assert_not_called()


def test_mergebw_sizesnotexists(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    sizes = 'sizes.txt'
    mb.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(mb.mergebw, ['-m', merge, '--sizes', sizes])
    assert result.exit_code != 0
    mb.merge_samples.assert_not_called()


def test_mergebw(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    mb.merge_sample = MagicMock()
    mb.merge_samples(merge, sizes)
    mb.merge_sample.assert_any_call('POLR2A', ['POLR2A_1', 'POLR2A_2'], sizes)
    mb.merge_sample.assert_any_call('ASDURF', ['ASDURF_1', 'ASDURF_2'], sizes)
    mb.merge_sample.assert_any_call('POLR1C', ['POLR1C_1', 'POLR1C_2'], sizes)


def test_mergebw_second(testdir, mock_testclass):
    merge = Path(__file__).parent.joinpath('merge.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    mb.merge_sample = MagicMock()
    mb.merge_samples(merge, sizes, 1)
    mb.merge_sample.assert_called_once_with('ASDURF', ['ASDURF_1', 'ASDURF_2'], sizes)


def test_merge_sample(testdir, mock_testclass):
    merge = 'POLR2A'
    merge_bw = merge + '.bw'
    sample1 = merge + '_1'
    sample1_bw = sample1 + '.bw'
    sample2 = merge + '_2'
    sample2_bw = sample2 + '.bw'
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    copyfile(Path(__file__).parent.joinpath('sample.bw'), sample1_bw)
    copyfile(Path(__file__).parent.joinpath('sample2.bw'), sample2_bw)
    Bed.sort = MagicMock()
    Bed.bedgraph_to_bigwig = MagicMock(side_effect=create_file(['-o', merge_bw]))
    os_remove = os.remove
    os.remove = MagicMock()
    mb.merge_sample(merge, [sample1, sample2], sizes)
    Bed.sort.assert_called_once_with(ANY, ANY)
    Bed.bedgraph_to_bigwig.assert_called_once_with(ANY, merge_bw, sizes)
    assert Bed.sort.call_args.args[1] == Bed.bedgraph_to_bigwig.call_args.args[0]
    assert os.path.exists(merge_bw)
    with open(Bed.sort.call_args.args[0], 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="POLR2A"\n'
        assert infile.readline() == 'chrI\t0\t1\t0\n'
        assert infile.readline() == 'chrI\t1\t2\t0\n'
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t2\t3\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.1, abs_tol=0.001), line.split('\t')[3]
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t3\t4\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.7, abs_tol=0.001), line.split('\t')[3]
        assert infile.readline() == 'chrI\t4\t5\t0\n'
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t5\t6\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.7, abs_tol=0.001), line.split('\t')[3]
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t6\t7\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.5, abs_tol=0.001), line.split('\t')[3]
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t7\t8\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.5, abs_tol=0.001), line.split('\t')[3]
        assert infile.readline() == 'chrI\t8\t9\t0\n'
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t9\t10\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.6, abs_tol=0.001), line.split('\t')[3]
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t10\t11\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.8, abs_tol=0.001), line.split('\t')[3]
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t11\t12\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.7, abs_tol=0.001), line.split('\t')[3]
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t12\t13\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.6, abs_tol=0.001), line.split('\t')[3]
        line = infile.readline()
        assert len(line.split('\t')) == 4
        assert line.startswith('chrI\t13\t14\t'), line
        assert math.isclose(float(line.split('\t')[3]), 0.6, abs_tol=0.001), line.split('\t')[3]
        assert infile.readline() == 'chrI\t14\t15\t0\n'
        assert infile.readline() == ''
    for remove_args in os.remove.call_args_list:
        os_remove(remove_args.args[0])
