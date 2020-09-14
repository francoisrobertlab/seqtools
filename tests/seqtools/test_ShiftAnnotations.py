import logging
import os
from pathlib import Path
from shutil import copyfile
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import ShiftAnnotations as sa
from seqtools.txt import Parser


@pytest.fixture
def mock_testclass():
    shift_annotations_samples = sa.shift_annotations_samples
    shift_annotations_sample = sa.shift_annotations_sample
    first = Parser.first
    yield
    sa.shift_annotations_samples = shift_annotations_samples
    sa.shift_annotations_sample = shift_annotations_sample
    Parser.first = first
    
    
def create_file(*args, **kwargs):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_shiftannotations(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    genome = 'sacCer3.chrom.sizes'
    copyfile(Path(__file__).parent.joinpath('sizes.txt'), genome)
    sa.shift_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(sa.shiftannotations, ['-s', samples, '-g', genome])
    assert result.exit_code == 0
    sa.shift_annotations_samples.assert_called_once_with(samples, genome, None, ())


def test_shiftannotations_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    genome = 'sacCer3.chrom.sizes'
    copyfile(Path(__file__).parent.joinpath('sizes.txt'), genome)
    minus = '2'
    plus = '-2'
    index = 1
    sa.shift_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(sa.shiftannotations, ['-s', samples, '-g', genome, '-m', minus, '-p', plus, '-i', index])
    assert result.exit_code == 0
    sa.shift_annotations_samples.assert_called_once_with(samples, genome, index, ('-m', minus, '-p', plus,))


def test_shiftannotations_second(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    genome = 'sacCer3.chrom.sizes'
    copyfile(Path(__file__).parent.joinpath('sizes.txt'), genome)
    index = 1
    sa.shift_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(sa.shiftannotations, ['-s', samples, '-g', genome, '-i', index])
    assert result.exit_code == 0
    sa.shift_annotations_samples.assert_called_once_with(samples, genome, index, ())


def test_shiftannotations_samplesnotexists(testdir, mock_testclass):
    samples = 'samples.txt'
    genome = 'sacCer3.chrom.sizes'
    copyfile(Path(__file__).parent.joinpath('sizes.txt'), genome)
    sa.shift_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(sa.shiftannotations, ['-s', samples, '-g', genome])
    assert result.exit_code != 0
    sa.shift_annotations_samples.assert_not_called()


def test_shiftannotations_genomenotexists(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    genome = 'sacCer3.chrom.sizes'
    sa.shift_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(sa.shiftannotations, ['-s', samples, '-g', genome])
    assert result.exit_code != 0
    sa.shift_annotations_samples.assert_not_called()


def test_shift_annotations_samples(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    genome = 'sacCer3.chrom.sizes'
    Parser.first = MagicMock(return_value=samples)
    sa.shift_annotations_sample = MagicMock()
    sa.shift_annotations_samples(samples_file)
    Parser.first.assert_called_once_with(samples_file)
    for sample in samples:
        sa.shift_annotations_sample.assert_any_call(sample, genome, ())


def test_shift_annotations_samples_parameters(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    genome = 'sacCer3.chrom.sizes'
    minus = '2'
    plus = '-2'
    Parser.first = MagicMock(return_value=samples)
    sa.shift_annotations_sample = MagicMock()
    sa.shift_annotations_samples(samples_file, bedtools_args=('-m', minus, '-p', plus,))
    Parser.first.assert_called_once_with(samples_file)
    for sample in samples:
        sa.shift_annotations_sample.assert_any_call(sample, genome, ('-m', minus, '-p', plus,))


def test_shift_annotations_samples_second(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    genome = 'sacCer3.chrom.sizes'
    Parser.first = MagicMock(return_value=samples)
    sa.shift_annotations_sample = MagicMock()
    sa.shift_annotations_samples(samples_file, index=1)
    Parser.first.assert_called_once_with(samples_file)
    sa.shift_annotations_sample.assert_any_call(samples[1], genome, ())

    
def test_shift_annotations_sample(testdir, mock_testclass):
    sample = 'POLR2A'
    genome = 'sacCer3.chrom.sizes'
    src_bed = Path(__file__).parent.joinpath('sample.bed')
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    copyfile(src_bed, bed)
    subprocess.run = MagicMock(side_effect=create_file)
    sa.shift_annotations_sample(sample, genome)
    subprocess.run.assert_called_once_with(['bedtools', 'shift', '-i', bed, '-g', genome], stdout=ANY, check=True)
    assert os.path.exists(forcov)
    with open(forcov, 'r') as infile:
        assert 'test' == infile.readline()


def test_shift_annotations_sample_parameters(testdir, mock_testclass):
    sample = 'POLR2A'
    genome = 'sacCer3.chrom.sizes'
    minus = '2'
    plus = '-2'
    src_bed = Path(__file__).parent.joinpath('sample.bed')
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    copyfile(src_bed, bed)
    subprocess.run = MagicMock(side_effect=create_file)
    sa.shift_annotations_sample(sample, genome, ('-m', minus, '-p', plus,))
    subprocess.run.assert_called_once_with(['bedtools', 'shift', '-i', bed, '-g', genome, '-m', minus, '-p', plus], stdout=ANY, check=True)
    assert os.path.exists(forcov)
    with open(forcov, 'r') as infile:
        assert 'test' == infile.readline()
