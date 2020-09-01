import logging
import os
from pathlib import Path
from shutil import copyfile
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest
from seqtools import MoveAnnotations as ma
from seqtools.txt import Parser


@pytest.fixture
def mock_testclass():
    moveannotations_samples = ma.moveannotations_samples
    moveannotations_sample = ma.moveannotations_sample
    first = Parser.first
    yield
    ma.moveannotations_samples = moveannotations_samples
    ma.moveannotations_sample = moveannotations_sample
    Parser.first = first
    
    
def test_moveannotations(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    distance = 2
    ma.moveannotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ma.moveannotations, ['-s', samples, '-d', distance])
    assert result.exit_code == 0
    ma.moveannotations_samples.assert_called_once_with(samples, distance, True, True, None)


def test_moveannotations_parameters(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    distance = 2
    ma.moveannotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ma.moveannotations, ['-s', samples, '-d', distance, '--keep-negatives', '--forward-negative'])
    assert result.exit_code == 0
    ma.moveannotations_samples.assert_called_once_with(samples, distance, False, False, None)


def test_moveannotations_second(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    distance = 2
    index = 1
    ma.moveannotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ma.moveannotations, ['-s', samples, '-d', 2, '-i', index])
    assert result.exit_code == 0
    ma.moveannotations_samples.assert_called_once_with(samples, distance, True, True, index)


def test_moveannotations_samplesnotexists(testdir, mock_testclass):
    samples = 'samples.txt'
    distance = 2
    ma.moveannotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ma.moveannotations, ['-s', samples, '-d', distance])
    assert result.exit_code != 0
    ma.moveannotations_samples.assert_not_called()


def test_moveannotations_samples(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    distance = 2
    Parser.first = MagicMock(return_value=samples)
    ma.moveannotations_sample = MagicMock()
    ma.moveannotations_samples(samples_file, distance)
    Parser.first.assert_called_once_with(samples_file)
    for sample in samples:
        ma.moveannotations_sample.assert_any_call(sample, distance, True, True)


def test_moveannotations_samples_second(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    distance = 2
    Parser.first = MagicMock(return_value=samples)
    ma.moveannotations_sample = MagicMock()
    ma.moveannotations_samples(samples_file, distance, index=1)
    Parser.first.assert_called_once_with(samples_file)
    ma.moveannotations_sample.assert_any_call(samples[1], distance, True, True)

    
def test_moveannotations_sample(testdir, mock_testclass):
    sample = 'POLR2A'
    distance = 32
    src_bed = Path(__file__).parent.joinpath('sample.bed')
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    copyfile(src_bed, bed)
    ma.moveannotations_sample(sample, distance)
    assert os.path.exists(forcov)
    with open(forcov, 'r') as infile:
        assert 'track name=test\n' == infile.readline()
        assert 'chr1\t132\t182\ttest1\t1\t+\n' == infile.readline()
        assert 'chr2\t432\t482\ttest2\t2\t+\n' == infile.readline()
        assert 'chr3\t532\t682\ttest3\t3\t+\n' == infile.readline()
        assert 'chr4\t832\t782\ttest4\t4\t+\n' == infile.readline()
        assert 'chr5\t68\t118\ttest5\t1\t-\n' == infile.readline()
        assert 'chr6\t368\t418\ttest6\t2\t-\n' == infile.readline()
        assert 'chr7\t468\t618\ttest7\t3\t-\n' == infile.readline()
        assert 'chr8\t768\t718\ttest8\t4\t-\n' == infile.readline()


def test_moveannotations_sample_negativedistance(testdir, mock_testclass):
    sample = 'POLR2A'
    distance = -32
    src_bed = Path(__file__).parent.joinpath('sample.bed')
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    copyfile(src_bed, bed)
    ma.moveannotations_sample(sample, distance)
    assert os.path.exists(forcov)
    with open(forcov, 'r') as infile:
        assert 'track name=test\n' == infile.readline()
        assert 'chr1\t68\t118\ttest1\t1\t+\n' == infile.readline()
        assert 'chr2\t368\t418\ttest2\t2\t+\n' == infile.readline()
        assert 'chr3\t468\t618\ttest3\t3\t+\n' == infile.readline()
        assert 'chr4\t768\t718\ttest4\t4\t+\n' == infile.readline()
        assert 'chr5\t132\t182\ttest5\t1\t-\n' == infile.readline()
        assert 'chr6\t432\t482\ttest6\t2\t-\n' == infile.readline()
        assert 'chr7\t532\t682\ttest7\t3\t-\n' == infile.readline()
        assert 'chr8\t832\t782\ttest8\t4\t-\n' == infile.readline()


def test_moveannotations_sample_discardnegatives(testdir, mock_testclass):
    sample = 'POLR2A'
    distance = 101
    src_bed = Path(__file__).parent.joinpath('sample.bed')
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    copyfile(src_bed, bed)
    ma.moveannotations_sample(sample, distance, True, True)
    assert os.path.exists(forcov)
    with open(forcov, 'r') as infile:
        assert 'track name=test\n' == infile.readline()
        assert 'chr1\t201\t251\ttest1\t1\t+\n' == infile.readline()
        assert 'chr2\t501\t551\ttest2\t2\t+\n' == infile.readline()
        assert 'chr3\t601\t751\ttest3\t3\t+\n' == infile.readline()
        assert 'chr4\t901\t851\ttest4\t4\t+\n' == infile.readline()
        assert 'chr6\t299\t349\ttest6\t2\t-\n' == infile.readline()
        assert 'chr7\t399\t549\ttest7\t3\t-\n' == infile.readline()
        assert 'chr8\t699\t649\ttest8\t4\t-\n' == infile.readline()


def test_moveannotations_sample_keepnegatives(testdir, mock_testclass):
    sample = 'POLR2A'
    distance = 101
    src_bed = Path(__file__).parent.joinpath('sample.bed')
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    copyfile(src_bed, bed)
    ma.moveannotations_sample(sample, distance, False, True)
    assert os.path.exists(forcov)
    with open(forcov, 'r') as infile:
        assert 'track name=test\n' == infile.readline()
        assert 'chr1\t201\t251\ttest1\t1\t+\n' == infile.readline()
        assert 'chr2\t501\t551\ttest2\t2\t+\n' == infile.readline()
        assert 'chr3\t601\t751\ttest3\t3\t+\n' == infile.readline()
        assert 'chr4\t901\t851\ttest4\t4\t+\n' == infile.readline()
        assert 'chr5\t-1\t49\ttest5\t1\t-\n' == infile.readline()
        assert 'chr6\t299\t349\ttest6\t2\t-\n' == infile.readline()
        assert 'chr7\t399\t549\ttest7\t3\t-\n' == infile.readline()
        assert 'chr8\t699\t649\ttest8\t4\t-\n' == infile.readline()


def test_moveannotations_sample_forwardnegatives(testdir, mock_testclass):
    sample = 'POLR2A'
    distance = 32
    src_bed = Path(__file__).parent.joinpath('sample.bed')
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    copyfile(src_bed, bed)
    ma.moveannotations_sample(sample, distance, True, False)
    assert os.path.exists(forcov)
    with open(forcov, 'r') as infile:
        assert 'track name=test\n' == infile.readline()
        assert 'chr1\t132\t182\ttest1\t1\t+\n' == infile.readline()
        assert 'chr2\t432\t482\ttest2\t2\t+\n' == infile.readline()
        assert 'chr3\t532\t682\ttest3\t3\t+\n' == infile.readline()
        assert 'chr4\t832\t782\ttest4\t4\t+\n' == infile.readline()
        assert 'chr5\t132\t182\ttest5\t1\t-\n' == infile.readline()
        assert 'chr6\t432\t482\ttest6\t2\t-\n' == infile.readline()
        assert 'chr7\t532\t682\ttest7\t3\t-\n' == infile.readline()
        assert 'chr8\t832\t782\ttest8\t4\t-\n' == infile.readline()


def test_moveannotations_sample_negativedistance_forwardnegatives(testdir, mock_testclass):
    sample = 'POLR2A'
    distance = -32
    src_bed = Path(__file__).parent.joinpath('sample.bed')
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    copyfile(src_bed, bed)
    ma.moveannotations_sample(sample, distance, True, False)
    assert os.path.exists(forcov)
    with open(forcov, 'r') as infile:
        assert 'track name=test\n' == infile.readline()
        assert 'chr1\t68\t118\ttest1\t1\t+\n' == infile.readline()
        assert 'chr2\t368\t418\ttest2\t2\t+\n' == infile.readline()
        assert 'chr3\t468\t618\ttest3\t3\t+\n' == infile.readline()
        assert 'chr4\t768\t718\ttest4\t4\t+\n' == infile.readline()
        assert 'chr5\t68\t118\ttest5\t1\t-\n' == infile.readline()
        assert 'chr6\t368\t418\ttest6\t2\t-\n' == infile.readline()
        assert 'chr7\t468\t618\ttest7\t3\t-\n' == infile.readline()
        assert 'chr8\t768\t718\ttest8\t4\t-\n' == infile.readline()

