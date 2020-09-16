import logging
from pathlib import Path
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import CenterAnnotations as ca
from seqtools import Split as sb
from seqtools.txt import Parser


@pytest.fixture
def mock_testclass():
    center_annotations_samples = ca.center_annotations_samples
    center_annotations_sample_splits = ca.center_annotations_sample_splits
    center_annotations_sample = ca.center_annotations_sample
    center_annotations = ca.center_annotations
    splits = sb.splits
    first = Parser.first
    yield
    ca.center_annotations_samples = center_annotations_samples
    ca.center_annotations_sample_splits = center_annotations_sample_splits
    ca.center_annotations_sample = center_annotations_sample
    ca.center_annotations = center_annotations
    sb.splits = splits
    Parser.first = first
    
    
def test_centerannotations(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    ca.center_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ca.centerannotations, ['-s', samples])
    assert result.exit_code == 0
    ca.center_annotations_samples.assert_called_once_with(samples, None)


def test_centerannotations_second(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    index = 1
    ca.center_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ca.centerannotations, ['-s', samples, '-i', index])
    assert result.exit_code == 0
    ca.center_annotations_samples.assert_called_once_with(samples, index)


def test_centerannotations_samplesnotexists(testdir, mock_testclass):
    samples = 'samples.txt'
    ca.center_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(ca.centerannotations, ['-s', samples])
    assert result.exit_code != 0
    ca.center_annotations_samples.assert_not_called()


def test_center_annotations_samples(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    Parser.first = MagicMock(return_value=samples)
    ca.center_annotations_sample_splits = MagicMock()
    ca.center_annotations_samples(samples_file)
    Parser.first.assert_called_once_with(samples_file)
    for sample in samples:
        ca.center_annotations_sample_splits.assert_any_call(sample)


def test_center_annotations_samples_second(testdir, mock_testclass):
    samples_file = Path(__file__).parent.joinpath('samples.txt')
    samples = ['POLR2A', 'ASDURF', 'POLR1C']
    Parser.first = MagicMock(return_value=samples)
    ca.center_annotations_sample_splits = MagicMock()
    ca.center_annotations_samples(samples_file, 1)
    Parser.first.assert_called_once_with(samples_file)
    ca.center_annotations_sample_splits.assert_called_once_with(samples[1])

    
def test_center_annotations_sample_splits(testdir, mock_testclass):
    sample = 'POLR2A'
    splits = ['POLR2A-100-110', 'POLR2A-120-130']
    ca.center_annotations_sample = MagicMock()
    sb.splits = MagicMock(return_value=splits)
    ca.center_annotations_sample_splits(sample)
    ca.center_annotations_sample.assert_any_call(sample)
    for split in splits:
        ca.center_annotations_sample.assert_any_call(split)

    
def test_center_annotations_sample_splits_notsplits(testdir, mock_testclass):
    sample = 'POLR2A'
    splits = []
    ca.center_annotations_sample = MagicMock()
    sb.splits = MagicMock(return_value=splits)
    ca.center_annotations_sample_splits(sample)
    ca.center_annotations_sample.assert_called_once_with(sample)

    
def test_center_annotations_sample(testdir, mock_testclass):
    sample = 'POLR2A'
    bed = sample + '.bed'
    forcov = sample + '-forcov.bed'
    ca.center_annotations = MagicMock()
    ca.center_annotations_sample(sample)
    ca.center_annotations.assert_called_once_with(bed, forcov)


def test_center_annotations_sample_split(testdir, mock_testclass):
    split = 'POLR2A-100-110'
    bed = split + '.bed'
    forcov = split + '-forcov.bed'
    ca.center_annotations = MagicMock()
    ca.center_annotations_sample(split)
    ca.center_annotations.assert_called_once_with(bed, forcov)

    
def test_center_annotations(testdir, mock_testclass):
    bed = Path(__file__).parent.joinpath('sample.bed')
    forcov = 'POLR2A-forcov.bed'
    ca.center_annotations(bed, forcov)
    with open(forcov, 'r') as infile:
        infile.readline() == 'chr1\t125\t126\ttest1\t1\t+\n'
        infile.readline() == 'chr2\t425\t426\ttest2\t2\t+\n'
        infile.readline() == 'chr3\t575\t576\ttest3\t3\t+\n'
        infile.readline() == 'chr4\t775\t775\ttest4\t4\t+\n'
        infile.readline() == 'chr5\t125\t126\ttest5\t1\t-\n'
        infile.readline() == 'chr6\t425\t426\ttest6\t2\t-\n'
        infile.readline() == 'chr7\t575\t576\ttest7\t3\t-\n'
        infile.readline() == 'chr8\t875\t776\ttest8\t4\t-\n'