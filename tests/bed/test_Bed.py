import logging
import os
from pathlib import Path
from shutil import copyfile
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools.bed import Bed


@pytest.fixture
def mock_testclass():
    run = subprocess.run
    yield 
    subprocess.run = run


def create_file(*args, **kwargs):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_count_bed(testdir):
    bed = Path(__file__).parent.parent.joinpath('sample.bed')
    assert 8 == Bed.count_bed(bed)


def test_count_bed_strand(testdir):
    bed = Path(__file__).parent.parent.joinpath('sample.bed')
    assert 4 == Bed.count_bed(bed, '+')


def test_empty_bed(testdir):
    bed = 'sample.bed'
    sample = 'POLR2A'
    Bed.empty_bed(bed, sample)
    assert os.path.exists(bed)
    with open(bed, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + '"\n'
        assert infile.readline() == ''


def test_empty_bed_plusstrand(testdir):
    bed = 'sample.bed'
    sample = 'POLR2A'
    strand = '+'
    Bed.empty_bed(bed, sample, strand)
    assert os.path.exists(bed)
    with open(bed, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + ' Plus"\n'
        assert infile.readline() == ''


def test_empty_bed_minusstrand(testdir):
    bed = 'sample.bed'
    sample = 'POLR2A'
    strand = '-'
    Bed.empty_bed(bed, sample, strand)
    assert os.path.exists(bed)
    with open(bed, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + ' Minus"\n'
        assert infile.readline() == ''


def test_sort(testdir):
    bed = Path(__file__).parent.parent.joinpath('sample.bed')
    output = 'test.bed'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort(bed, output)
    subprocess.run.assert_called_with(['bedtools', 'sort', '-i', bed], stdout=ANY, check=True)
    assert os.path.exists(output)


def test_sort_bysize(testdir):
    bed = Path(__file__).parent.parent.joinpath('sample.bed')
    output = 'test.bed'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort_bysize(bed, output)
    subprocess.run.assert_called_with(['bedtools', 'sort', '-sizeA', '-i', bed], stdout=ANY, check=True)
    assert os.path.exists(output)


def test_bedgraph_to_bigwig(testdir):
    bed = Path(__file__).parent.parent.joinpath('sample.bed')
    sizes = Path(__file__).parent.parent.joinpath('sizes.txt')
    output = 'test.bw'
    subprocess.run = MagicMock(side_effect=create_file(['-o', output]))
    Bed.bedgraph_to_bigwig(bed, output, sizes)
    subprocess.run.assert_called_with(['bedGraphToBigWig', bed, sizes, output], check=True)
    assert os.path.exists(output)
