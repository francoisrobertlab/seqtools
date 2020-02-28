import logging
import os
from pathlib import Path
from shutil import copyfile
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
from more_itertools.more import side_effect
import pytest

from seqtools import GenomeCoverage as gc
from seqtools import SplitBed as sb
from seqtools.bed import Bed

BASE_SCALE = 1000000


@pytest.fixture
def mock_testclass():
    coverage = gc.coverage
    genome_coverage = gc.genome_coverage
    sample_splits_genome_coverage = gc.sample_splits_genome_coverage
    yield coverage, genome_coverage, sample_splits_genome_coverage
    gc.coverage = coverage
    gc.genome_coverage = genome_coverage
    gc.sample_splits_genome_coverage = sample_splits_genome_coverage
    
    
def create_file(*args, **kwargs):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_main(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = 'sacCer3.chrom.sizes'
    copyfile(Path(__file__).parent.joinpath('sizes.txt'), sizes)
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.main, ['-s', samples ])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, None, None)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, None, None)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, None, None)


def test_main_all_scale(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.main, ['-s', samples, '--sizes', sizes, '--scale', scale])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, scale, None)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, scale, None)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, scale, None)


def test_main_second_scale(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.main, ['-s', samples, '--sizes', sizes, '--scale', scale, '-i', 1])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_called_once_with('ASDURF', sizes, scale, None)


def test_main_all_scale_negativestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '-'
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.main, ['-s', samples, '--sizes', sizes, '--scale', scale, '--strand', strand])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, scale, strand)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, scale, strand)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, scale, strand)


def test_main_second_scale_negativestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '-'
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.main, ['-s', samples, '--sizes', sizes, '--scale', scale, '--strand', strand, '-i', 1])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_called_once_with('ASDURF', sizes, scale, strand)


def test_main_all_scale_positivestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '+'
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.main, ['-s', samples, '--sizes', sizes, '--scale', scale, '--strand', strand])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, scale, strand)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, scale, strand)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, scale, strand)


def test_main_second_scale_positivestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '+'
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.main, ['-s', samples, '--sizes', sizes, '--scale', scale, '--strand', strand, '-i', 1])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_called_once_with('ASDURF', sizes, scale, strand)

    
def test_sample_splits_genome_coverage(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, None, None)
    gc.genome_coverage.assert_any_call(split1, sizes, None, None)
    gc.genome_coverage.assert_any_call(split2, sizes, None, None)


def test_sample_splits_genome_coverage_scale(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    scale = 1.5
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, scale)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, scale, None)
    gc.genome_coverage.assert_any_call(split1, sizes, scale, None)
    gc.genome_coverage.assert_any_call(split2, sizes, scale, None)


def test_sample_splits_genome_coverage_negativestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    scale = 1.5
    strand = '-'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, None, strand)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, None, strand)
    gc.genome_coverage.assert_any_call(split1, sizes, None, strand)
    gc.genome_coverage.assert_any_call(split2, sizes, None, strand)


def test_sample_splits_genome_coverage_scale_negativestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    scale = 1.5
    strand = '-'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, scale, strand)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, scale, strand)
    gc.genome_coverage.assert_any_call(split1, sizes, scale, strand)
    gc.genome_coverage.assert_any_call(split2, sizes, scale, strand)


def test_sample_splits_genome_coverage_positivestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    strand = '+'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, None, strand)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, None, strand)
    gc.genome_coverage.assert_any_call(split1, sizes, None, strand)
    gc.genome_coverage.assert_any_call(split2, sizes, None, strand)


def test_sample_splits_genome_coverage_scale_positivestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    scale = 1.5
    strand = '+'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, scale, strand)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, scale, strand)
    gc.genome_coverage.assert_any_call(split1, sizes, scale, strand)
    gc.genome_coverage.assert_any_call(split2, sizes, scale, strand)

    
def test_genome_coverage(testdir, mock_testclass):
    sample = 'POLR2A'
    forcov = sample + '-forcov.bed'
    with open(forcov, 'w') as outfile:
        outfile.write('test')
    cov = sample + '-cov.bed'
    bw = sample + '-cov.bw'
    sizes = 'human.sizes'
    count = 2000000
    Bed.count_bed = MagicMock(return_value=count)
    gc.coverage = MagicMock()
    Bed.bedgraph_to_bigwig = MagicMock()
    gc.genome_coverage(sample, sizes)
    Bed.count_bed.assert_called_once_with(forcov)
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, BASE_SCALE / count, None)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_genome_coverage_noforcov(testdir, mock_testclass):
    sample = 'POLR2A'
    bed = sample + '.bed'
    cov = sample + '-cov.bed'
    bw = sample + '-cov.bw'
    sizes = 'human.sizes'
    count = 2000000
    Bed.count_bed = MagicMock(return_value=count)
    gc.coverage = MagicMock()
    Bed.bedgraph_to_bigwig = MagicMock()
    gc.genome_coverage(sample, sizes)
    Bed.count_bed.assert_called_once_with(bed)
    gc.coverage.assert_called_once_with(bed, cov, sizes, sample, BASE_SCALE / count, None)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_genome_coverage_scale(testdir, mock_testclass):
    sample = 'POLR2A'
    forcov = sample + '-forcov.bed'
    with open(forcov, 'w') as outfile:
        outfile.write('test')
    cov = sample + '-cov.bed'
    bw = sample + '-cov.bw'
    sizes = 'human.sizes'
    scale = 1.5
    Bed.count_bed = MagicMock()
    gc.coverage = MagicMock()
    Bed.bedgraph_to_bigwig = MagicMock()
    gc.genome_coverage(sample, sizes, scale)
    Bed.count_bed.assert_not_called()
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, scale, None)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_genome_coverage_negativestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    forcov = sample + '-forcov.bed'
    with open(forcov, 'w') as outfile:
        outfile.write('test')
    cov = sample + '-cov-neg.bed'
    bw = sample + '-cov-neg.bw'
    sizes = 'human.sizes'
    strand = '-'
    count = 2000000
    Bed.count_bed = MagicMock(return_value=count)
    gc.coverage = MagicMock()
    Bed.bedgraph_to_bigwig = MagicMock()
    gc.genome_coverage(sample, sizes, None, strand)
    Bed.count_bed.assert_called_once_with(forcov)
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, BASE_SCALE / count, strand)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_genome_coverage_scale_negativestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    forcov = sample + '-forcov.bed'
    with open(forcov, 'w') as outfile:
        outfile.write('test')
    cov = sample + '-cov-neg.bed'
    bw = sample + '-cov-neg.bw'
    sizes = 'human.sizes'
    scale = 1.5
    strand = '-'
    Bed.count_bed = MagicMock()
    gc.coverage = MagicMock()
    Bed.bedgraph_to_bigwig = MagicMock()
    gc.genome_coverage(sample, sizes, scale, strand)
    Bed.count_bed.assert_not_called()
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, scale, strand)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_genome_coverage_positivestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    forcov = sample + '-forcov.bed'
    with open(forcov, 'w') as outfile:
        outfile.write('test')
    cov = sample + '-cov-pos.bed'
    bw = sample + '-cov-pos.bw'
    sizes = 'human.sizes'
    strand = '+'
    count = 2000000
    Bed.count_bed = MagicMock(return_value=count)
    gc.coverage = MagicMock()
    Bed.bedgraph_to_bigwig = MagicMock()
    gc.genome_coverage(sample, sizes, None, strand)
    Bed.count_bed.assert_called_once_with(forcov)
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, BASE_SCALE / count, strand)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_genome_coverage_scale_positivestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    forcov = sample + '-forcov.bed'
    with open(forcov, 'w') as outfile:
        outfile.write('test')
    cov = sample + '-cov-pos.bed'
    bw = sample + '-cov-pos.bw'
    sizes = 'human.sizes'
    scale = 1.5
    strand = '+'
    Bed.count_bed = MagicMock()
    gc.coverage = MagicMock()
    Bed.bedgraph_to_bigwig = MagicMock()
    gc.genome_coverage(sample, sizes, scale, strand)
    Bed.count_bed.assert_not_called()
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, scale, strand)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_coverage(testdir):
    sample = 'POLR2A'
    bed = sample + '.bed'
    coverage_output = bed + '.cov'
    sort_output = bed + '.sort'
    output = sample + '-out.bed'
    sizes = 'human.sizes'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort = MagicMock(side_effect=create_file(['-o', sort_output]))
    gc.coverage(bed, output, sizes, sample)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-5', '-i', bed, '-g', sizes], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + '"\n'


def test_coverage_scale(testdir):
    sample = 'POLR2A'
    bed = sample + '.bed'
    coverage_output = bed + '.cov'
    sort_output = bed + '.sort'
    output = sample + '-out.bed'
    sizes = 'human.sizes'
    scale = 1.5
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort = MagicMock(side_effect=create_file(['-o', sort_output]))
    gc.coverage(bed, output, sizes, sample, scale)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-5', '-i', bed, '-g', sizes, '-scale', str(scale)], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + '"\n'


def test_coverage_negativestrand(testdir):
    sample = 'POLR2A'
    bed = sample + '.bed'
    coverage_output = bed + '.cov'
    sort_output = bed + '.sort'
    output = sample + '-out.bed'
    sizes = 'human.sizes'
    strand = '-'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort = MagicMock(side_effect=create_file(['-o', sort_output]))
    gc.coverage(bed, output, sizes, sample, None, strand)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-5', '-i', bed, '-g', sizes, '-strand', strand], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + ' Minus"\n'


def test_coverage_scale_negativestrand(testdir):
    sample = 'POLR2A'
    bed = sample + '.bed'
    coverage_output = bed + '.cov'
    sort_output = bed + '.sort'
    output = sample + '-out.bed'
    sizes = 'human.sizes'
    scale = 1.5
    strand = '-'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort = MagicMock(side_effect=create_file(['-o', sort_output]))
    gc.coverage(bed, output, sizes, sample, scale, strand)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-5', '-i', bed, '-g', sizes, '-scale', str(scale), '-strand', strand], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + ' Minus"\n'


def test_coverage_positivestrand(testdir):
    sample = 'POLR2A'
    bed = sample + '.bed'
    coverage_output = bed + '.cov'
    sort_output = bed + '.sort'
    output = sample + '-out.bed'
    sizes = 'human.sizes'
    strand = '+'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort = MagicMock(side_effect=create_file(['-o', sort_output]))
    gc.coverage(bed, output, sizes, sample, None, strand)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-5', '-i', bed, '-g', sizes, '-strand', strand], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + ' Plus"\n'


def test_coverage_scale_positivestrand(testdir):
    sample = 'POLR2A'
    bed = sample + '.bed'
    coverage_output = bed + '.cov'
    sort_output = bed + '.sort'
    output = sample + '-out.bed'
    sizes = 'human.sizes'
    scale = 1.5
    strand = '+'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort = MagicMock(side_effect=create_file(['-o', sort_output]))
    gc.coverage(bed, output, sizes, sample, scale, strand)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-5', '-i', bed, '-g', sizes, '-scale', str(scale), '-strand', strand], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + ' Plus"\n'
