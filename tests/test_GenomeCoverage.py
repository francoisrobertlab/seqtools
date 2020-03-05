import logging
import os
from pathlib import Path
from shutil import copyfile
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import GenomeCoverage as gc
from seqtools import Split as sb
from seqtools.bed import Bed

BASE_SCALE = 1000000


@pytest.fixture
def mock_testclass():
    coverage = gc.coverage
    genome_coverage = gc.genome_coverage
    sample_splits_genome_coverage = gc.sample_splits_genome_coverage
    splits = sb.splits
    sort = Bed.sort
    count_bed = Bed.count_bed
    bedgraph_to_bigwig = Bed.bedgraph_to_bigwig
    run = subprocess.run
    yield
    gc.coverage = coverage
    gc.genome_coverage = genome_coverage
    gc.sample_splits_genome_coverage = sample_splits_genome_coverage
    sb.splits = splits
    Bed.sort = sort
    Bed.count_bed = count_bed
    Bed.bedgraph_to_bigwig = bedgraph_to_bigwig
    subprocess.run = run
    
    
def create_file(*args, **kwargs):
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_genomecov(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = 'sacCer3.chrom.sizes'
    copyfile(Path(__file__).parent.joinpath('sizes.txt'), sizes)
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples ])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, False, False, None, None)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, False, False, None, None)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, False, False, None, None)


def test_genomecov_all_five(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--five' ])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, True, False, None, None)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, True, False, None, None)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, True, False, None, None)


def test_genomecov_second_five(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--five' , '-i', 1])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_called_once_with('ASDURF', sizes, True, False, None, None)


def test_genomecov_all_three(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--three' ])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, False, True, None, None)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, False, True, None, None)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, False, True, None, None)


def test_genomecov_second_three(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--three' , '-i', 1])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_called_once_with('ASDURF', sizes, False, True, None, None)


def test_genomecov_all_scale(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--scale', scale])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, False, False, scale, None)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, False, False, scale, None)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, False, False, scale, None)


def test_genomecov_second_scale(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--scale', scale, '-i', 1])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_called_once_with('ASDURF', sizes, False, False, scale, None)


def test_genomecov_all_scale_negativestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '-'
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--scale', scale, '--strand', strand])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, False, False, scale, strand)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, False, False, scale, strand)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, False, False, scale, strand)


def test_genomecov_second_scale_negativestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '-'
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--scale', scale, '--strand', strand, '-i', 1])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_called_once_with('ASDURF', sizes, False, False, scale, strand)


def test_genomecov_all_scale_positivestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '+'
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--scale', scale, '--strand', strand])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_any_call('POLR2A', sizes, False, False, scale, strand)
    gc.sample_splits_genome_coverage.assert_any_call('ASDURF', sizes, False, False, scale, strand)
    gc.sample_splits_genome_coverage.assert_any_call('POLR1C', sizes, False, False, scale, strand)


def test_genomecov_second_scale_positivestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '+'
    gc.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(gc.genomecov, ['-s', samples, '--sizes', sizes, '--scale', scale, '--strand', strand, '-i', 1])
    assert result.exit_code == 0
    gc.sample_splits_genome_coverage.assert_called_once_with('ASDURF', sizes, False, False, scale, strand)

    
def test_sample_splits_genome_coverage(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, False, False, None, None)
    gc.genome_coverage.assert_any_call(split1, sizes, False, False, None, None)
    gc.genome_coverage.assert_any_call(split2, sizes, False, False, None, None)


def test_sample_splits_genome_coverage_five(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, five=True)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, True, False, None, None)
    gc.genome_coverage.assert_any_call(split1, sizes, True, False, None, None)
    gc.genome_coverage.assert_any_call(split2, sizes, True, False, None, None)


def test_sample_splits_genome_coverage_three(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, three=True)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, False, True, None, None)
    gc.genome_coverage.assert_any_call(split1, sizes, False, True, None, None)
    gc.genome_coverage.assert_any_call(split2, sizes, False, True, None, None)


def test_sample_splits_genome_coverage_scale(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    scale = 1.5
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, False, False, scale)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, False, False, scale, None)
    gc.genome_coverage.assert_any_call(split1, sizes, False, False, scale, None)
    gc.genome_coverage.assert_any_call(split2, sizes, False, False, scale, None)


def test_sample_splits_genome_coverage_negativestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    scale = 1.5
    strand = '-'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, False, False, None, strand)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, False, False, None, strand)
    gc.genome_coverage.assert_any_call(split1, sizes, False, False, None, strand)
    gc.genome_coverage.assert_any_call(split2, sizes, False, False, None, strand)


def test_sample_splits_genome_coverage_scale_negativestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    scale = 1.5
    strand = '-'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, False, False, scale, strand)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, False, False, scale, strand)
    gc.genome_coverage.assert_any_call(split1, sizes, False, False, scale, strand)
    gc.genome_coverage.assert_any_call(split2, sizes, False, False, scale, strand)


def test_sample_splits_genome_coverage_positivestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    strand = '+'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, False, False, None, strand)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, False, False, None, strand)
    gc.genome_coverage.assert_any_call(split1, sizes, False, False, None, strand)
    gc.genome_coverage.assert_any_call(split2, sizes, False, False, None, strand)


def test_sample_splits_genome_coverage_scale_positivestrand(testdir, mock_testclass):
    sample = 'POLR2A'
    sizes = 'human.sizes'
    split1 = sample + '-100-110'
    split2 = sample + '-110-120'
    scale = 1.5
    strand = '+'
    sb.splits = MagicMock(return_value=[split1, split2])
    gc.genome_coverage = MagicMock()
    gc.sample_splits_genome_coverage(sample, sizes, False, False, scale, strand)
    sb.splits.assert_called_once_with(sample)
    gc.genome_coverage.assert_any_call(sample, sizes, False, False, scale, strand)
    gc.genome_coverage.assert_any_call(split1, sizes, False, False, scale, strand)
    gc.genome_coverage.assert_any_call(split2, sizes, False, False, scale, strand)

    
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
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, False, False, BASE_SCALE / count, None)
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
    gc.coverage.assert_called_once_with(bed, cov, sizes, sample, False, False, BASE_SCALE / count, None)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_genome_coverage_five(testdir, mock_testclass):
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
    gc.genome_coverage(sample, sizes, five=True)
    Bed.count_bed.assert_called_once_with(forcov)
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, True, False, BASE_SCALE / count, None)
    Bed.bedgraph_to_bigwig.assert_called_once_with(cov, bw, sizes)


def test_genome_coverage_three(testdir, mock_testclass):
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
    gc.genome_coverage(sample, sizes, three=True)
    Bed.count_bed.assert_called_once_with(forcov)
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, False, True, BASE_SCALE / count, None)
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
    gc.genome_coverage(sample, sizes, False, False, scale)
    Bed.count_bed.assert_not_called()
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, False, False, scale, None)
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
    gc.genome_coverage(sample, sizes, False, False, None, strand)
    Bed.count_bed.assert_called_once_with(forcov)
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, False, False, BASE_SCALE / count, strand)
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
    gc.genome_coverage(sample, sizes, False, False, scale, strand)
    Bed.count_bed.assert_not_called()
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, False, False, scale, strand)
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
    gc.genome_coverage(sample, sizes, False, False, None, strand)
    Bed.count_bed.assert_called_once_with(forcov)
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, False, False, BASE_SCALE / count, strand)
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
    gc.genome_coverage(sample, sizes, False, False, scale, strand)
    Bed.count_bed.assert_not_called()
    gc.coverage.assert_called_once_with(forcov, cov, sizes, sample, False, False, scale, strand)
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
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-i', bed, '-g', sizes], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + '"\n'


def test_coverage_five(testdir):
    sample = 'POLR2A'
    bed = sample + '.bed'
    coverage_output = bed + '.cov'
    sort_output = bed + '.sort'
    output = sample + '-out.bed'
    sizes = 'human.sizes'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort = MagicMock(side_effect=create_file(['-o', sort_output]))
    gc.coverage(bed, output, sizes, sample, five=True)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-i', bed, '-g', sizes, '-5'], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + '"\n'


def test_coverage_three(testdir):
    sample = 'POLR2A'
    bed = sample + '.bed'
    coverage_output = bed + '.cov'
    sort_output = bed + '.sort'
    output = sample + '-out.bed'
    sizes = 'human.sizes'
    subprocess.run = MagicMock(side_effect=create_file)
    Bed.sort = MagicMock(side_effect=create_file(['-o', sort_output]))
    gc.coverage(bed, output, sizes, sample, three=True)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-i', bed, '-g', sizes, '-3'], stdout=ANY, check=True)
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
    gc.coverage(bed, output, sizes, sample, False, False, scale)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-i', bed, '-g', sizes, '-scale', str(scale)], stdout=ANY, check=True)
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
    gc.coverage(bed, output, sizes, sample, False, False, None, strand)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-i', bed, '-g', sizes, '-strand', strand], stdout=ANY, check=True)
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
    gc.coverage(bed, output, sizes, sample, False, False, scale, strand)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-i', bed, '-g', sizes, '-scale', str(scale), '-strand', strand], stdout=ANY, check=True)
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
    gc.coverage(bed, output, sizes, sample, False, False, None, strand)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-i', bed, '-g', sizes, '-strand', strand], stdout=ANY, check=True)
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
    gc.coverage(bed, output, sizes, sample, False, False, scale, strand)
    subprocess.run.assert_called_once_with(['bedtools', 'genomecov', '-bg', '-i', bed, '-g', sizes, '-scale', str(scale), '-strand', strand], stdout=ANY, check=True)
    Bed.sort.assert_called_once_with(coverage_output, sort_output)
    assert not os.path.exists(coverage_output)
    assert not os.path.exists(sort_output)
    with open(output, 'r') as infile:
        assert infile.readline() == 'track type=bedGraph name="' + sample + ' Plus"\n'
