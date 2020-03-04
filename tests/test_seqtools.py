import logging
from pathlib import Path
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

import pandas as pd
from seqtools import seqtools, Bam2Bed, Bowtie2, Bwa, DownloadSample, FilterBam, GenomeCoverage, Intersect, Merge, MergeBigwigs, Plot2do, SlowSplit, Split, Statistics, Vap


@pytest.fixture
def mock_testclass():
    read_csv = pd.read_csv
    yield
    pd.read_csv = read_csv 


def test_seqtools_bam2bed(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    Bam2Bed.bam2bed_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bam2bed', '--samples', samples, '--unpaired', '--threads', threads, '--index', index])
    assert result.exit_code == 0
    Bam2Bed.bam2bed_sample.assert_called_once_with('POLR1C', False, threads)


def test_seqtools_bowtie2(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    Bowtie2.bowtie_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bowtie2', '--samples', samples, '--threads', threads, '--index', index])
    assert result.exit_code == 0
    Bowtie2.bowtie_sample.assert_called_once_with('POLR1C', threads, ())


def test_seqtools_bwa(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    fasta = Path(__file__).parent.joinpath('sacCer3.fa')
    threads = 3
    index = 2
    Bwa.bwa_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bwa', '--samples', samples, '--fasta', fasta, '--threads', threads, '--index', index])
    assert result.exit_code == 0
    Bwa.bwa_sample.assert_called_once_with('POLR1C', fasta, threads, ())


def test_seqtools_download(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    fasta = Path(__file__).parent.joinpath('sacCer3.fa')
    mem = '200MB'
    threads = 3
    index = 2
    DownloadSample.download_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['download', '--samples', samples, '--slow', '--threads', threads, '--mem', mem, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    DownloadSample.download_sample.assert_called_once_with('POLR1C', 'SRR8518915', False, threads, mem)


def test_seqtools_filterbam(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    FilterBam.filterbam_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['filterbam', '--samples', samples, '--unpaired', '--threads', threads, '--index', index])
    assert result.exit_code == 0
    FilterBam.filterbam_sample.assert_called_once_with('POLR1C', False, threads)


def test_seqtools_genomecov(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '-'
    index = 2
    GenomeCoverage.sample_splits_genome_coverage = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['genomecov', '--samples', samples, '--sizes', sizes, '--scale', scale, '--strand', strand, '--index', index])
    assert result.exit_code == 0
    GenomeCoverage.sample_splits_genome_coverage.assert_called_once_with('POLR1C', sizes, False, False, scale, strand)


def test_seqtools_intersect(testdir):
    samples = Path(__file__).parent.joinpath('intersect.txt')
    annotations = Path(__file__).parent.joinpath('annotations.bed')
    annotations_length = 4
    index = 2
    Intersect.annotations_length = MagicMock(return_value=annotations_length)
    Intersect.intersect_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['intersect', '--samples', samples, '--annotations', annotations, '--index', index])
    assert result.exit_code == 0
    Intersect.intersect_sample.assert_called_once_with('POLR1C', 'POLR1C-inter', annotations, annotations_length)


def test_seqtools_merge(testdir):
    samples = Path(__file__).parent.joinpath('merge.txt')
    index = 2
    Merge.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['merge', '--merge', samples, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Merge.merge_samples.assert_called_once_with('POLR1C', ['POLR1C_1', 'POLR1C_2'])


def test_seqtools_mergebw(testdir):
    samples = Path(__file__).parent.joinpath('merge.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    index = 2
    MergeBigwigs.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['mergebw', '--merge', samples, '--sizes', sizes, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    MergeBigwigs.merge_samples.assert_called_once_with('POLR1C', ['POLR1C_1', 'POLR1C_2'], sizes)

 
def test_seqtools_plot2do(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    index = 2
    Plot2do.plot2do_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['plot2do', '--file', samples, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Plot2do.plot2do_sample.assert_called_once_with(Path(__file__).parent.joinpath('POLR1C'), None, None, None, None, None, None, None, None, None, None, None, None, None)


def test_seqtools_slowsplit(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    binlength = 20
    binminlength = 50
    binmaxlength = 150
    index = 2
    SlowSplit.split_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['slowsplit', '--samples', samples, '--binLength', binlength, '--binMinLength', binminlength, '--binMaxLength', binmaxlength, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    SlowSplit.split_sample.assert_called_once_with('POLR1C', binlength, binminlength, binmaxlength)


def test_seqtools_split(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    binlength = 20
    binminlength = 50
    binmaxlength = 150
    index = 2
    Split.split_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['split', '--samples', samples, '--binLength', binlength, '--binMinLength', binminlength, '--binMaxLength', binmaxlength, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Split.split_sample.assert_called_once_with('POLR1C', binlength, binminlength, binmaxlength)


def test_seqtools_statistics(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sample_pandas = ['samples-pandas']
    merge = Path(__file__).parent.joinpath('merge.txt')
    merge_pandas = ['merge-pandas']
    output = 'stats.txt'
    Statistics.all_statistics = MagicMock()
    pd.read_csv = MagicMock(side_effect=[sample_pandas, merge_pandas])
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['statistics', '--samples', samples, '--merge', merge, '--output', output])
    logging.warning(result.output)
    assert result.exit_code == 0
    Statistics.all_statistics.assert_called_once_with(sample_pandas[0], merge_pandas[0], output)


def test_seqtools_vap(testdir):
    samples = Path(__file__).parent.joinpath('samples.txt')
    parameters = Path(__file__).parent.joinpath('parameters.txt')
    index = 2
    Vap.vap_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['vap', '--samples', samples, '--parameters', parameters, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Vap.vap_sample.assert_called_once_with('POLR1C', parameters)
