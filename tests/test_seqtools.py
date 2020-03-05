import logging
from pathlib import Path
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

import pandas as pd
from seqtools import seqtools, Bam2Bed, Bowtie2, Bwa, Download, FilterBam, GenomeCoverage, Intersect, Merge, MergeBigwigs, Plot2do, SlowSplit, Split, Statistics, Vap


@pytest.fixture
def mock_testclass():
    bam2bed_sample = Bam2Bed.bam2bed_sample
    bowtie_sample = Bowtie2.bowtie_sample
    bwa_sample = Bwa.bwa_sample
    download_sample = Download.download_sample
    filterbam_sample = FilterBam.filterbam_sample
    sample_splits_genome_coverage = GenomeCoverage.sample_splits_genome_coverage
    annotations_length = Intersect.annotations_length
    intersect_sample = Intersect.intersect_sample
    merge_samples = Merge.merge_samples
    merge_samples_bw = MergeBigwigs.merge_samples
    plot2do_sample = Plot2do.plot2do_sample
    slow_split_sample = SlowSplit.split_sample
    split_sample = Split.split_sample
    all_statistics = Statistics.all_statistics
    vap_sample = Vap.vap_sample
    read_csv = pd.read_csv
    yield
    Bam2Bed.bam2bed_sample = bam2bed_sample
    Bowtie2.bowtie_sample = bowtie_sample
    Bwa.bwa_sample = bwa_sample
    Download.download_sample = download_sample
    FilterBam.filterbam_sample = filterbam_sample
    GenomeCoverage.sample_splits_genome_coverage = sample_splits_genome_coverage
    Intersect.annotations_length = annotations_length
    Intersect.intersect_sample = intersect_sample
    Merge.merge_samples = merge_samples
    MergeBigwigs.merge_samples = merge_samples_bw
    Plot2do.plot2do_sample = plot2do_sample
    SlowSplit.split_sample = slow_split_sample
    Split.split_sample = split_sample
    Statistics.all_statistics = all_statistics
    Vap.vap_sample = vap_sample
    pd.read_csv = read_csv 


def test_seqtools_bam2bed(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    Bam2Bed.bam2bed_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bam2bed', '--samples', samples, '--unpaired', '--threads', threads, '--index', index])
    assert result.exit_code == 0
    Bam2Bed.bam2bed_sample.assert_called_once_with('POLR1C', False, threads)


def test_seqtools_bowtie2(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    Bowtie2.bowtie_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bowtie2', '--samples', samples, '--threads', threads, '--index', index])
    assert result.exit_code == 0
    Bowtie2.bowtie_sample.assert_called_once_with('POLR1C', threads, ())


def test_seqtools_bwa(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    fasta = Path(__file__).parent.joinpath('sacCer3.fa')
    threads = 3
    index = 2
    Bwa.bwa_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bwa', '--samples', samples, '--fasta', fasta, '--threads', threads, '--index', index])
    assert result.exit_code == 0
    Bwa.bwa_sample.assert_called_once_with('POLR1C', fasta, threads, ())


def test_seqtools_download(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    fasta = Path(__file__).parent.joinpath('sacCer3.fa')
    mem = '200MB'
    threads = 3
    index = 2
    Download.download_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['download', '--samples', samples, '--slow', '--threads', threads, '--mem', mem, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Download.download_sample.assert_called_once_with('POLR1C', 'SRR8518915', False, threads, mem)


def test_seqtools_filterbam(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    FilterBam.filterbam_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['filterbam', '--samples', samples, '--unpaired', '--threads', threads, '--index', index])
    assert result.exit_code == 0
    FilterBam.filterbam_sample.assert_called_once_with('POLR1C', False, threads)


def test_seqtools_genomecov(testdir, mock_testclass):
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


def test_seqtools_intersect(testdir, mock_testclass):
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


def test_seqtools_merge(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('merge.txt')
    index = 2
    Merge.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['merge', '--merge', samples, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Merge.merge_samples.assert_called_once_with('POLR1C', ['POLR1C_1', 'POLR1C_2'])


def test_seqtools_mergebw(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('merge.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    index = 2
    MergeBigwigs.merge_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['mergebw', '--merge', samples, '--sizes', sizes, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    MergeBigwigs.merge_samples.assert_called_once_with('POLR1C', ['POLR1C_1', 'POLR1C_2'], sizes)

 
def test_seqtools_plot2do(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    index = 2
    Plot2do.plot2do_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['plot2do', '--file', samples, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Plot2do.plot2do_sample.assert_called_once_with(Path(__file__).parent.joinpath('POLR1C'), None, None, None, None, None, None, None, None, None, None, None, None, None)


def test_seqtools_slowsplit(testdir, mock_testclass):
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


def test_seqtools_split(testdir, mock_testclass):
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


def test_seqtools_vap(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    parameters = Path(__file__).parent.joinpath('parameters.txt')
    index = 2
    Vap.vap_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['vap', '--samples', samples, '--parameters', parameters, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Vap.vap_sample.assert_called_once_with('POLR1C', parameters)
