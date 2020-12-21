import logging
from pathlib import Path
import pytest
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner

from seqtools import seqtools, Bam2Bed, Bowtie2, Bwa, CenterAnnotations, ChipexoQual, Download, FilterBam, Fixmd5, GenomeCoverage, IgnoreStrand, Intersect, Merge, MergeBam, MergeBigwigs, Plot2do, RemoveSecondMate, Rename, ShiftAnnotations, SlowSplit, Split, Statistics, Vap


@pytest.fixture
def mock_testclass():
    bam2bed_samples = Bam2Bed.bam2bed_samples
    bowtie_samples = Bowtie2.bowtie_samples
    bwa_samples = Bwa.bwa_samples
    center_annotations_samples = CenterAnnotations.center_annotations_samples
    chipexoqual_datasets = ChipexoQual.chipexoqual_datasets
    download_samples = Download.download_samples
    filter_bam = FilterBam.filter_bam
    fixmd5 = Fixmd5.fixmd5_
    genome_coverage_samples = GenomeCoverage.genome_coverage_samples
    ignore_strand_samples = IgnoreStrand.ignore_strand_samples
    intersect_samples = Intersect.intersect_samples
    merge_datasets = Merge.merge_datasets
    merge_datasets_bam = MergeBam.merge_datasets
    merge_datasets_bw = MergeBigwigs.merge_datasets
    plot2do_samples = Plot2do.plot2do_samples
    removesecondmate_samples = RemoveSecondMate.removesecondmate_samples
    rename = Rename.rename_
    shift_annotations_samples = ShiftAnnotations.shift_annotations_samples
    slow_split_samples = SlowSplit.split_samples
    split_samples = Split.split_samples
    statistics_samples = Statistics.statistics_samples
    vap_samples = Vap.vap_samples
    yield
    Bam2Bed.bam2bed_samples = bam2bed_samples
    Bowtie2.bowtie_samples = bowtie_samples
    Bwa.bwa_samples = bwa_samples
    CenterAnnotations.center_annotations_samples = center_annotations_samples
    ChipexoQual.chipexoqual_datasets = chipexoqual_datasets
    Download.download_samples = download_samples
    FilterBam.filter_bam = filter_bam
    Fixmd5.fixmd5_ = fixmd5
    GenomeCoverage.genome_coverage_samples = genome_coverage_samples
    IgnoreStrand.ignore_strand_samples = ignore_strand_samples
    Intersect.intersect_samples = intersect_samples
    Merge.merge_datasets = merge_datasets
    MergeBam.merge_datasets = merge_datasets_bam
    MergeBigwigs.merge_datasets = merge_datasets_bw
    Plot2do.plot2do_samples = plot2do_samples
    RemoveSecondMate.removesecondmate_samples = removesecondmate_samples
    Rename.rename_ = rename
    ShiftAnnotations.shift_annotations_samples = shift_annotations_samples
    SlowSplit.split_samples = slow_split_samples
    Split.split_samples = split_samples
    Statistics.statistics_samples = statistics_samples
    Vap.vap_samples = vap_samples


def test_seqtools_bam2bed(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    Bam2Bed.bam2bed_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bam2bed', '--samples', samples, '--unpaired', '--threads', threads, '--index', index])
    assert result.exit_code == 0
    Bam2Bed.bam2bed_samples.assert_called_once_with(samples, False, threads, '-dedup', index)


def test_seqtools_bowtie2(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    Bowtie2.bowtie_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bowtie2', '--samples', samples, '--threads', threads, '--index', index, '-x', 'sacCer3.fa'])
    assert result.exit_code == 0
    Bowtie2.bowtie_samples.assert_called_once_with(samples, threads, '', index, ('-x', 'sacCer3.fa'))


def test_seqtools_bwa(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    fasta = Path(__file__).parent.joinpath('sacCer3.fa')
    threads = 3
    index = 2
    Bwa.bwa_sample = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['bwa', '--samples', samples, '--fasta', fasta, '--threads', threads, '--index', index])
    assert result.exit_code == 0
    Bwa.bwa_sample.assert_called_once_with('POLR1C', fasta, threads, '', ())


def test_seqtools_centerannotations(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    CenterAnnotations.center_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['centerannotations', '--samples', samples])
    assert result.exit_code == 0
    CenterAnnotations.center_annotations_samples.assert_called_once_with(samples, '', '-forcov', None)


def test_seqtools_chipexoqual(testdir, mock_testclass):
    datasets = Path(__file__).parent.joinpath('dataset.txt')
    ChipexoQual.chipexoqual_datasets = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['chipexoqual', '--datasets', datasets])
    assert result.exit_code == 0
    ChipexoQual.chipexoqual_datasets.assert_called_once_with(datasets, '', None, ())


def test_seqtools_download(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    fasta = Path(__file__).parent.joinpath('sacCer3.fa')
    mem = '200MB'
    threads = 3
    index = 2
    Download.download_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['download', '--samples', samples, '--slow', '--threads', threads, '--mem', mem, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Download.download_samples.assert_called_once_with(samples, False, threads, mem, index)


def test_seqtools_filterbam(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 3
    index = 2
    FilterBam.filter_bam = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['filterbam', '--samples', samples, '--unpaired', '--threads', threads, '--index', index])
    assert result.exit_code == 0
    FilterBam.filter_bam.assert_called_once_with(samples, False, True, threads, '', '', index)


def test_seqtools_fixmd5(testdir, mock_testclass):
    Fixmd5.fixmd5_ = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['fixmd5'])
    assert result.exit_code == 0
    Fixmd5.fixmd5_.assert_called_once_with('*.md5', False)


def test_seqtools_genomecov(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    scale = 1.5
    strand = '-'
    input_suffix = '-forcov'
    output_suffix = '-cov'
    index = 2
    GenomeCoverage.genome_coverage_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['genomecov', '--samples', samples, '-g', sizes, '-5', '-scale', scale, '-strand', strand, '-is', input_suffix, '-os', output_suffix, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    GenomeCoverage.genome_coverage_samples.assert_called_once_with(samples, sizes, scale, strand, input_suffix, output_suffix, None, None, index, ('-5',))


def test_seqtools_ignorestrand(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    IgnoreStrand.ignore_strand_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['ignorestrand', '--samples', samples])
    assert result.exit_code == 0
    IgnoreStrand.ignore_strand_samples.assert_called_once_with(samples, '', '-forcov', None)


def test_seqtools_intersect(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('intersect.txt')
    annotations = Path(__file__).parent.joinpath('annotations.bed')
    index = 2
    Intersect.intersect_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['intersect', '--samples', samples, '--annotations', annotations, '--index', index])
    assert result.exit_code == 0
    Intersect.intersect_samples.assert_called_once_with(samples, annotations, index)


def test_seqtools_merge(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('dataset.txt')
    index = 2
    Merge.merge_datasets = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['merge', '--datasets', samples, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Merge.merge_datasets.assert_called_once_with(samples, index)


def test_seqtools_mergebam(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('dataset.txt')
    index = 2
    MergeBam.merge_datasets = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['mergebam', '--datasets', samples])
    logging.warning(result.output)
    assert result.exit_code == 0
    MergeBam.merge_datasets.assert_called_once_with(samples, '', 1, None)


def test_seqtools_mergebw(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('dataset.txt')
    sizes = Path(__file__).parent.joinpath('sizes.txt')
    index = 2
    MergeBigwigs.merge_datasets = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['mergebw', '--datasets', samples, '--sizes', sizes, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    MergeBigwigs.merge_datasets.assert_called_once_with(samples, sizes, index)

 
def test_seqtools_plot2do(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    index = 2
    type = 'dyads'
    Plot2do.plot2do_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['plot2do', '--file', samples, '--type', type, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Plot2do.plot2do_samples.assert_called_once_with(samples, '', index, ('--type', type,))


def test_seqtools_removesecondmate(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    RemoveSecondMate.removesecondmate_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['removesecondmate', '--samples', samples])
    assert result.exit_code == 0
    RemoveSecondMate.removesecondmate_samples.assert_called_once_with(samples, '-dedup', '-mate1', 1, None)


def test_seqtools_rename(testdir, mock_testclass):
    names = Path(__file__).parent.joinpath('names.txt')
    Rename.rename_ = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['rename', '--names', names])
    logging.warning(result.output)
    assert result.exit_code == 0
    Rename.rename_.assert_called_once_with(names, True, False, False)


def test_seqtools_shiftannotations(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    ShiftAnnotations.shift_annotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['shiftannotations', '--samples', samples])
    assert result.exit_code == 0
    ShiftAnnotations.shift_annotations_samples.assert_called_once_with(samples, '', '-forcov', None, ())


def test_seqtools_slowsplit(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    binlength = 20
    binminlength = 50
    binmaxlength = 150
    index = 2
    SlowSplit.split_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['slowsplit', '--samples', samples, '--binLength', binlength, '--binMinLength', binminlength, '--binMaxLength', binmaxlength, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    SlowSplit.split_samples.assert_called_once_with(samples, index, binlength, binminlength, binmaxlength)


def test_seqtools_split(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    binlength = 20
    binminlength = 50
    binmaxlength = 150
    index = 2
    Split.split_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['split', '--samples', samples, '--binLength', binlength, '--binMinLength', binminlength, '--binMaxLength', binmaxlength, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Split.split_samples.assert_called_once_with(samples, index, binlength, binminlength, binmaxlength)


def test_seqtools_statistics(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    datasets = Path(__file__).parent.joinpath('dataset.txt')
    output = 'stats.txt'
    Statistics.statistics_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['statistics', '--samples', samples, '--datasets', datasets, '--output', output])
    logging.warning(result.output)
    assert result.exit_code == 0
    Statistics.statistics_samples.assert_called_once_with(samples, datasets, output)


def test_seqtools_vap(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    parameters = Path(__file__).parent.joinpath('parameters.txt')
    index = 2
    Vap.vap_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(seqtools.seqtools, ['vap', '--samples', samples, '--parameters', parameters, '--index', index])
    logging.warning(result.output)
    assert result.exit_code == 0
    Vap.vap_samples.assert_called_once_with(samples, parameters, None, index)
