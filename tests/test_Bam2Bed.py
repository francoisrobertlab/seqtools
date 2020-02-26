import logging
import os
from pathlib import Path
from shutil import copyfile
import subprocess
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
import pytest

from seqtools import Bam2Bed as bb
from seqtools.bed import Bed


@pytest.fixture
def mock_testclass():
    bam2bed = bb.bam2bed
    bam2bedpe = bb.bam2bedpe
    bedpe2bed = bb.bedpe2bed
    yield bam2bed, bam2bedpe, bedpe2bed
    bb.bam2bed = bam2bed
    bb.bam2bedpe = bam2bedpe
    bb.bedpe2bed = bedpe2bed
    
    
def create_file(*args, **kwargs):
    logging.warning('create_file with: {}'.format(args))
    logging.warning('create_file with kw: {}'.format(kwargs))
    if 'stdout' in kwargs:
        outfile = kwargs['stdout']
        outfile.write('test')
    elif '-o' in args[0]:
        output = args[0][args[0].index('-o') + 1]
        with open(output, 'w') as outfile:
            outfile.write('test')


def test_main(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    threads = 1
    bb.bam2bed = MagicMock()
    runner = CliRunner()
    result = runner.invoke(bb.main, ['-s', samples ])
    assert result.exit_code == 0
    bb.bam2bed.assert_any_call('POLR2A', threads)
    bb.bam2bed.assert_any_call('ASDURF', threads)
    bb.bam2bed.assert_any_call('POLR1C', threads)

    
def test_bam2bed(testdir, mock_testclass):
    sample = 'POLR2A'
    bam = sample + '.bam'
    bedpe = sample + '.bedpe'
    bed = sample + '.bed'
    threads = 2
    bb.bam2bedpe = MagicMock(side_effect=create_file(['-o', bedpe]))
    bb.bedpe2bed = MagicMock()
    bb.bam2bed(sample, threads)
    bb.bam2bedpe.assert_called_with(bam, bedpe, threads)
    bb.bedpe2bed.assert_called_with(bedpe, bed)
    assert not os.path.exists(bedpe)
    
    
def test_bam2bedpe(testdir):
    bam = 'POLR2A.bam'
    bam_sort = bam + '.sort'
    bedpe = 'POLR2A.bedpe'
    threads = 2
    subprocess.run = MagicMock(side_effect=create_file)
    bb.bam2bedpe(bam, bedpe, threads)
    subprocess.run.assert_any_call(['samtools', 'sort', '-n', '--threads', str(threads - 1), '-o', bam_sort, bam], check=True)
    subprocess.run.assert_any_call(['bedtools', 'bamtobed', '-bedpe', '-mate1', '-i', bam_sort], stdout=ANY, check=True)
    assert not os.path.exists(bam_sort)
    assert os.path.exists(bedpe)


def test_bam2bedpe_singlethread(testdir):
    bam = 'POLR2A.bam'
    bam_sort = bam + '.sort'
    bedpe = 'POLR2A.bedpe'
    threads = 1
    subprocess.run = MagicMock(side_effect=create_file)
    bb.bam2bedpe(bam, bedpe, threads)
    subprocess.run.assert_any_call(['samtools', 'sort', '-n', '-o', bam_sort, bam], check=True)
    subprocess.run.assert_any_call(['bedtools', 'bamtobed', '-bedpe', '-mate1', '-i', bam_sort], stdout=ANY, check=True)
    assert not os.path.exists(bam_sort)
    assert os.path.exists(bedpe)


def test_bedpe2bed(testdir):
    bedpe = 'POLR2A.bedpe'
    sample_bedpe = Path(__file__).parent.joinpath('sample.bedpe')
    copyfile(sample_bedpe, bedpe)
    merge = bedpe + '-merge.bed'
    bed = 'POLR2A.bed'
    Bed.sort = MagicMock(side_effect=copyfile)
    bb.bedpe2bed(bedpe, bed)
    Bed.sort.assert_called_with(merge, bed)
    assert not os.path.exists(merge)
    assert os.path.exists(bed)
    with open(bed, 'r') as infile:
        assert infile.readline() == 'chr1\t100\t250\ttest1\t1\t+\t0\ttest0\n'
        assert infile.readline() == 'chr2\t300\t450\ttest2\t2\t+\t1\ttest1\n'
        assert infile.readline() == 'chr3\t500\t650\ttest3\t3\t+\t2\ttest2\n'
        assert infile.readline() == 'chr4\t700\t850\ttest4\t4\t+\t3\ttest3\n'
        assert infile.readline() == 'chr5\t100\t250\ttest5\t1\t-\t4\ttest4\n'
        assert infile.readline() == 'chr6\t300\t450\ttest6\t2\t-\t5\ttest5\n'
        assert infile.readline() == 'chr7\t500\t650\ttest7\t3\t-\t6\ttest6\n'
        assert infile.readline() == 'chr8\t700\t850\ttest8\t4\t-\t7\ttest7\n'
