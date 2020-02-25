import os
import shutil
import subprocess
import tempfile
import unittest

import click
from click.testing import CliRunner
from pathlib import Path
from seqtools import DownloadSample as ds
from unittest.mock import MagicMock


class TestDownloadSample(unittest.TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()
        os.chdir(self.test_dir)

    def tearDown(self):
        # Remove the directory after the test
        os.chdir(Path(__file__).parent)
        shutil.rmtree(self.test_dir)

    def create_file(*args):
        for name in args[1:]:
            with open(name, 'w') as outfile:
                outfile.write('test')

    def test_main(self):
        samples = Path(__file__).parent.joinpath('samples.txt')
        mem = '100MB'
        threads = 6
        subprocess.run = MagicMock(side_effect=[self.create_file('SRR8518913_1.fastq'), self.create_file('SRX5322424_1.fastq', 'SRX5322424_2.fastq'), self.create_file('SRR8518915_1.fastq', 'SRR8518915_2.fastq')])
        runner = CliRunner()
        result = runner.invoke(ds.main, ['-s', samples ])
        self.assertEqual(0, result.exit_code)
        subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRR8518913'], check=True)
        subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRX5322424'], check=True)
        subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRR8518915'], check=True)
        self.assertTrue(os.path.exists('POLR2A_1.fastq'))
        self.assertFalse(os.path.exists('POLR2A_2.fastq'))
        self.assertTrue(os.path.exists('ASDURF_1.fastq'))
        self.assertTrue(os.path.exists('ASDURF_2.fastq'))
        self.assertTrue(os.path.exists('POLR1C_1.fastq'))
        self.assertTrue(os.path.exists('POLR1C_2.fastq'))

    def test_main_allfast(self):
        samples = Path(__file__).parent.joinpath('samples.txt')
        mem = '200MB'
        threads = 2
        subprocess.run = MagicMock(side_effect=[self.create_file('SRR8518913_1.fastq'), self.create_file('SRX5322424_1.fastq', 'SRX5322424_2.fastq'), self.create_file('SRR8518915_1.fastq', 'SRR8518915_2.fastq')])
        runner = CliRunner()
        result = runner.invoke(ds.main, ['-s', samples, '--fast', '-e', str(threads), '-m', mem])
        self.assertEqual(0, result.exit_code)
        subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRR8518913'], check=True)
        subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRX5322424'], check=True)
        subprocess.run.assert_any_call(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRR8518915'], check=True)
        self.assertTrue(os.path.exists('POLR2A_1.fastq'))
        self.assertFalse(os.path.exists('POLR2A_2.fastq'))
        self.assertTrue(os.path.exists('ASDURF_1.fastq'))
        self.assertTrue(os.path.exists('ASDURF_2.fastq'))
        self.assertTrue(os.path.exists('POLR1C_1.fastq'))
        self.assertTrue(os.path.exists('POLR1C_2.fastq'))

    def test_main_allslow(self):
        samples = Path(__file__).parent.joinpath('samples.txt')
        mem = '200MB'
        threads = 2
        subprocess.run = MagicMock(side_effect=[self.create_file('SRR8518913_1.fastq'), self.create_file('SRX5322424_1.fastq', 'SRX5322424_2.fastq'), self.create_file('SRR8518915_1.fastq', 'SRR8518915_2.fastq')])
        runner = CliRunner()
        result = runner.invoke(ds.main, ['-s', samples, '--slow', '-e', str(threads), '-m', mem])
        self.assertEqual(0, result.exit_code)
        subprocess.run.assert_any_call(['fastq-dump', '--split-files', 'SRR8518913'], check=True)
        subprocess.run.assert_any_call(['fastq-dump', '--split-files', 'SRX5322424'], check=True)
        subprocess.run.assert_any_call(['fastq-dump', '--split-files', 'SRR8518915'], check=True)
        self.assertTrue(os.path.exists('POLR2A_1.fastq'))
        self.assertFalse(os.path.exists('POLR2A_2.fastq'))
        self.assertTrue(os.path.exists('ASDURF_1.fastq'))
        self.assertTrue(os.path.exists('ASDURF_2.fastq'))
        self.assertTrue(os.path.exists('POLR1C_1.fastq'))
        self.assertTrue(os.path.exists('POLR1C_2.fastq'))

    def test_main_second(self):
        samples = Path(__file__).parent.joinpath('samples.txt')
        mem = '100MB'
        threads = 6
        subprocess.run = MagicMock(side_effect=self.create_file('SRX5322424_1.fastq', 'SRX5322424_2.fastq'))
        runner = CliRunner()
        result = runner.invoke(ds.main, ['-s', samples, '-i', 1])
        self.assertEqual(0, result.exit_code)
        subprocess.run.assert_called_once_with(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, 'SRX5322424'], check=True)
        self.assertTrue(os.path.exists('ASDURF_1.fastq'))
        self.assertTrue(os.path.exists('ASDURF_2.fastq'))

    def test_download_singleend(self):
        srr = 'SRR8518913'
        mem = '100MB'
        threads = 1
        subprocess.run = MagicMock(side_effect=self.create_file(srr + '_1.fastq'))
        ds.download('POLR2A', srr, True, threads, mem)
        subprocess.run.assert_called_with(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, srr], check=True)
        self.assertTrue(os.path.exists('POLR2A_1.fastq'))

    def test_download_pairedend(self):
        srr = 'SRR8518913'
        mem = '100MB'
        threads = 1
        subprocess.run = MagicMock(side_effect=self.create_file(srr + '_1.fastq', srr + '_2.fastq'))
        ds.download('POLR2A', srr, True, threads, mem)
        subprocess.run.assert_called_with(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, srr], check=True)
        self.assertTrue(os.path.exists('POLR2A_1.fastq'))
        self.assertTrue(os.path.exists('POLR2A_2.fastq'))

    def test_download_slow_singleend(self):
        srr = 'SRR8518913'
        mem = '100MB'
        threads = 1
        subprocess.run = MagicMock(side_effect=self.create_file(srr + '_1.fastq'))
        ds.download('POLR2A', srr, False, threads, mem)
        subprocess.run.assert_called_with(['fastq-dump', '--split-files', srr], check=True)
        self.assertTrue(os.path.exists('POLR2A_1.fastq'))

    def test_download_slow_pairedend(self):
        srr = 'SRR8518913'
        mem = '100MB'
        threads = 1
        subprocess.run = MagicMock(side_effect=self.create_file(srr + '_1.fastq', srr + '_2.fastq'))
        ds.download('POLR2A', srr, False, threads, mem)
        subprocess.run.assert_called_with(['fastq-dump', '--split-files', srr], check=True)
        self.assertTrue(os.path.exists('POLR2A_1.fastq'))
        self.assertTrue(os.path.exists('POLR2A_2.fastq'))

    def test_download_failed(self):
        srr = 'SRR8518913'
        mem = '100MB'
        threads = 1
        subprocess.run = MagicMock(side_effect=subprocess.CalledProcessError('Could not download file', ['test']))
        self.assertRaises(subprocess.CalledProcessError, lambda:ds.download('POLR2A', srr, True, threads, mem))
        subprocess.run.assert_called_with(['fasterq-dump', '--split-files', '--threads', str(threads), '--mem', mem, srr], check=True)
        self.assertFalse(os.path.exists('POLR2A_1.fastq'))
        self.assertFalse(os.path.exists('POLR2A_2.fastq'))


if __name__ == '__main__':
    unittest.main()
