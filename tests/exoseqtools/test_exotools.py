import logging
from pathlib import Path
import pytest
from unittest.mock import MagicMock, ANY

import click
from click.testing import CliRunner
from exoseqtools import exotools, MoveAnnnotations, RemoveSecondMateBam


@pytest.fixture
def mock_testclass():
    moveannotations_samples = MoveAnnnotations.moveannotations_samples
    removesecondmate_samples = RemoveSecondMateBam.removesecondmate_samples
    yield
    MoveAnnnotations.moveannotations_samples = moveannotations_samples
    RemoveSecondMateBam.removesecondmate_samples = removesecondmate_samples


def test_moveannotations(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    distance = 32
    MoveAnnnotations.moveannotations_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(exotools.exotools, ['moveannotations', '--samples', samples, '--distance', distance])
    assert result.exit_code == 0
    MoveAnnnotations.moveannotations_samples.assert_called_once_with(samples, distance, True, True, None)


def test_removesecondmate(testdir, mock_testclass):
    samples = Path(__file__).parent.joinpath('samples.txt')
    RemoveSecondMateBam.removesecondmate_samples = MagicMock()
    runner = CliRunner()
    result = runner.invoke(exotools.exotools, ['removesecondmate', '--samples', samples])
    assert result.exit_code == 0
    RemoveSecondMateBam.removesecondmate_samples.assert_called_once_with(samples, 1, None)
