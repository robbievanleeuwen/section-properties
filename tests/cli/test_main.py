"""Test cases for the __main__ module."""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from sectionproperties import __main__


@pytest.fixture
def runner() -> CliRunner:
    """Fixture for invoking command-line interfaces.

    Returns:
        CliRunner object
    """
    return CliRunner()


def test_main_succeeds(runner: CliRunner) -> None:
    """It exits with a status code of zero.

    Args:
        runner: CliRunner object
    """
    result = runner.invoke(__main__.main)
    assert result.exit_code == 0
