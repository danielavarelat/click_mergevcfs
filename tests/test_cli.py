"""=click_mergevcfs cli tests."""

from click.testing import CliRunner
import pytest

from click_mergevcfs import cli

def test_main():
    """Sample test for main command."""
    with pytest.raises(SystemExit) as _:
        cli.main()

