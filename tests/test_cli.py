"""click_mergevcfs cli tests."""

from click.testing import CliRunner
from click_mergevcfs import cli

def test_main():
    """Test for the main command."""
    runner = CliRunner()
    result = runner.invoke(cli.main, ["--help"])
    assert result.exit_code == 0
