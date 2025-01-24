import pytest
import click
from click.testing import CliRunner

from bin import cli
from bin.utils import split_shell_command


@pytest.mark.skip(reason="The CLI function uses sys.argv to check inputes,"
                  " which doesn't work well with invoking directly with click runner.invoke")
def test_cli():
    # Test the cli function
    
    args = "call annotate distill format-kegg"
    runner = CliRunner()

    result = runner.invoke(cli.cli, args)
    click.echo(result.output)
    assert result.exception
    assert "format-kegg must be run alone." in result.output