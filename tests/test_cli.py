from typer.testing import CliRunner

from masif_mimicry import __version__
from masif_mimicry.cli import app


def test_version_is_a_top_level_option():
    result = CliRunner().invoke(app, ["--version"])

    assert result.exit_code == 0
    assert result.stdout.strip() == __version__
