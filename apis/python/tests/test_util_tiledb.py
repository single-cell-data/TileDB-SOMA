import pytest

import tiledbsoma


def test_stats(capfd: pytest.CaptureFixture[str]):
    """Make sure these exist, don't throw, and write correctly."""
    tiledbsoma.stats_enable()
    tiledbsoma.stats_dump()
    tiledbsoma.stats_reset()
    tiledbsoma.stats_disable()
    stdout, stderr = capfd.readouterr()
    assert stdout != ""
    assert stderr == ""
