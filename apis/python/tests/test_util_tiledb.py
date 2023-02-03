import tiledbsoma


def test_stats():
    """Make sure these exist and don't throw"""
    tiledbsoma.stats_enable()
    tiledbsoma.stats_dump()
    tiledbsoma.stats_reset()
    tiledbsoma.stats_disable()
