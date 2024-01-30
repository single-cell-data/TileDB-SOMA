import numpy as np

from tiledbsoma import SOMATileDBContext, build_index


def test_reindexer_api_thread_count():
    keys = np.arange(3, 10, 2)
    ids = np.arange(3, 10, 2)
    expected = np.array([0, 1, 2, 3])
    indexer = build_index(keys)
    result = indexer.get_indexer(ids)
    assert np.equal(result.all(), expected.all())


def test_reindexer_api_context():
    context = SOMATileDBContext()
    keys = np.arange(3, 10, 2)
    ids = np.arange(3, 10, 2)
    expected = np.array([0, 1, 2, 3])
    indexer = build_index(keys, context)
    result = indexer.get_indexer(ids)
    assert np.equal(result.all(), expected.all())
