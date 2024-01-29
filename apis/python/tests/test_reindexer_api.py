import numpy as np

from tiledbsoma import build_index


def main():
    keys = np.arange(3, 10, 2)
    ids = np.arange(3, 10, 2)
    expected = np.array([0, 1, 2, 3])
    indexer = build_index(keys)
    result = indexer.get_indexer(ids)
    assert np.equal(result.all(), expected.all())


main()
