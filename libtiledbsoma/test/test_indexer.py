import numpy as np
import pandas as pd

from tiledbsoma.IntIndexer import IntIndexer

# tiledbsoma.pytiledbsoma.config_logging("debug")


# 1d array to list


def indexer_test(keys: np.array, lookups: np.array, fail: bool):
    if fail:
        indexer_test_fail(keys, lookups)
    else:
        indexer_test_pass(keys, lookups)


def indexer_test_fail(keys: np.array, lookups: np.array):
    reindexing_exception: bool = False
    try:
        indexer = IntIndexer.map_locations(keys)
        indexer.get_indexer(lookups)
    except pd.errors.InvalidIndexError:
        reindexing_exception = True

    try:
        panda_indexer = pd.Index(keys)
        panda_indexer.get_indexer(lookups)
    except pd.errors.InvalidIndexError as exception:
        if reindexing_exception:
            return
        raise exception


def indexer_test_pass(keys: np.array, lookups: np.array):
    indexer = IntIndexer.map_locations(keys)
    results = indexer.get_indexer(lookups)
    panda_indexer = pd.Index(keys)
    panda_results = panda_indexer.get_indexer(lookups)
    assert np.equal(results.all(), panda_results.all())


test_data = [
    {
        "keys": [-1, -1, -1, 0, 0, 0],
        "lookups": [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5],
        "pass": False,
    },
    {
        "keys": [
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
        ],
        "lookups": [-10000, 1, 2, 3, 5, 6],
        "pass": False,
    },
    {
        "keys": [-1, 1, 2, 3, 4, 5],
        "lookups": [
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
        ],
        "pass": True,
    },
    {
        "keys": [-10000, -100000, 200000, 5, 1, 7],
        "lookups": [
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
        ],
        "pass": True,
    },
    {
        "keys": [-10000, -200000, 1000, 3000, 1, 2],
        "lookups": [
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
            -1,
            1,
            2,
            3,
            4,
            5,
        ],
        "pass": True,
    },
    {
        "keys": [i for i in range(1, 10000)],
        "lookups": [i for i in range(1, 10)],
        "pass": True,
    },
    {
        "keys": [1, 2, 3, 4, 5, 2**63 - 1],
        "lookups": [i for i in range(2**63 - 1000, 2**63 - 1)],
        "pass": True,
    },
]


def main():
    for data in test_data:
        indexer_test(data["keys"], data["lookups"], not data["pass"])


if __name__ == "__main__":
    main()
