from typing import List, Union

import numpy as np
import pandas as pd
import pytest

from tiledbsoma._index_util import tiledbsoma_build_index
from tiledbsoma.options import SOMATileDBContext
from tiledbsoma.options._soma_tiledb_context import _validate_soma_tiledb_context


@pytest.mark.parametrize(
    "keys, lookups",
    [
        ([-1, -1, -1, 0, 0, 0], np.tile(np.arange(1, 6), 4)),
        (np.tile(np.array([-1, 1, 2, 3, 4, 5]), 4), [-10000, 1, 2, 3, 5, 6]),
    ],
)
def test_duplicate_key_indexer_error(
    keys: Union[np.array, List[int]], lookups: np.array
):
    context = _validate_soma_tiledb_context(SOMATileDBContext())
    with pytest.raises(RuntimeError, match="There are duplicate keys."):
        tiledbsoma_build_index(keys, context=context)

    pd_index = pd.Index(keys)
    with pytest.raises(pd.errors.InvalidIndexError):
        pd_index.get_indexer(lookups)


@pytest.mark.parametrize(
    "keys, lookups",
    [
        ([1], [1, 1, 1, 1]),
        ([-1, 1, 2, 3, 4, 5], np.tile(np.array([-1, 1, 2, 3, 4, 5]), 4)),
        ([-10000, -100000, 200000, 5, 1, 7], np.tile(np.array([-1, 1, 2, 3, 4, 5]), 4)),
        (
            [-10000, -200000, 1000, 3000, 1, 2],
            np.tile(np.array([-1, 1, 2, 3, 4, 5]), 4),
        ),
        (list(range(1, 10000)), list(range(1, 10))),
        (
            list(range(1, 10000)),
            [
                525,
                1293,
                1805,
                5802,
                7636,
                7754,
                7791,
                7957,
                7959,
                8067,
                8340,
                8736,
                8806,
                9329,
                9377,
                9653,
            ],
        ),
        (list(range(1, 10000)), list(range(1, 10000))),
    ],
)
def test_indexer(keys: np.array, lookups: np.array):
    context = _validate_soma_tiledb_context(SOMATileDBContext())
    indexer = tiledbsoma_build_index(keys, context=context)
    results = indexer.get_indexer(lookups)
    panda_indexer = pd.Index(keys)
    panda_results = panda_indexer.get_indexer(lookups)
    np.testing.assert_equal(results.all(), panda_results.all())
