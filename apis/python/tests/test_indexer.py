from __future__ import annotations

import threading

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

from tiledbsoma import SOMAError
from tiledbsoma._indexer import IntIndexer
from tiledbsoma.options import SOMATileDBContext
from tiledbsoma.options._soma_tiledb_context import _validate_soma_tiledb_context


@pytest.mark.parametrize(
    "keys, lookups",
    [
        ([-1, -1, -1, 0, 0, 0], np.tile(np.arange(1, 6), 4)),
        (np.tile(np.array([-1, 1, 2, 3, 4, 5]), 4), [-10000, 1, 2, 3, 5, 6]),
    ],
)
def test_duplicate_key_indexer_error(keys: np.array | list[int], lookups: np.array):
    context = _validate_soma_tiledb_context(SOMATileDBContext())
    with pytest.raises(SOMAError, match=r"There are duplicate keys."):
        IntIndexer(keys, context=context)

    pd_index = pd.Index(keys)
    with pytest.raises(pd.errors.InvalidIndexError):
        pd_index.get_indexer(lookups)


@pytest.mark.parametrize(
    "contextual",
    [True, False],
)
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
        (np.array(range(1, 10000)), np.array(range(1, 10000))),
        (pa.array(range(1, 10000)), pa.array(range(1, 10000))),
        (pd.array(range(1, 10000)), pd.array(range(1, 10000))),
        (
            pa.chunked_array(
                [
                    list(range(1, 10000)),
                    list(range(10000, 20000)),
                    list(range(30000, 40000)),
                ],
            ),
            pa.chunked_array(
                [
                    list(range(1, 10000)),
                    list(range(10000, 20000)),
                    list(range(30000, 40000)),
                ],
            ),
        ),
        (
            pd.Series(list(range(1, 10000)), copy=False),
            pd.Series(list(range(1, 10000)), copy=False),
        ),
    ],
)
def test_indexer(contextual: bool, keys: np.array, lookups: np.array):
    context = _validate_soma_tiledb_context(SOMATileDBContext()) if contextual else None
    all_results = []
    num_threads = 10

    def target():
        try:
            indexer = IntIndexer(keys, context=context)
            results = indexer.get_indexer(lookups)
            all_results.append(results)
        except Exception as e:
            all_results.append(e)

    for t in range(num_threads):
        thread = threading.Thread(target=target, args=())
        thread.start()
        thread.join()
    panda_indexer = pd.Index(keys)
    panda_results = panda_indexer.get_indexer(lookups)
    for i in range(num_threads):
        if isinstance(all_results[i], Exception):
            raise all_results[i]
        np.testing.assert_equal(all_results[i].all(), panda_results.all())


def test_expected_errors() -> None:
    context = SOMATileDBContext()

    # ndim != 1
    with pytest.raises(ValueError):
        IntIndexer(np.array([[0, 1], [2, 3]], dtype=np.int64), context=context)
    with pytest.raises(ValueError):
        IntIndexer(np.array([0, 1, 2, 3], dtype=np.int64), context=context).get_indexer(
            np.array([[0, 1]], dtype=np.int64),
        )

    # non-native byte order (one of the two must fail)
    with pytest.raises(TypeError):
        IntIndexer(np.array([0, 1, 2, 3], dtype=np.dtype("<i8")), context=context)
        IntIndexer(np.array([0, 1, 2, 3], dtype=np.dtype(">i8")), context=context)
    with pytest.raises(TypeError):
        IntIndexer(np.array([0, 1, 2, 3], dtype=np.int64), context=context).get_indexer(
            np.array([0], dtype=np.dtype("<i8")),
        )
        IntIndexer(np.array([0, 1, 2, 3], dtype=np.int64), context=context).get_indexer(
            np.array([0], dtype=np.dtype(">i8")),
        )

    # unsupported dtype
    with pytest.raises(TypeError):
        IntIndexer(np.array([0, 1, 2], dtype=np.uint64), context=context)
    with pytest.raises(TypeError):
        IntIndexer(np.array([0, 1, 2], dtype=np.int64), context=context).get_indexer(
            np.array([0, 1, 2], dtype=np.float64),
        )


def test_arrow_offset_handling() -> None:
    """Ensure accurate handling of slice offset/length"""
    assert np.array_equal(
        IntIndexer([0, 1]).get_indexer(pa.array([0, 1, 2]).slice(offset=1, length=1)),
        np.array([1]),
    )
    assert np.array_equal(
        IntIndexer([0, 1, 2, 3, 4, 5]).get_indexer(pa.chunked_array([[0, 1], [2, 3]]).slice(offset=1, length=2)),
        np.array([1, 2]),
    )
