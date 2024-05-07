from time import perf_counter

import numpy as np
import pandas as pd
import pyarrow as pa

from tiledbsoma._indexer import tiledbsoma_build_index
from tiledbsoma.options import SOMATileDBContext
from tiledbsoma.options._soma_tiledb_context import _validate_soma_tiledb_context

"""
Performance test evaluating the reindexer performance compared to pandas.Index for different data types
"""


def build(keys, pandas):
    context = _validate_soma_tiledb_context(SOMATileDBContext())
    if pandas:
        indexer = pd.Index(keys)
        indexer.get_indexer([1])
    else:
        perf_counter()
        indexer = tiledbsoma_build_index(keys, context=context.native_context)
    return indexer


def lookup(indexer, lookups):
    results = indexer.get_indexer(lookups)
    return results


def run(data_type, keys, lookups, pandas):
    start_time = perf_counter()
    indexer = build(keys, pandas)
    build_time = perf_counter()
    lookup(indexer, lookups)
    lookup_time = perf_counter()
    if pandas:
        name = "pandas"
    else:
        name = "reindexer"
    print(
        f"{data_type}: Setup time: {name}: {build_time - start_time}: Lookup time: {lookup_time - build_time}"
    )


def indexer_test_build(data_type, keys, lookups):
    run(data_type, keys, lookups, False)
    run(data_type, keys, lookups, True)


def main():
    # Creating keys and lookup values of different types (np.array, pa.array, pa.chunked_array, pa.array, pd.Series,
    # pd.array and python list)  and try the re-indexer and pandas.Index
    keys = np.array(list(range(1, 100000000)), dtype=np.int64)
    lookups = np.random.randint(0, 1000, 1000000000)
    indexer_test_build("np.array", keys, lookups)

    keys = pa.array(list(range(1, 100000000)))
    lookups = pa.array(np.random.randint(0, 1000, 1000000000).astype(np.int64))
    indexer_test_build("pa.array", keys, lookups)

    keys = pa.chunked_array(
        [
            list(range(1, 10000000)),
            list(range(10000000, 20000000)),
            list(range(20000000, 30000000)),
            list(range(30000000, 40000000)),
            list(range(40000000, 50000000)),
            list(range(50000000, 6000000)),
            list(range(60000000, 70000000)),
            list(range(70000000, 80000000)),
            list(range(80000000, 90000000)),
            list(range(90000000, 100000000)),
        ]
    )
    lookups = []
    for x in range(10):
        lookups.append(list(np.random.randint(0, 1000, 100000000, dtype=np.int64)))
    lookups = pa.chunked_array(lookups)
    indexer_test_build("pa.chunked_array", keys, lookups)

    keys = pd.Series(list(range(1, 100000000)))
    lookups = pd.Series(np.random.randint(0, 1000, 1000000000))
    indexer_test_build("pd.Series", keys, lookups)

    keys = pd.array(list(range(1, 100000000)))
    lookups = pd.array(np.random.randint(0, 1000, 1000000000, dtype=np.int64))
    indexer_test_build("pd.array", keys, lookups)

    keys = list(range(1, 100000000))
    lookups = list(np.random.randint(0, 1000, 1000000000))


#   Commented out as it takes a very long time for pandas
#   indexer_test_build("python list", keys, lookups)


main()
