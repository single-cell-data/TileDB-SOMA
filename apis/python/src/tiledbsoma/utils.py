import numpy as np
import pandas as pd

from tiledbsoma import pytiledbsoma as clib

from .options import SOMATileDBContext

"""Builds an indexer object compatible with :meth:`pd.Index.get_indexer`."""


def build_index(
    keys: np.typing.NDArray[np.int64], context: SOMATileDBContext
) -> clib.IntIndexer:
    if len(np.unique(keys)) != len(keys):
        raise pd.errors.InvalidIndexError(
            "Reindexing only valid with uniquely valued Index objects"
        )
    tdb_concurrency = int(
        context.tiledb_ctx.config().get("sm.compute_concurrency_level", 10)
    )
    thread_count = tdb_concurrency // 2
    reindexer = clib.IntIndexer()
    reindexer.map_locations(keys, thread_count)
    return reindexer
