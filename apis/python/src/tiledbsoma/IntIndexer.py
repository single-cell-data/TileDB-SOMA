from typing import Any

import numpy as np
import pandas as pd
import tiledb
import numpy.typing as npt

from . import pytiledbsoma as clib
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context

'''External API'''
class IntIndexer:

    @staticmethod
    def map_locations(keys: np.ndarray[np.int64, Any]) -> clib.IntIndexer:
        return indexer_map_locations(keys)

'''For internal use'''
def indexer_map_locations(keys: np.ndarray[np.int64, Any]) -> clib.IntIndexer:
    if len(np.unique(keys)) != len(keys):
        raise pd.errors.InvalidIndexError(
            "Reindexing only valid with uniquely valued Index objects"
        )
    context: SOMATileDBContext = _validate_soma_tiledb_context(
        SOMATileDBContext(tiledb.default_ctx())
    )
    compute_concurrency: int = 5
    if context._tiledb_ctx:
        compute_concurrency = int(
            int(context._tiledb_ctx.config()["sm.compute_concurrency_level"]) / 2
        )
    compute_concurrency = 8
    thread_counts = int(compute_concurrency)
    reindexer = clib.IntIndexer()
    reindexer.map_locations(keys, thread_counts)
    return reindexer
