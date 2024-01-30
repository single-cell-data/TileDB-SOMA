"""
This file is separate from _util.py, due to a circular-import issue with
SOMATileDBContext which would otherwise ensue.
"""

from typing import Optional

import numpy as np
import pandas as pd
from somacore.query.types import IndexLike

from tiledbsoma import pytiledbsoma as clib

from .options import SOMATileDBContext


def tiledbsoma_build_index(
    keys: np.typing.NDArray[np.int64],
    *,
    context: Optional[SOMATileDBContext] = None,
    thread_count: int = 4,
) -> IndexLike:
    """Builds an indexer object compatible with :meth:`pd.Index.get_indexer`."""
    if len(np.unique(keys)) != len(keys):
        raise pd.errors.InvalidIndexError(
            "Reindexing only valid with uniquely valued Index objects"
        )
    if context is not None:
        tdb_concurrency = int(
            context.tiledb_ctx.config().get("sm.compute_concurrency_level", 10)
        )
        thread_count = tdb_concurrency // 2

    reindexer = clib.IntIndexer()
    reindexer.map_locations(keys, thread_count)
    return reindexer  # type: ignore[no-any-return]
