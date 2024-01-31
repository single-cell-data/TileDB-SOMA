"""
This file is separate from _util.py, due to a circular-import issue with
SOMATileDBContext which would otherwise ensue.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np
from somacore.query.types import IndexLike

from tiledbsoma import pytiledbsoma as clib

if TYPE_CHECKING:
    from .options import SOMATileDBContext


def tiledbsoma_build_index(
    keys: np.typing.NDArray[np.int64],
    *,
    context: Optional["SOMATileDBContext"] = None,
    thread_count: int = 4,
) -> IndexLike:
    """
    Returns a IndexLike re-indexer
    The reindxer has an API similar to :meth:`pd.Index.get_indexer`

    Args:
        keys:
           integer keys used to build the index (hash) table.
       context:
           tiledbcontext object containing concurrecy level (exclusive with thread_count).
       thread_count:
           concurrecy level when the user does not want to use `context`.

    Lifecycle:
        Experimental.
    """
    if context is not None:
        tdb_concurrency = int(
            context.tiledb_ctx.config().get("sm.compute_concurrency_level", 10)
        )
    tdb_concurrency = 10
    if context is not None:
        tdb_concurrency = int(
            context.tiledb_ctx.config().get("sm.compute_concurrency_level", 10)
        )
    thread_count = tdb_concurrency // 2
    reindexer = clib.IntIndexer()
    reindexer.map_locations(keys, thread_count)
    return reindexer  # type: ignore[no-any-return]
