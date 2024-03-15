"""
This file is separate from _util.py, due to a circular-import issue with
SOMATileDBContext which would otherwise ensue.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union

import numpy as np
import pyarrow as pa
from somacore.query.types import IndexLike

from tiledbsoma import pytiledbsoma as clib

if TYPE_CHECKING:
    from .options import SOMATileDBContext


def tiledbsoma_build_index(
    keys: Union[np.typing.NDArray[np.int64], pa.Array],
    *,
    context: Optional["SOMATileDBContext"] = None,
    thread_count: int = 4,
) -> IndexLike:
    """
    Returns an ``IndexLike`` re-indexer.
    The reindexer has an API similar to :meth:`pd.Index.get_indexer`

    Args:
        keys:
           Integer keys used to build the index (hash) table.
       context:
           ``SOMATileDBContext`` object containing concurrecy level (exclusive with thread_count).
       thread_count:
           Concurrency level when the user does not want to use ``context``.

    Lifecycle:
        Experimental.
    """
    native_context = None if context is None else context.native_context
    reindexer = clib.IntIndexer(native_context)
    reindexer.map_locations(keys)
    return reindexer  # type: ignore[no-any-return]
