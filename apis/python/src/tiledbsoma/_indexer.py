"""
This file is separate from _util.py, due to a circular-import issue with
SOMATileDBContext which would otherwise ensue.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from somacore.query.types import IndexLike

from tiledbsoma import pytiledbsoma as clib

if TYPE_CHECKING:
    from .options import SOMATileDBContext


def tiledbsoma_build_index(
    keys: Union[  # type: ignore[type-arg]
        npt.NDArray[np.int64],
        pa.Array,
        pa.IntegerArray,
        pd.Series,
        pd.arrays.IntegerArray,
        pa.ChunkedArray,
        list[int],
    ],
    *,
    context: Optional["SOMATileDBContext"] = None,
) -> IndexLike:
    """
    Returns an ``IndexLike`` re-indexer.
    The reindexer has an API similar to :meth:`pd.Index.get_indexer`

    Args:
       keys:
           Integer keys used to build the index (hash) table.
       context:
           ``SOMATileDBContext`` object containing concurrecy level (exclusive with thread_count).

    Lifecycle:
        Experimental.
    """
    return IntIndex(keys, context=context)


class IntIndex:
    """TODO: Add IntIndex class docstring before merging."""

    def __init__(
        self,
        data: Union[  # type: ignore[type-arg]
            npt.NDArray[np.int64],
            pa.Array,
            pa.IntegerArray,
            pd.Series,
            pd.arrays.IntegerArray,
            pa.ChunkedArray,
            list[int],
        ],
        *,
        context: Optional["SOMATileDBContext"] = None,
    ):
        """TODO: Add IntIndex.__init__ docstring before merging.
            Returns an ``IndexLike`` re-indexer.
            The reindexer has an API similar to :meth:`pd.Index.get_indexer`

        Args:
           keys:
               Integer keys used to build the index (hash) table.
           context:
               ``SOMATileDBContext`` object containing concurrecy level (exclusive with
               thread_count).

        Lifecycle:
            Experimental.
        """
        self._context = context
        self._reindexer = clib.IntIndexer(
            None if self._context is None else self._context.native_context
        )
        self._reindexer.map_locations(data)

    def get_indexer(
        self,
        target: Union[  # type: ignore[type-arg]
            npt.NDArray[np.int64],
            pa.Array,
            pa.IntegerArray,
            pd.Series,
            pd.arrays.IntegerArray,
            pa.ChunkedArray,
            list[int],
        ],
    ) -> Any:
        """Something compatible with Pandas' Index.get_indexer method."""
        return (
            self._reindexer.get_indexer_pyarrow(target)
            if isinstance(target, (pa.Array, pa.ChunkedArray))
            else self._reindexer.get_indexer_general(target)
        )
