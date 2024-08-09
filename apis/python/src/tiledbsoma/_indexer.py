from __future__ import annotations

from typing import TYPE_CHECKING, Any, List, Optional, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from somacore.query.types import IndexLike

from tiledbsoma import pytiledbsoma as clib

from ._types import PDSeries

if TYPE_CHECKING:
    from .options import SOMATileDBContext

IndexerDataType = Union[
    npt.NDArray[np.int64],
    pa.Array,
    pa.IntegerArray,
    PDSeries,
    pd.arrays.IntegerArray,
    pa.ChunkedArray,
    List[int],
]


def tiledbsoma_build_index(
    data: IndexerDataType, *, context: Optional["SOMATileDBContext"] = None
) -> IndexLike:
    """Initialize re-indexer for provided indices.

    Deprecated. Provides the same functionality as the``IntIndexer`` class.

    Args:
       data:
           Integer keys used to build the index (hash) table.
       context:
           ``SOMATileDBContext`` object containing concurrecy level.

    Lifecycle:
        Deprecated.
    """

    return IntIndexer(data, context=context)


class IntIndexer:
    """A re-indexer for unique integer indices.

    Lifecycle:
        Maturing.
    """

    def __init__(
        self, data: IndexerDataType, *, context: Optional["SOMATileDBContext"] = None
    ):
        """Initialize re-indexer for provided indices.

        Args:
           data:
               Integer keys used to build the index (hash) table.
           context:
               ``SOMATileDBContext`` object containing concurrecy level.

        Lifecycle:
            Maturing.
        """
        self._context = context
        self._reindexer = clib.IntIndexer(
            None if self._context is None else self._context.native_context
        )
        self._reindexer.map_locations(data)

    def get_indexer(self, target: IndexerDataType) -> Any:
        """Compute underlying indices of index for target data.

        Compatible with Pandas' Index.get_indexer method.

        Args:
            target: Data to return re-index data for.
        """
        return (
            self._reindexer.get_indexer_pyarrow(target)
            if isinstance(target, (pa.Array, pa.ChunkedArray))
            else self._reindexer.get_indexer_general(target)
        )
