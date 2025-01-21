# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

from typing import List, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from somacore.query.types import IndexLike

from tiledbsoma import pytiledbsoma as clib

from ._types import PDSeries
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
    data: IndexerDataType, *, context: SOMATileDBContext | None = None
) -> IndexLike:
    """Initialize re-indexer for provided indices (deprecated).

    Provides the same functionality as the``IntIndexer`` class.

    Args:
       data:
           Integer keys used to build the index (hash) table.
       context:
           ``SOMATileDBContext`` object containing concurrency level.

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
        self, data: IndexerDataType, *, context: SOMATileDBContext | None = None
    ):
        """Initialize re-indexer for provided indices.

        Args:
           data:
               Integer keys used to build the index (hash) table.
           context:
               ``SOMATileDBContext`` object containing concurrency level.

        Lifecycle:
            Maturing.
        """
        self._context = context
        self._reindexer = (
            clib.IntIndexer()
            if self._context is None
            else clib.IntIndexer(self._context.native_context)
        )

        # TODO: the map_locations interface does not accept chunked arrays. It would
        # save a copy (reduce memory usage) if they were natively supported.
        if isinstance(data, (pa.Array, pa.ChunkedArray)):
            data = data.to_numpy()
        elif isinstance(data, list):
            data = np.array(data, dtype=np.int64)
        elif isinstance(data, (pd.arrays.IntegerArray, pd.Series)):
            data = data.to_numpy(dtype=np.int64, copy=False)

        self._reindexer.map_locations(data)

    def get_indexer(self, target: IndexerDataType) -> npt.NDArray[np.intp]:
        """Compute underlying indices of index for target data.

        Compatible with Pandas' Index.get_indexer method.

        Args:
            target: Data to return re-index data for.
        """
        if isinstance(target, (pa.ChunkedArray, pa.Array)):
            return self._reindexer.get_indexer_pyarrow(target)

        if isinstance(target, (pd.arrays.IntegerArray, pd.Series)):
            target = target.to_numpy(dtype=np.int64, copy=False)
        elif isinstance(target, list):
            target = np.array(target, dtype=np.int64)

        return self._reindexer.get_indexer_general(target)
