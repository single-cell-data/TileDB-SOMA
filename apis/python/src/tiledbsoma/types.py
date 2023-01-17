import pathlib
from typing import TYPE_CHECKING, Any, List, Mapping, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from typing_extensions import Literal

if TYPE_CHECKING:
    NPInteger = np.integer[npt.NBitBase]
    NPFloating = np.floating[npt.NBitBase]
    NDArray = npt.NDArray[np.number[npt.NBitBase]]
else:
    NPInteger = np.integer
    NPFloating = np.floating
    NDArray = np.ndarray


Path = Union[str, pathlib.Path]

Ids = Union[List[str], List[bytes], List[int]]

Labels = Union[Sequence[str], pd.Index]

NTuple = Tuple[int, ...]

BatchFormat = Literal["dense", "coo", "csr", "csc", "record-batch", "table"]
ReadPartitions = Literal["IofN"]
BatchSize = Literal["count", "size", "auto"]
ResultOrder = Literal["row-major", "column-major", "auto"]
IngestMode = Literal["write", "schema_only", "resume"]  # for static-analysis checks
INGEST_MODES = ("write", "schema_only", "resume")  # for run-time checks

ArrowReadResult = Union[
    pa.Table,
    pa.RecordBatch,
    pa.Tensor,
    pa.SparseCOOTensor,
    pa.SparseCSRMatrix,
    pa.SparseCSCMatrix,
]

DenseCoordinates = Union[int, slice]
DenseNdCoordinates = Tuple[DenseCoordinates, ...]

# TODO: add support for non-ints once the libtiledbsoma SOMAReader class has supports
# for non-ints. See also:
# https://github.com/single-cell-data/TileDB-SOMA/issues/418
# https://github.com/single-cell-data/TileDB-SOMA/issues/419
#
# Note: we intentionally use `Union[None, ...]` in place of `Optional[...]` since
# we choose to emphasize that the argument-slots this is used in are not "optional"
# arguments -- they're required arguments, which can take the `None` value.
SparseDataFrameCoordinate = Union[
    None,
    int,
    slice,
    Sequence[int],
    pa.Array,
    pa.ChunkedArray,
    npt.NDArray[NPInteger],
]
SparseDataFrameCoordinates = Sequence[SparseDataFrameCoordinate]

# Note: we intentionally use `Union[None, ...]` in place of `Optional[...]` since
# we choose to emphasize that the argument-slots this is used in are not "optional"
# arguments -- they're required arguments, which can take the `None` value.

SparseNdCoordinates = Union[
    None,
    Sequence[
        Union[
            None,
            DenseCoordinates,
            Sequence[int],
            npt.NDArray[NPInteger],
            pa.IntegerArray,
        ]
    ],
]

PlatformConfig = Mapping[str, Any]
"""The platform-configuration dictionary. May contain a ``tiledb`` key."""
