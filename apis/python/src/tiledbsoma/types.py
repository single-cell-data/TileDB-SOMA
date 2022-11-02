import pathlib
from typing import List, Literal, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import pyarrow as pa

Path = Union[str, pathlib.Path]

Ids = Union[List[str], List[bytes], List[int]]

Labels = Union[Sequence[str], pd.Index]

NTuple = Tuple[int, ...]

BatchFormat = Literal["dense", "coo", "csr", "csc", "record-batch", "table"]
ReadPartitions = Literal["IofN"]
BatchSize = Literal["count", "size", "auto"]
ResultOrder = Literal["row-major", "column-major", "unordered", "rowid-ordered"]

ArrowReadResult = Union[
    pa.Table,
    pa.RecordBatch,
    pa.Tensor,
    pa.SparseCOOTensor,
    pa.SparseCSRMatrix,
    pa.SparseCSCMatrix,
]

DenseCoordinates = Union[int, slice, pa.Array]
DenseNdCoordinates = Tuple[DenseCoordinates, ...]

SparseCoordinates = Union[int, slice, Tuple[int, ...], List[int], pa.IntegerArray]
SparseNdCoordinates = Tuple[DenseCoordinates, ...]

# TODO: add support for non-ints once the libtiledbsoma SOMAReader class has supports
# for non-ints. See also:
# https://github.com/single-cell-data/TileDB-SOMA/issues/418
# https://github.com/single-cell-data/TileDB-SOMA/issues/419
SparseIndexedDataFrameCoordinate = Optional[
    Union[int, slice, Sequence[int], pa.Array, pa.ChunkedArray, np.ndarray]
]
SparseIndexedDataFrameCoordinates = Sequence[SparseIndexedDataFrameCoordinate]
