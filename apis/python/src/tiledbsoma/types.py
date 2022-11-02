import pathlib
from typing import List, Literal, Sequence, Tuple, Union

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
SparseNdCoordinates = Sequence[DenseCoordinates]
