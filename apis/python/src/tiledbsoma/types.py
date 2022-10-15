import pathlib
from typing import List, Literal, Sequence, Tuple, Union

import pandas as pd
import pyarrow as pa

Path = Union[str, pathlib.Path]

Ids = Union[List[str], List[bytes], List[int]]

Labels = Union[Sequence[str], pd.Index]

NTuple = Tuple[int, ...]

SOMABatchFormat = Literal["dense", "coo", "csr", "csc", "record-batch", "table"]
SOMAReadPartitions = Literal["IofN"]
SOMABatchSize = Literal["count", "size", "auto"]
SOMAResultOrder = Literal["row-major", "column-major", "unordered", "rowid-ordered"]

ArrowReadResult = Union[
    pa.Table,
    pa.RecordBatch,
    pa.Tensor,
    pa.SparseCOOTensor,
    pa.SparseCSRMatrix,
    pa.SparseCSCMatrix,
]

SOMADenseCoordinates = Union[int, slice]
SOMADenseNdCoordinates = Tuple[SOMADenseCoordinates, ...]

SOMASparseCoordinates = Union[int, slice, Tuple[int, ...], List[int], pa.IntegerArray]
SOMASparseNdCoordinates = Tuple[SOMADenseCoordinates, ...]
