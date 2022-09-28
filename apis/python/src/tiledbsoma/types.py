import pathlib
from typing import List, Literal, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp

Path = Union[str, pathlib.Path]

Ids = Union[List[str], List[bytes], List[int]]

Labels = Union[Sequence[str], pd.Index]

Matrix = Union[np.ndarray, sp.csr_matrix, sp.csc_matrix]

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
