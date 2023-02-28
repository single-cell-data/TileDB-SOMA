import pathlib
from typing import TYPE_CHECKING, Any, List, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from somacore import types
from typing_extensions import Literal

if TYPE_CHECKING:
    NPInteger = np.integer[npt.NBitBase]
    NPFloating = np.floating[npt.NBitBase]
    NPNDArray = npt.NDArray[np.number[npt.NBitBase]]
    PDSeries = pd.Series[Any]
else:
    NPInteger = np.integer
    NPFloating = np.floating
    NPNDArray = np.ndarray
    PDSeries = pd.Series


Path = Union[str, pathlib.Path]

Ids = Union[List[str], List[bytes], List[int]]

Labels = Union[Sequence[str], pd.Index]

NTuple = Tuple[int, ...]

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

# Re-exporting things from the somacore types namespace here.
Slice = types.Slice
is_nonstringy_sequence = types.is_nonstringy_sequence
is_slice_of = types.is_slice_of
