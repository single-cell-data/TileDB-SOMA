# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

from __future__ import annotations

import datetime
import pathlib
from typing import TYPE_CHECKING, Any, List, Sequence, Tuple, Union, get_args

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from somacore import types
from typing_extensions import Literal

if TYPE_CHECKING:
    # `pd.{Series,Index}` require type parameters iff `pandas>=2`. Our pandas dependency (in `setup.py`) is unpinned,
    # which generally resolves to `pandas>=2`, but may be pandas<2 if something else in the user's environment requires
    # that. For type-checking purposes, `.pre-commit-config.yaml` specifies `pandas-stubs>=2`, and we type-check against
    # the `pandas>=2` types here.
    PDSeries = pd.Series[Any]
    PDIndex = pd.Index[Any]

    NPInteger = np.integer[npt.NBitBase]
    NPFloating = np.floating[npt.NBitBase]
    NPNDArray = npt.NDArray[np.number[npt.NBitBase]]
else:
    # When not-type-checking, but running with `pandas>=2`, the "missing" type-params don't affect anything.
    PDSeries = pd.Series
    PDIndex = pd.Index

    # Type subscription requires python >= 3.9, and we currently only type-check against 3.11.
    # TODO: remove these (and unify around subscripted types above) when we drop support for 3.8.
    NPInteger = np.integer
    NPFloating = np.floating
    # This alias likely needs to remain special-cased, even in Python ≥3.11, as tests pass the `Matrix` type alias
    # (which includes `NPNDArray` via `DenseMatrix`) to `isinstance`, causing error "argument 2 cannot be a
    # parameterized generic".
    NPNDArray = np.ndarray

Path = Union[str, pathlib.Path]

Ids = Union[List[str], List[bytes], List[int]]

Labels = Union[Sequence[str], PDIndex]

NTuple = Tuple[int, ...]

IngestMode = Literal["write", "schema_only", "resume"]  # for static-analysis checks
INGEST_MODES = get_args(IngestMode)  # for run-time checks

# Internal version of ``IngestMode`` that includes "update"; see ``IngestionParams``.
_IngestMode = Union[IngestMode, Literal["update"]]
_INGEST_MODES = INGEST_MODES + ("update",)


OpenTimestamp = Union[int, datetime.datetime]
"""Types that can be used as a timestamp to open a TileDB object.

Integers are treated as milliseconds since the Unix epoch.
"""

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

Metadatum = Union[bytes, float, int, str]
METADATA_TYPES = (bytes, float, int, str)
