# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Common constants and types used during ingestion/outgestion."""

from typing import Union

import h5py
import scipy.sparse as sp
from anndata._core.sparse_dataset import SparseDataset

from tiledbsoma._types import NPNDArray

SparseMatrix = Union[sp.csr_matrix, sp.csc_matrix, SparseDataset]
DenseMatrix = Union[NPNDArray, h5py.Dataset]
Matrix = Union[DenseMatrix, SparseMatrix]

# Arrays of strings from AnnData's uns are stored in SOMA as SOMADataFrame,
# since SOMA ND arrays are necessarily arrays *of numbers*. This is okay since
# the one and only job of SOMA uns is to faithfully ingest from AnnData and
# outgest back. These are parameters common to ingest and outgest of these.
_UNS_OUTGEST_HINT_KEY = "soma_uns_outgest_hint"
_UNS_OUTGEST_HINT_1D = "array_1d"
_UNS_OUTGEST_HINT_2D = "array_2d"
_UNS_OUTGEST_COLUMN_NAME_1D = "values"
_UNS_OUTGEST_COLUMN_PREFIX_2D = "values_"

_TILEDBSOMA_TYPE = "soma_tiledbsoma_type"
_DATAFRAME_ORIGINAL_INDEX_NAME_JSON = "soma_dataframe_original_index_name"
