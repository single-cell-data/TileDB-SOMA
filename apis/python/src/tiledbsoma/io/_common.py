# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Common constants and types used during ingestion/outgestion."""

from __future__ import annotations

from typing import Dict, List, Mapping, Union

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

# As of anndata 0.11 we get a warning importing anndata.experimental.
# But anndata.abc doesn't exist in anndata 0.10. And anndata 0.11 doesn't
# exist for Python 3.9. And we have not yet dropped support for Python 3.9.
try:
    from anndata.abc import CSCDataset, CSRDataset
except (AttributeError, ModuleNotFoundError):
    from anndata.experimental import CSCDataset, CSRDataset

from tiledbsoma._types import Metadatum, NPNDArray

SparseMatrix = Union[sp.csr_matrix, sp.csc_matrix, CSCDataset, CSRDataset]
DenseMatrix = Union[NPNDArray, h5py.Dataset]
Matrix = Union[DenseMatrix, SparseMatrix]

UnsScalar = Union[str, int, float, np.generic]
# TODO: support sparse matrices in `uns`
UnsLeaf = Union[UnsScalar, List[UnsScalar], pd.DataFrame, NPNDArray]
UnsNode = Union[UnsLeaf, Mapping[str, "UnsNode"]]
UnsMapping = Mapping[str, UnsNode]
# Specialize `UnsNode` to `Dict` instead of `Mapping`
# `Mapping` doesn't expose `__set__`, so this is useful for building `uns` dictionaries.
UnsDictNode = Union[UnsLeaf, Dict[str, "UnsDictNode"]]
UnsDict = Dict[str, UnsDictNode]

AdditionalMetadata = Union[Dict[str, Metadatum], None]

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
