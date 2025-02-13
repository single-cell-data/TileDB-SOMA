# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Conversion utility methods.
"""

from __future__ import annotations

from typing import TypeVar, cast

import numpy as np
import pandas as pd
import pandas._typing as pdt
import pandas.api.types
import pyarrow as pa
import scipy.sparse as sp

from .._fastercsx import CompressedMatrix
from .._funcs import typeguard_ignore
from .._types import NPNDArray, PDSeries
from ..options._soma_tiledb_context import SOMATileDBContext

_DT = TypeVar("_DT", bound=pdt.Dtype)
_MT = TypeVar("_MT", NPNDArray, sp.spmatrix, PDSeries)
_str_to_type = {"boolean": bool, "string": str, "bytes": bytes}

COLUMN_DECAT_THRESHOLD = 32767
"""
For enum-of-string columns with a cardinality higher than this, we convert from
enum-of-string in the AnnData ``obs``/``var``, to plain string in TileDB-SOMA
``obs``/``var``. However, if we're appending to existing storage, we follow the
schema there. Reasoning behind this choice: accommodate signed 16-bit index type.
See also https://github.com/single-cell-data/TileDB-SOMA/pull/3415.
"""


def obs_or_var_to_tiledb_supported_array_type(obs_or_var: pd.DataFrame) -> pd.DataFrame:
    """
    Performs a typecast into types that TileDB can persist.  This includes, as a
    performance improvement, converting high-cardinality categorical-of-string
    columns (cardinality > 4096) to plain string.
    """
    if len(obs_or_var.columns) == 0:
        return obs_or_var.copy()

    return pd.DataFrame.from_dict(
        {
            str(k): to_tiledb_supported_array_type(str(k), v)
            for k, v in obs_or_var.items()
        },
    )


@typeguard_ignore
def _to_tiledb_supported_dtype(dtype: _DT) -> _DT:
    """A handful of types are cast into the TileDB type system."""
    # TileDB has no float16 -- cast up to float32
    return cast(_DT, np.dtype("float32")) if dtype == np.dtype("float16") else dtype


def to_tiledb_supported_array_type(name: str, x: _MT) -> _MT:
    """Converts datatypes unrepresentable by TileDB into datatypes it can represent,
    e.g., float16 -> float32.
    """
    if isinstance(x, (np.ndarray, sp.spmatrix)) or not isinstance(
        x.dtype, pd.CategoricalDtype
    ):
        target_dtype = _to_tiledb_supported_dtype(x.dtype)
        return x if target_dtype == x.dtype else x.astype(target_dtype)

    # If the column is categorical-of-string of high cardinality, we declare
    # this is likely a mistake, and it will definitely lead to performance
    # issues in subsequent processing.
    if isinstance(x, pd.Series) and isinstance(x.dtype, pd.CategoricalDtype):
        # Heuristic number
        if len(x.cat.categories) > COLUMN_DECAT_THRESHOLD:
            return x.astype(x.cat.categories.dtype)

    return x


def csr_from_coo_table(
    tbl: pa.Table, num_rows: int, num_cols: int, context: SOMATileDBContext
) -> sp.csr_matrix:
    """Given an Arrow Table containing COO data, return a ``scipy.sparse.csr_matrix``."""
    s = CompressedMatrix.from_soma(
        tbl, (num_rows, num_cols), "csr", True, context
    ).to_scipy()
    return s


def df_to_arrow(df: pd.DataFrame) -> pa.Table:
    """
    Handle special cases where pa.Table.from_pandas is not sufficient.
    """
    nullable_fields = set()
    # Not for name, col in df.items() since we need df[k] on the left-hand sides
    for key in df:
        # Make attributes nullable. Context:
        # * df_to_arrow is _solely_ for use of tiledbsoma.io
        #   o Anyone calling the SOMA API directly has user-provided Arrow
        #     schema which must be respected
        #   o Anyone calling tiledbsoma.io -- including from_h5ad/from_anndata,
        #     and update_obs/update_var -- does not provide an Arrow schema
        #     explicitly.  We compute an Arrow schema for them here.
        # * Even when the _initial_ data is all non-null down a particular
        #   string column, there are two ways a _subsequent_ write can provide
        #   nulls: append-mode ingest, or, update_obs/update_var wherein the new
        #   data has nulls even when the data used at schema-create time was
        #   non-null.
        # * We have no way of knowing at initial ingest time whether or not
        #   users will later be appending, or updating, with null data.
        # * Note that Arrow has a per-field nullable flag in its schema metadata
        #   -- and so do TileDB array schemas.
        #
        # Note in particular this is for the use of tiledbsoma.io:
        #
        # * In the tiledbsoma API (e.g. DataFrame.create) the user passes an
        #   Arrow schema and we respect it as-is. They specify nullability, or
        #   not, as they wish.
        # * In tiledbsoma.io, the user-provided inputs are AnnData objects.
        #   We compute the Arrow schema _for_ them. And we must accommodate
        #   reasonable/predictable needs.

        nullable_fields.add(key)

        # Handle special cases for all null columns where the dtype is "object"
        # or "category" and must be explicitly casted to the correct pandas
        # extension dtype.
        #
        # Note: with
        #   anndata.obs['new_col'] = pd.Series(data=np.nan, dtype=np.dtype(str))
        # the dtype comes in to us via `tiledbsoma.io.from_anndata` not
        # as `pd.StringDtype()` but rather as `object`.
        #
        # Note: we're working around a subtle issue involving Pandas, and Arrow's
        # from_pandas, and categoricals.
        #
        # * If you do this:
        #     pd.Series(["a", "b", "c", "d"], dtype=pd.CategoricalDtype())
        #   then you get Pandas categorical of string with no nulls -- as desired.
        # * If you do this:
        #     pd.Series(["a", "b", None, "d"], dtype=pd.CategoricalDtype())
        #   or
        #     pd.Series(["a", "b", np.nan, "d"], dtype=pd.CategoricalDtype())
        #   then you get Pandas categorical of string, with some nulls -- as desired
        # * If you do this:
        #     pd.Series([None] * 4, dtype=pd.CategoricalDtype())
        #   or
        #     pd.Series([np.nan] * 4, dtype=pd.CategoricalDtype())
        #   then you get Pandas categorical of double -- with NaN values -- not as desired.
        if df[key].isnull().all():
            if df[key].dtype.name == "object":
                df[key] = pd.Series([None] * df.shape[0], dtype=pd.StringDtype())
            elif df[key].dtype.name == "category":
                # This is a trick to avoid getting float64 value type in the Pandas categorical.
                # That's the good news. The bad news is that pa.Table.from_pandas() of this
                # will result in Arrow value-type of pa.null().  Part two, to deal with
                # this, is below.
                df[key] = pd.Series(
                    ["X"] * df.shape[0], dtype=pd.CategoricalDtype()
                ).cat.remove_categories(["X"])

    # For categoricals, it's possible to get
    #   TypeError: Object of type bool_ is not JSON serializable
    # deep within library functions. Debugging reveals that this happens when
    # the df[key].values.ordered is of type np.bool_ rather than Python bool.
    # So, we cast and reconstruct.
    for key in df:
        column = df[key]
        if isinstance(column.dtype, pd.CategoricalDtype):
            if hasattr(column.values, "categories"):
                categories = column.values.categories

            if hasattr(column.values, "ordered"):
                ordered = bool(column.values.ordered)

            df[key] = pd.Categorical(
                values=column, categories=categories, ordered=ordered
            )

    arrow_table = pa.Table.from_pandas(df)

    md = arrow_table.schema.metadata
    md.update(dict.fromkeys(nullable_fields, "nullable"))
    arrow_table = arrow_table.replace_schema_metadata(md)

    # For tiledbsoma.io (for which this method exists) _any_ dataset can be appended to
    # later on. This means that on fresh ingest we must use a larger bit-width than
    # the bare minimum necessary.
    new_map = {}
    for field in arrow_table.schema:
        if field.name == "__index_level_0__":
            continue
        elif pa.types.is_dictionary(field.type):
            old_index_type = field.type.index_type
            new_index_type = (
                pa.int32()
                if old_index_type in [pa.int8(), pa.int16()]
                else old_index_type
            )
            # This is part two of what we need to do to get null-filled Pandas
            # categorical-of-string conveyed to Arrow. An entirely null-filled
            # Pandas categorical-of-string series, after py.Table.from_pandas(),
            # will have type pa.null.
            old_value_type = field.type.value_type
            new_value_type = (
                pa.large_string() if old_value_type == pa.null() else old_value_type
            )
            new_map[field.name] = pa.dictionary(
                new_index_type,
                new_value_type,
                field.type.ordered,
            )
        else:
            new_map[field.name] = field.type
    new_schema = pa.schema(new_map, metadata=arrow_table.schema.metadata)

    arrow_table = pa.Table.from_pandas(df, schema=new_schema)

    return arrow_table
