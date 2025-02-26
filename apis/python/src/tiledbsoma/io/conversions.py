# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Conversion utility methods.
"""

from __future__ import annotations

from typing import Dict, TypeVar, Union, cast

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


# Metadata indicating a SOMA DataFrame's original index column name, serialized as a
# JSON string or `null`. SOMA DataFrames are always given a `soma_joinid` index, but
# we want to be able to outgest a `pd.DataFrame` that is identical to the one we
# ingested, so we store an "original index name" in the DataFrame's metadata.
OriginalIndexMetadata = Union[None, str]


def _string_dict_from_arrow_schema(schema: pa.Schema) -> Dict[str, str]:
    """Converts an Arrow schema to a string/string dict.

    This is easier on the eyes, easier to convert from/to JSON for distributed logging,
    and easier to do del-key on.
    """

    _EQUIVALENCES = {
        "large_string": "string",
        "large_binary": "binary",
    }

    def _stringify_type(arrow_type: pa.DataType) -> str:
        """Turns an Arrow data type into a string.

        Note: Arrow string and large_string must map to TileDB string, which is
        large-only. Similarly for Arrow binary and large_binary.
        """
        if pa.types.is_dictionary(arrow_type):
            arrow_type = arrow_type.index_type
        str_type = str(arrow_type)
        return _EQUIVALENCES.get(str_type, str_type)

    # Stringify types skipping the soma_joinid field (it is specific to SOMA data
    # and does not exist in AnnData/H5AD).
    arrow_columns = {
        name: _stringify_type(schema.field(name).type)
        for name in schema.names
        if name != "soma_joinid"
    }

    return arrow_columns


def _prepare_df_for_ingest(df: pd.DataFrame, id_column_name: str | None) -> str | None:
    """Prepare a `pd.DataFrame` for persisting as a SOMA DataFrame.

    Demote its index to a column (to make way for a required `soma_joinid` index), and
    compute and return metadata for restoring the index column and name later (on
    outgest).

    If `df.index` has a name (and it's not "index", which is taken to be a default/unset
    value):
      - `df.index.name` takes precedence over the `id_column_name` arg: the index
      will be reset to an eponymous column.
      - That original `df.index.name` will be logged as `OriginalIndexMetadata` (for
      promotion back to index on outgest).

    In this case, the overall round trip is basically just:
    - `reset_index` on ingest (demote index to eponymous column).
    - `set_index` on outgest (restore column to index, with its original name).

    Otherwise (index name is `None` or "index"):
    - A fallback name (`id_column_name` if provided, "index" otherwise) is used for the
    column that the index becomes.
    - The returned `OriginalIndexMetadata` will be `None`.

    There are several edge cases (detailed below and in
    `test_dataframe_io_roundtrips.py` and
    https://github.com/single-cell-data/TileDB-SOMA/issues/2829) where the index, its
    name, or a specific column are not restored properly on outgest. For now, all such
    behavior is preserved, for backwards compatibility, but we should look into ways
    of improving these "round-trip mutation" cases.

    See https://github.com/single-cell-data/TileDB-SOMA/issues/2829 for more info.
    """
    use_existing_index = df.index.name is not None and df.index.name != "index"

    original_index_name = None
    if use_existing_index:
        original_index_name = df.index.name

    df.reset_index(inplace=True)
    if id_column_name is not None:
        if id_column_name in df:
            if "index" in df:
                # The assumption here is that the column named "index" was previously
                # an unnamed `df.index`, and `id_column_name` was already a column
                # (per the grandparent `if` above). In this case, we drop the
                # original unnamed `df.index`.
                # TODO: This prevents outgesting the same DataFrame we ingested. We
                # should fix it; see
                #  https://github.com/single-cell-data/TileDB-SOMA/issues/2829.
                #
                # Also note: if the DataFrame already had columns named "index" and
                # `id_column_name`, the original `df.index` will have been "reset" to
                # a column named `level_0`, and we end up just dropping the column
                # named "index" here.
                #
                # Another version of this occurs when the original DataFrame has
                # `df.index.name == id_column_name` and a column named "index". In this
                # case, the index will have been "reset" to a column named
                # `id_column_name` above, which then satisfies the grendparent `if`'s
                # predicate, and causes us to drop the column named "index" here.
                df.drop(columns=["index"], inplace=True)
        else:
            # If `id_column_name` was passed, and is not already a column in the
            # DataFrame, we assume the original index was "reset" to a column named
            # "index" (by `reset_index` above), and we rename that column to
            # `id_column_name`, so that `id_column_name` matches the name of a column
            # representing the original
            # DataFrame's index.
            #
            # NOTE: the assumption above can break in a few ways:
            # 1. The original DataFrame index has a name other than "index" or
            #    `id_column_name`…
            #    a. and there is a column named "index" ⇒ that column will be renamed to
            #       `id_column_name`
            #    b. and there is no column named "index" ⇒ the rename below is a no-op
            #       (outgest currently restores the original DataFrame in this case)
            # 2. The original DataFrame has a column named "index":
            #    - That column will become `df.index` on outgest, and acquire the
            #      original `df.index.name` as its name.
            #    - The original index will end up as a column, on outgest:
            #      - If it had a name, the column will have that name.
            #      - Otherwise, it will end up as a column named e.g. `level_0` (or
            #        `level_1`, if a column named `level_0` already exists, etc.)
            #
            # See https://github.com/single-cell-data/TileDB-SOMA/issues/2829 for more info.
            df.rename(columns={"index": id_column_name}, inplace=True)

    return original_index_name


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


def df_to_arrow_table(df: pd.DataFrame) -> pa.Table:
    """
    Handle special cases where pa.Table.from_pandas is not sufficient.
    """
    nullable_fields = set()
    # Not for name, col in df.items() since we need df[k] on the left-hand sides
    for key in df:
        # Make attributes nullable. Context:
        # * this is _solely_ for use of tiledbsoma.io
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


def df_to_arrow_schema(df: pd.DataFrame, default_index_name: str) -> pa.Schema:
    """Creates the arrow schema from a pandas dataframe.

    This function does not mutate the input ``pandas.DataFrame``.
    """
    df = df.head(1).copy()  # since reset_index can be expensive on full data
    _prepare_df_for_ingest(df, default_index_name)
    arrow_table = df_to_arrow_table(df)
    arrow_schema = arrow_table.schema.remove_metadata()
    return arrow_schema
