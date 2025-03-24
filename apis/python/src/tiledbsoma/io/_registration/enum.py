from __future__ import annotations

import numpy as np
import pandas as pd
import pyarrow as pa

from tiledbsoma import DataFrame


def _get_enumeration(df: DataFrame, column_name: str) -> pd.CategoricalDtype:
    import tiledb

    with tiledb.open(
        df.uri,
        timestamp=df.tiledb_timestamp_ms,
        mode="r",
        config=df.context.tiledb_config,
    ) as A:
        if column_name not in A.schema.attr_names:
            raise ValueError("Column not found in schema")

        enum_label = A.schema.attr(column_name).enum_label
        if enum_label is None:
            raise ValueError("Column is not an enumeration")

        enum = A.enum(enum_label)
        return pd.CategoricalDtype(categories=enum.values(), ordered=enum.ordered)


def _extend_enumeration(df: DataFrame, column_name: str, values: pa.Array) -> None:

    # first confirm we are a dictionary. If we have been decategorical-ized, we
    # will just be a plain-ol-array of value type.
    if not pa.types.is_dictionary(df.schema.field(column_name).type):
        return

    import tiledb

    # determine if we have any new enum values
    existing_dtype = _get_enumeration(df, column_name)
    new_enum_values = pd.Index(values).difference(existing_dtype.categories, sort=False)
    if len(new_enum_values) == 0:
        return

    with tiledb.open(df.uri, config=df.context.tiledb_config) as A:
        enum_label = A.schema.attr(column_name).enum_label
        # work around sc-64488
        enum = A.enum(enum_label)
        if np.issubdtype(enum.dtype, np.number):
            enum = enum.extend(new_enum_values.to_numpy())
        else:
            enum = enum.extend(new_enum_values.to_list())

    ctx = tiledb.Ctx(config=df.context.tiledb_config)
    se = tiledb.ArraySchemaEvolution(ctx=ctx)
    se.extend_enumeration(enum)
    se.array_evolve(df.uri)
