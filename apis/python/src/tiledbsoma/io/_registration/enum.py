from __future__ import annotations

import pandas as pd
import pyarrow as pa

from tiledbsoma import DataFrame


def _get_enumeration(df: DataFrame, column_name: str) -> pd.CategoricalDtype:
    values = df.get_enumeration_values((column_name,))[column_name]
    ordered = df.schema.field(column_name).type.ordered
    return pd.CategoricalDtype(categories=values, ordered=ordered)


def _extend_enumeration(df: DataFrame, column_name: str, values: pa.Array) -> None:

    # first confirm we are a dictionary. If we have been decategorical-ized, we
    # will just be an array of value type.
    if not pa.types.is_dictionary(df.schema.field(column_name).type):
        return

    # determine if we have any new enum values
    existing_dtype = _get_enumeration(df, column_name)
    new_enum_values = pd.Index(values).difference(existing_dtype.categories, sort=False)
    if len(new_enum_values) == 0:
        return

    # if there are new values, extend the array schema enum
    new_enum_values = pa.array(new_enum_values.to_numpy())
    if df.mode == "w":
        df.extend_enumeration_values({column_name: new_enum_values}, deduplicate=False)
    else:
        with df.reopen(mode="w") as wdf:
            wdf.extend_enumeration_values(
                {column_name: new_enum_values}, deduplicate=False
            )
