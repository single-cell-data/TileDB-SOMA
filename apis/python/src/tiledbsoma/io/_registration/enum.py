from __future__ import annotations

from typing import Sequence

import pandas as pd
import pyarrow as pa

from tiledbsoma import DataFrame


def get_enumerations(
    df: DataFrame, column_names: Sequence[str]
) -> dict[str, pd.CategoricalDtype]:
    """Look up enum info in schema, and return as a Pandas CategoricalDType. This
    is a convenience wrapper around ``DataFrame.get_enumeration_values``, for use
    in the registration module.
    """
    # skip columns which are not of type dictionary
    column_names = [
        c for c in column_names if pa.types.is_dictionary(df.schema.field(c).type)
    ]
    return {
        k: pd.CategoricalDtype(categories=v, ordered=df.schema.field(k).type.ordered)
        for k, v in df.get_enumeration_values(column_names).items()
    }


def extend_enumerations(df: DataFrame, columns: dict[str, pd.CategoricalDtype]) -> None:
    """Extend enumerations as needed, starting with a CategoricalDType for each
    cat/enum/dict column. A convenience wrapper around ``DataFrame.extend_enumeration_values``,
    for use in the registration module.

    DataFrame must be open for write.
    """
    current_enums = get_enumerations(df, list(columns.keys()))
    columns_to_extend = {}
    for column_name, cat_dtype in columns.items():

        # first confirm this is a dictionary. If it has been decategorical-ized, i.e.,
        # are an array of the value type, don't extend.
        if column_name not in current_enums:
            assert not pa.types.is_dictionary(df.schema.field(column_name).type)
            continue

        # determine if we have any new enum values in this column
        existing_dtype = current_enums[column_name]
        new_enum_values = pd.Index(cat_dtype.categories).difference(
            existing_dtype.categories, sort=False
        )
        if len(new_enum_values) == 0:
            continue

        # if there are new values, extend the array schema enum
        new_enum_values = pa.array(new_enum_values.to_numpy())
        columns_to_extend[column_name] = new_enum_values

    # and evolve the schema
    df.extend_enumeration_values(columns_to_extend, deduplicate=False)
