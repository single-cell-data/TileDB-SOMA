import re
from typing import Optional, TypeVar

import numpy as np
import pandas as pd
import scipy.sparse as sp

T = TypeVar("T", np.ndarray, pd.Series, pd.DataFrame, sp.spmatrix)

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


def is_tiledb_creation_uri(uri: str) -> bool:
    return bool(re.match("^tiledb://.*s3://.*$", uri))


def tiledb_result_order_from_soma_result_order_non_indexed(
    soma_result_order: Optional[str],
) -> Optional[str]:
    """
    Maps SOMA-spec `result_order` syntax to TileDB-specific syntax, for non-indexed dataframes.
    """
    # :param order: 'C', 'F', or 'G' (row-major, col-major, tiledb global order)
    if soma_result_order is None:
        return None  # use tiledb default
    elif soma_result_order == "rowid-ordered":
        return "C"
    elif soma_result_order == "unordered":
        return "U"
    else:
        raise Exception(f'result-order value unrecognized: "{soma_result_order}"')


def tiledb_result_order_from_soma_result_order_indexed(
    soma_result_order: Optional[str],
) -> Optional[str]:
    """
    Maps SOMA-spec `result_order` syntax to TileDB-specific syntax, for indexed dataframes.
    """
    # :param order: 'C', 'F', or 'G' (row-major, col-major, tiledb global order)
    if soma_result_order is None:
        return None  # use tiledb default
    elif soma_result_order == "row-major":
        return "C"
    elif soma_result_order == "col-major":
        return "F"
    elif soma_result_order == "unordered":
        return "U"
    else:
        raise Exception(f'result-order value unrecognized: "{soma_result_order}"')


def to_tiledb_supported_dtype(dtype: np.dtype) -> np.dtype:
    """A handful of types are cast into the TileDB type system."""
    # TileDB has no float16 -- cast up to float32
    if dtype == np.dtype("float16"):
        return np.dtype("float32")
    return dtype


def to_tiledb_supported_array_type(x: T) -> T:
    """
    Converts datatypes unrepresentable by TileDB into datatypes it can represent.
    E.g., categorical strings -> string.

    See also https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html

    Preferentially converts to the underlying primitive type, as TileDB does not
    support most complex types. NOTE: this does not support `datetime64` conversion.

    Categoricals are a special case. If the underlying categorical type is a
    primitive, convert to that. If the array contains NA/NaN (i.e. not in the
    category, code == -1), raise error unless it is a float or string.
    """
    if isinstance(x, pd.DataFrame):
        return pd.DataFrame.from_dict(
            {k: to_tiledb_supported_array_type(v) for k, v in x.items()}
        )

    # If a Pandas categorical, use the type of the underlying category.
    # If the array contains NaN/NA, and the primitive is unable to represent
    # a reasonable facsimile, i.e. not string or float, raise.
    if pd.api.types.is_categorical_dtype(x.dtype):
        categories = x.cat.categories
        cat_dtype = categories.dtype
        if cat_dtype.kind in ["f", "u", "i"]:
            if x.hasnans and cat_dtype.kind == "i":
                raise ValueError(
                    "Categorical array contains NaN -- unable to convert to TileDB array."
                )

            return x.astype(to_tiledb_supported_dtype(cat_dtype))

        # Into the weirdness. See if Pandas can help with edge cases.
        inferred = pd.api.types.infer_dtype(categories)
        if inferred == "boolean":
            if x.hasnans:
                raise ValueError(
                    "Categorical array contains NaN -- unable to convert to TileDB array."
                )
            return x.astype(bool)

        if inferred == "string":
            return x.astype(str)

        if inferred == "bytes":
            if x.hasnans:
                raise ValueError(
                    "Categorical array contains NaN -- unable to convert to TileDB array."
                )
            return x.astype(bytes)

        return x.astype("O")

    target_dtype = to_tiledb_supported_dtype(x.dtype)
    return x if target_dtype == x.dtype else x.astype(target_dtype)
