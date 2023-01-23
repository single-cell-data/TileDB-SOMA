from typing import TypeVar, cast

import numpy as np
import pandas._typing as pdt
import scipy.sparse as sp
import tiledb
from pandas.api.types import infer_dtype, is_categorical_dtype

from .types import NPNDArray, PDSeries

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


_DT = TypeVar("_DT", bound=pdt.Dtype)
_MT = TypeVar("_MT", NPNDArray, sp.spmatrix, PDSeries)
_str_to_type = {"boolean": bool, "string": str, "bytes": bytes}


def _to_tiledb_supported_dtype(dtype: _DT) -> _DT:
    """A handful of types are cast into the TileDB type system."""
    # TileDB has no float16 -- cast up to float32
    return cast(_DT, np.dtype("float32")) if dtype == np.dtype("float16") else dtype


def to_tiledb_supported_array_type(x: _MT) -> _MT:
    """
    Converts datatypes unrepresentable by TileDB into datatypes it can represent.
    E.g., categorical strings -> string.

    See also [https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html](https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html).

    Preferentially converts to the underlying primitive type, as TileDB does not support
    most complex types. NOTE: this does not support ``datetime64`` conversion.

    Categoricals are a special case. If the underlying categorical type is a primitive,
    convert to that. If the array contains NA/NaN (i.e. not in the category, code == -1),
    raise error unless it is a float or string.
    """
    if isinstance(x, (np.ndarray, sp.spmatrix)) or not is_categorical_dtype(x):
        target_dtype = _to_tiledb_supported_dtype(x.dtype)
        return x if target_dtype == x.dtype else x.astype(target_dtype)

    categories = x.cat.categories
    cat_dtype = categories.dtype
    if cat_dtype.kind in ("f", "u", "i"):
        if x.hasnans and cat_dtype.kind == "i":
            raise ValueError(
                "Categorical array contains NaN -- unable to convert to TileDB array."
            )
        target_dtype = _to_tiledb_supported_dtype(cat_dtype)
    else:
        # Into the weirdness. See if Pandas can help with edge cases.
        inferred = infer_dtype(categories)
        if x.hasnans and inferred in ("boolean", "bytes"):
            raise ValueError(
                "Categorical array contains NaN -- unable to convert to TileDB array."
            )
        target_dtype = np.dtype(_str_to_type.get(inferred, object))

    return x.astype(target_dtype)


def is_does_not_exist_error(e: tiledb.TileDBError) -> bool:
    """ "
    Given a TileDBError, return true if it indicates the object does not exist.

    Example
    -------

    try:
        with tiledb.open(uri):
            ...
    except tiledb.TileDBError as e:
        if is_does_not_exist_error(e):
            ...
        raise e
    """
    stre = str(e)
    # Local-disk/S3 does-not-exist exceptions say 'Group does not exist'; TileDB Cloud
    # does-not-exist exceptions are worded less clearly.
    if "does not exist" in stre or "HTTP code 401" in stre:
        return True

    return False


def is_duplicate_group_key_error(e: tiledb.TileDBError) -> bool:
    """
    Given a TileDBError, return try if it indicates a duplicate member
    add request in a tiledb.Group.
    """
    stre = str(e)
    if "member already exists in group" in stre:
        return True

    return False
