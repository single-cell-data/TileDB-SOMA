import functools
from typing import TypeVar, cast

import numpy as np
import pandas as pd
import pandas._typing as pdt
import scipy.sparse as sp
import tiledb

from .types import NPNDArray, PDSeries

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


_DT = TypeVar("_DT", bound=pdt.Dtype)


def _to_tiledb_supported_dtype(dtype: _DT) -> _DT:
    """A handful of types are cast into the TileDB type system."""
    # TileDB has no float16 -- cast up to float32
    return cast(_DT, np.dtype("float32")) if dtype == np.dtype("float16") else dtype


# Code-coverage note: the syntax in this file with @functools.singledispatch and @foo.register is
# likely to be unfamiliar to anyone but the most expert pythonistas. Suffice it to say that this
# function -- which appears to only ever raise unconditionally -- is invoked from the @foo.register
# functions via Python magic.
@functools.singledispatch
def to_tiledb_supported_array_type(x: object) -> object:
    """
    Converts datatypes unrepresentable by TileDB into datatypes it can represent.  E.g., categorical strings -> string.

    See also [https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html](https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html).

    Preferentially converts to the underlying primitive type, as TileDB does not support most complex types. NOTE: this does not support ``datetime64`` conversion.

    Categoricals are a special case. If the underlying categorical type is a primitive, convert to that. If the array contains NA/NaN (i.e. not in the category, code == -1), raise error unless it is a float or string.
    """
    raise TypeError(x.__class__)


@to_tiledb_supported_array_type.register
def _to_supported_dataframe(x: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame.from_dict(
        {k: to_tiledb_supported_array_type(v) for k, v in x.items()}
    )


@to_tiledb_supported_array_type.register
def _to_supported_series(x: PDSeries) -> PDSeries:
    if not pd.api.types.is_categorical_dtype(x):
        return _to_supported_base(x)

    categories = x.cat.categories
    cat_dtype = categories.dtype
    if cat_dtype.kind in ["f", "u", "i"]:
        if x.hasnans and cat_dtype.kind == "i":
            raise ValueError(
                "Categorical array contains NaN -- unable to convert to TileDB array."
            )

        return x.astype(_to_tiledb_supported_dtype(cat_dtype))

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


_MT = TypeVar("_MT", NPNDArray, sp.spmatrix, PDSeries)


def _to_supported_base(x: _MT) -> _MT:
    target_dtype = _to_tiledb_supported_dtype(x.dtype)
    return x if target_dtype == x.dtype else x.astype(target_dtype)


to_tiledb_supported_array_type.register(np.ndarray, _to_supported_base)
to_tiledb_supported_array_type.register(sp.spmatrix, _to_supported_base)


# ----------------------------------------------------------------
def list_fragments(array_uri: str) -> None:
    print(f"Listing fragments for array: '{array_uri}'")
    vfs = tiledb.VFS()

    fragments = []
    fi = tiledb.fragment.FragmentInfoList(array_uri=array_uri)

    for f in fi:
        f_dict = {
            "array_schema_name": f.array_schema_name,
            "num": f.num,
            "cell_num": f.cell_num,
            "size": vfs.dir_size(f.uri),
        }

        # parse nonempty domains into separate columns
        for d in range(len(f.nonempty_domain)):
            f_dict[f"d{d}"] = f.nonempty_domain[d]

        fragments.append(f_dict)

    frags_df = pd.DataFrame(fragments)
    print(frags_df)


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
