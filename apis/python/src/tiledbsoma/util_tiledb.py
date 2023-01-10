import re
from typing import Dict, List, Optional, Sequence, Tuple, TypeVar, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp
import tiledb

from .types import ResultOrder

T = TypeVar("T", np.ndarray, pd.Series, pd.DataFrame, sp.spmatrix)

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


def is_tiledb_creation_uri(uri: str) -> bool:
    return bool(re.match("^tiledb://.*s3://.*$", uri))


def tiledb_result_order_from_soma_result_order(
    soma_result_order: Optional[ResultOrder], accept: List[ResultOrder]
) -> Optional[str]:
    """
    Given a ResultOrder, return a TileDB result order.  Raise an error if
    the ``soma_result_order`` is not present in the acceptable values, as
    defined by ``accept``.
    """
    OrderMap: Dict[ResultOrder, str] = {
        "column-major": "F",
        "row-major": "C",
        "auto": "U",
    }
    if not soma_result_order:
        return None
    if soma_result_order not in accept or soma_result_order not in OrderMap:
        raise ValueError("result_order unsupported supported.")
    return OrderMap[soma_result_order]


def to_tiledb_supported_dtype(dtype: np.dtype) -> np.dtype:
    """A handful of types are cast into the TileDB type system."""
    # TileDB has no float16 -- cast up to float32
    if dtype == np.dtype("float16"):
        return np.dtype("float32")
    return dtype


def to_tiledb_supported_array_type(x: T) -> T:
    """
    Converts datatypes unrepresentable by TileDB into datatypes it can represent.  E.g., categorical strings -> string.

    See also [https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html](https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html).

    Preferentially converts to the underlying primitive type, as TileDB does not support most complex types. NOTE: this does not support ``datetime64`` conversion.

    Categoricals are a special case. If the underlying categorical type is a primitive, convert to that. If the array contains NA/NaN (i.e. not in the category, code == -1), raise error unless it is a float or string.
    """
    if isinstance(x, pd.DataFrame):
        return pd.DataFrame.from_dict(
            {k: to_tiledb_supported_array_type(v) for k, v in x.items()}
        )

    # If a Pandas categorical, use the type of the underlying category.
    # If the array contains NaN/NA, and the primitive is unable to represent
    # a reasonable facsimile, i.e. not string or float, raise.
    if pd.api.types.is_categorical_dtype(x.dtype):
        assert isinstance(x, pd.Series)
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


def split_column_names(
    array_schema: tiledb.ArraySchema, column_names: Optional[Sequence[str]]
) -> Tuple[Union[Sequence[str], None], Union[Sequence[str], None]]:
    """
    Given a tiledb ArraySchema and a list of dim or attr names, split
    them into a tuple of (dim_names, attr_names).

    This helper is used to turn the SOMA ``column_names`` parameter into a
    form that can be natively used by tiledb.Array.query, which requires
    that the list of names be separated into ``dims`` and ``attrs``.

    Parameters
    ----------
    array_schema : tiledb.Array
        An array schema which will be used to determine whether a
        column name is a dim or attr.
    column_names : Optional[Sequence[str]]
        List of column names to split into ``dim`` and ``attr`` names.

    Returns
    -------
    Tuple[Union[Sequence[str], None], Union[Sequence[str], None]]
        If column_names is ``None``, the tuple ``(None, None)`` will be returned.
        Otherwise, returns a tuple of (dim_names, attr_names), with any unknown
        names, i.e., not present in the array schema, ignored (dropped).
    """
    if column_names is None:
        return (None, None)

    dim_names = [
        array_schema.domain.dim(i).name for i in range(array_schema.domain.ndim)
    ]
    attr_names = [array_schema.attr(i).name for i in range(array_schema.nattr)]
    return (
        [c for c in column_names if c in dim_names],
        [c for c in column_names if c in attr_names],
    )


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
