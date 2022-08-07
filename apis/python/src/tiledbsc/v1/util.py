import re
import time
from typing import TypeVar

import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


def tiledb_type_from_arrow_type(t: pa.DataType) -> type:
    """
    TODO
    """
    if t == pa.string():
        # pyarrow's to_pandas_dtype maps pa.string() to dtype object which
        # isn't acceptable to tiledb -- we must say str.
        return str
    else:
        # mypy says:
        # Returning Any from function declared to return "type"  [no-any-return]
        return t.to_pandas_dtype()  # type: ignore


def is_tiledb_creation_uri(uri: str) -> bool:
    return bool(re.match("^tiledb://.*s3://.*$", uri))


def get_start_stamp() -> float:
    """
    Returns information about start time of an event. Nominally float seconds since the epoch,
    but articulated here as being compatible with the format_elapsed function.
    """
    return time.time()


def format_elapsed(start_stamp: float, message: str) -> str:
    """
    Returns the message along with an elapsed-time indicator, with end time relative to start
    start from `get_start_stamp`. Used for annotating elapsed time of a task.
    """
    return "%s TIME %.3f seconds" % (message, time.time() - start_stamp)


# ----------------------------------------------------------------
def is_local_path(path: str) -> bool:
    """
    Returns information about start time of an event. Nominally float seconds since the epoch,
    but articulated here as being compatible with the format_elapsed function.
    """
    if path.startswith("file://"):
        return True
    if "://" in path:
        return False
    return True


# ----------------------------------------------------------------
def _to_tiledb_supported_dtype(dtype: np.dtype) -> np.dtype:
    """A handful of types are cast into the TileDB type system."""
    # TileDB has no float16 -- cast up to float32
    if dtype == np.dtype("float16"):
        return np.dtype("float32")
    return dtype


# ----------------------------------------------------------------

T = TypeVar("T", np.ndarray, pd.Series, pd.DataFrame, sp.spmatrix)


def _to_tiledb_supported_array_type(x: T) -> T:
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
            {k: _to_tiledb_supported_array_type(v) for k, v in x.items()}
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

    target_dtype = _to_tiledb_supported_dtype(x.dtype)
    return x if target_dtype == x.dtype else x.astype(target_dtype)


# ================================================================
# ================================================================
# ================================================================


def _ascii_to_unicode_dataframe_readback(df: pd.DataFrame) -> pd.DataFrame:
    """
    Implements the 'decode on read' part of our ASCII/Unicode logic
    """
    for k in df:
        dfk = df[k]
        if len(dfk) > 0 and type(dfk.iat[0]) == bytes:
            df[k] = dfk.map(lambda e: e.decode())
    return df


# ----------------------------------------------------------------
def _find_csr_chunk_size(
    mat: sp.csr_matrix,
    start_row_index: int,
    goal_chunk_nnz: int,
) -> int:
    """
    Given a CSR matrix and a start row index, returns the number of rows with cumulative NNZ as
    desired. Context is chunked-COO ingest of larger CSR matrices: if mat is say 8000x9000 but
    sparse, maybe we'll read rows 0:45 as one chunk and convert that to COO and ingest, then maybe
    rows 46:78 as a second chunk and convert that to COO and ingest, and so on.
    :param mat: The input CSR matrix.
    :param start_row_index: the row index at which to start a chunk.
    :param goal_chunk_nnz: Desired number of non-zero array entries for the chunk.
    """
    chunk_size = 1
    sum_nnz = 0
    for row_index in range(start_row_index, mat.shape[0]):
        sum_nnz += mat[row_index].nnz
        if sum_nnz > goal_chunk_nnz:
            break
        chunk_size += 1

    return chunk_size


# ----------------------------------------------------------------
# This function is very similar to _find_csr_chunk_size. The code is largely repeated, and this is
# intentional.  Here we err on the side of increased readability, at the expense of line-count.
def _find_csc_chunk_size(
    mat: sp.csc_matrix,
    start_col_index: int,
    goal_chunk_nnz: int,
) -> int:
    """
    Given a CSC matrix and a start column index, returns the number of columns with cumulative nnz as
    desired. Context is chunked-COO ingest of larger CSC matrices: if mat is say 8000x9000 but
    sparse, maybe we'll read columns 0:45 as one chunk and convert that to COO and ingest, then maybe
    columns 46:78 as a second chunk and convert that to COO and ingest, and so on.
    :param mat: The input CSC matrix.
    :param start_col_index: the column index at which to start a chunk.
    :param goal_chunk_nnz: Desired number of non-zero array entries for the chunk.
    """
    chunk_size = 1
    sum_nnz = 0
    for col_index in range(start_col_index, mat.shape[1]):
        sum_nnz += mat[:, col_index].nnz
        if sum_nnz > goal_chunk_nnz:
            break
        chunk_size += 1

    return chunk_size
