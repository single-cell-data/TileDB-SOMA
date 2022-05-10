import tiledb

import numpy
import scipy
import pandas as pd

import time
from typing import Optional

# ----------------------------------------------------------------
def get_start_stamp():
    """
    Returns information about start time of an event. Nominally float seconds since the epoch,
    but articulated here as being compatible with the format_elapsed function.
    """
    return time.time()

# ----------------------------------------------------------------
def format_elapsed(start_stamp, message: str):
    """
    Returns the message along with an elapsed-time indicator, with end time relative to start
    start from get_start_stamp. Used for annotating elapsed time of a task.
    """
    return "%s TIME %.3f" % (message, time.time() - start_stamp)

# ----------------------------------------------------------------
def find_csr_chunk_size(mat: scipy.sparse._csr.csr_matrix, permutation: list, start_row_index: int, goal_chunk_nnz: int):
    """
    Given a CSR matrix and a start row index, returns the number of rows with cumulative nnz as
    desired. Context is chunked-COO ingest of larger CSR matrices: if mat is say 8000x9000 but
    sparse, maybe we'll read rows 0:45 as one chunk and convert that to COO and ingest, then maybe
    rows 46:78 as a second chunk and convert that to COO and ingest, and so on.
    :param mat: The input CSR matrix.
    :param permutation: Cursor-indices to access the CSR matrix, so it will be traversed in sort order.
    :param start_row_index: the row index at which to start a chunk.
    :param goal_chunk_nnz: Desired number of non-zero array entries for the chunk. 
    """
    chunk_size = 1
    sum_nnz = 0
    for row_index in range(start_row_index, mat.shape[0]):
        sum_nnz += mat[permutation[row_index]].nnz
        if sum_nnz > goal_chunk_nnz:
            break
        chunk_size += 1

    return chunk_size

# ----------------------------------------------------------------
def get_sort_and_permutation(lst: list):
    """
    Sorts a list, returned the sorted list along with a permutation-index list which can be used for
    cursored access to data which was indexed by the unsorted list. Nominally for chunking of CSR
    matrices into TileDB which needs sorted string dimension-values for efficient fragmentation.
    """
    # Example input: x=['E','A','C','D','B']

    # e.g. [('E', 0), ('A', 1), ('C', 2), ('D', 3), ('B', 4)]
    lst_and_indices = [(e, i) for i,e in enumerate(lst)]

    # e.g. [('A', 1), ('B', 4), ('C', 2), ('D', 3), ('E', 0)]
    lst_and_indices.sort(key=lambda pair: pair[0])

    # e.g. ['A' 'B' 'C' 'D' 'E']
    # and  [1, 4, 2, 3, 0]
    lst_sorted  = [e for e,i in lst_and_indices]
    permutation = [i for e,i in lst_and_indices]
    return (lst_sorted, permutation)

# ----------------------------------------------------------------
def _to_tiledb_supported_array_type(x):
    """
    Converts datatypes unrepresentable by TileDB into datatypes it can represent.
    Eg, categorical strings -> string; bool -> uint8, etc.

    See also https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html

    Preferentially converts to the underlying primitive type, as TileDB does not
    support most complex types. NOTE: this does not support `datetime64` conversion.

    Categoricals are a special case. If the underlying categorical type is a
    primitive, convert to that. If the array contains NA/NaN (i.e. not in the
    category, code == -1), raise error unless it is a float or string.
    """

    def _to_tiledb_supported_dtype(dtype):
        """A handful of types are cast into the TileDB type system."""
        # TileDB has no bool type -- instead cast to uint8
        if dtype == numpy.dtype('bool'):
            return numpy.dtype('uint8')

        # TileDB has no float16 -- cast up to float32
        if dtype == numpy.dtype('float16'):
            return numpy.dtype('float32')

        return dtype

    # If a Pandas categorical, use the type of the underlying category.
    # If the array contains NaN/NA, and the primitive is unable to represent
    # a reasonable facsimile, i.e. not string or float, raise.
    if pd.api.types.is_categorical_dtype(x.dtype):
        categories = x.cat.categories
        cat_dtype = categories.dtype
        if cat_dtype.kind in ['f', 'u', 'i']:
            if x.hasnans and cat_dtype.kind == 'i':
                raise ValueError("Categorical array contains NaN -- unable to convert to TileDB array.")

            return x.astype(_to_tiledb_supported_dtype(cat_dtype))

        # Into the weirdness. See if Pandas can help with edge cases.
        inferred = pd.api.types.infer_dtype(categories)
        if inferred == "boolean":
            if x.hasnans:
                raise ValueError("Categorical array contains NaN -- unable to convert to TileDB array.")
            return x.astype('uint8')

        if inferred == "string":
            return x.astype(str)

        if inferred == "bytes":
            if x.hasnans:
                raise ValueError("Categorical array contains NaN -- unable to convert to TileDB array.")
            return x.astype(bytes)

        return x.astype('O')

    target_dtype = _to_tiledb_supported_dtype(x.dtype)
    return x if target_dtype == x.dtype else x.astype(target_dtype)

# ----------------------------------------------------------------
def is_numpyable_object(obj):
    """
    Checks if the argument is of numpy type, nominally as a gate before a call to
    numpyable_object_to_tiledb_array. Additionally, for the benefit of unit test, supports
    native Python types.
    """

    # Unit test
    if _is_numpyable_scalar(obj):
        return True
    if isinstance(obj, list):
        return True

    # anndata .h5ad contents
    return 'numpy' in str(type(obj))

# ----------------------------------------------------------------
def _is_numpyable_scalar(obj):
    """
    Check if the object is a scalar we'll be able to wrap in a list and then
    turn that into a 1D array for tiledb.numpy. Nominally for unit-test data.
    """
    if isinstance(obj, int):
        return True
    if isinstance(obj, float):
        return True
    if isinstance(obj, str):
        return True
    return False

# ----------------------------------------------------------------
def numpyable_object_to_tiledb_array(obj, uri: str, ctx: Optional[tiledb.Ctx] = None):
    """
    Nominally for ingest of `uns` nested data from anndata objects. Handles scalar or array values
    -- the former, by wrapping in a 1D array. Maps to TileDB / tiledb.from_numpy storage semantics,
    including UTF-8 handling. Supports dtypes like
    """

    if isinstance(obj, numpy.ndarray):
        obj = _to_tiledb_supported_array_type(obj)
        _write_numpy_ndarray_to_tiledb_array(arr=obj, uri=uri, ctx=ctx)

    elif isinstance(obj, list):
        arr = numpy.asarray(obj)
        _write_numpy_ndarray_to_tiledb_array(arr, uri, ctx)

    elif isinstance(obj, numpy.str_):
        # Needs explicit cast from numpy.str_ to str for tiledb.from_numpy
        arr = numpy.asarray([obj]).astype(str)
        _write_numpy_ndarray_to_tiledb_array(arr, uri, ctx)

    else:
        arr = numpy.asarray([obj])
        arr = _to_tiledb_supported_array_type(arr)
        _write_numpy_ndarray_to_tiledb_array(arr, uri, ctx)

# ----------------------------------------------------------------
def _write_numpy_ndarray_to_tiledb_array(arr: numpy.ndarray, uri: str, ctx: Optional[tiledb.Ctx] = None):
    """
    Writes a numpy.ndarray to a TileDB array, nominally for ingest of `uns` nested data from anndata
    objects. Mostly tiledb.from_numpy, but with some necessary handling for data with UTF-8 values.
    """

    if 'numpy' in str(type(arr)) and str(arr.dtype).startswith('<U'):
        # Note arr.astype('str') does not lead to a successfuly tiledb.from_numpy.
        arr = numpy.array(arr, dtype='O')

    tiledb.from_numpy(uri=uri, array=arr, ctx=ctx)
