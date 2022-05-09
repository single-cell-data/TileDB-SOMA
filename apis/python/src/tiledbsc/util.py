import numpy as np
import scipy
import pandas as pd
import time

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
    primitive, convert to that. If the array contains NA/nan (ie, not in the
    category, code == -1), raise error unless it is a float or string.
    """

    def _to_tiledb_supported_dtype(dtype):
        """A handful of types are cast into the TileDB type system."""
        # TileDB has no bool type - instead cast to uint8
        if dtype == np.dtype('bool'):
            return np.dtype('uint8')

        # TileDB has no float16 - cast up to float32
        if dtype == np.dtype('float16'):
            return np.dtype('float32')

        return dtype

    # if a Pandas categorical, use the type of the underlying category.
    # If the array contains nan/NA, and the primitive is unable to represent
    # a reasonable facsimile, ie, not string or float, raise.
    if pd.api.types.is_categorical_dtype(x.dtype):
        categories = x.cat.categories
        cat_dtype = categories.dtype
        if cat_dtype.kind in ['f', 'u', 'i']:
            if x.hasnans and cat_dtype.kind == 'i':
                raise ValueError("Categorical array contains nan - unable to convert to TileDB array.")

            return x.astype(_to_tiledb_supported_dtype(cat_dtype))

        # Into the weirdness. See if Pandas can help with edge cases
        inferred = pd.api.types.infer_dtype(categories)
        if inferred == "boolean":
            if x.hasnans:
                raise ValueError("Categorical array contains nan - unable to convert to TileDB array.")
            return x.astype('uint8')

        if inferred == "string":
            return x.astype(str)

        if inferred == "bytes":
            if x.hasnans:
                raise ValueError("Categorical array contains nan - unable to convert to TileDB array.")
            return x.astype(bytes)

        return x.astype('O')

    target_dtype = _to_tiledb_supported_dtype(x.dtype)
    return x if target_dtype == x.dtype else x.astype(target_dtype)
