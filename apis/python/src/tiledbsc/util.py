import tiledb

import anndata as ad

import numpy
import scipy
import pandas as pd

import time
from typing import Optional, List

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
    return "%s TIME %.3f seconds" % (message, time.time() - start_stamp)


# ----------------------------------------------------------------
def _find_csr_chunk_size(
    mat: scipy.sparse._csr.csr_matrix,
    permutation: list,
    start_row_index: int,
    goal_chunk_nnz: int,
):
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
# This function is very similar to _find_csr_chunk_size. The code is largely repeated, and this is
# intentional.  Here we err on the side of increased readability, at the expense of line-count.
def _find_csc_chunk_size(
    mat: scipy.sparse._csc.csc_matrix,
    permutation: list,
    start_col_index: int,
    goal_chunk_nnz: int,
):
    """
    Given a CSC matrix and a start column index, returns the number of columns with cumulative nnz as
    desired. Context is chunked-COO ingest of larger CSC matrices: if mat is say 8000x9000 but
    sparse, maybe we'll read columns 0:45 as one chunk and convert that to COO and ingest, then maybe
    columns 46:78 as a second chunk and convert that to COO and ingest, and so on.
    :param mat: The input CSC matrix.
    :param permutation: Cursor-indices to access the CSC matrix, so it will be traversed in sort order.
    :param start_col_index: the column index at which to start a chunk.
    :param goal_chunk_nnz: Desired number of non-zero array entries for the chunk.
    """
    chunk_size = 1
    sum_nnz = 0
    for col_index in range(start_col_index, mat.shape[1]):
        sum_nnz += mat[:, permutation[col_index]].nnz
        if sum_nnz > goal_chunk_nnz:
            break
        chunk_size += 1

    return chunk_size


# ----------------------------------------------------------------
def _get_sort_and_permutation(lst: list):
    """
    Sorts a list, returned the sorted list along with a permutation-index list which can be used for
    cursored access to data which was indexed by the unsorted list. Nominally for chunking of CSR
    matrices into TileDB which needs sorted string dimension-values for efficient fragmentation.
    """
    # Example input: x=['E','A','C','D','B']

    # e.g. [('E', 0), ('A', 1), ('C', 2), ('D', 3), ('B', 4)]
    lst_and_indices = [(e, i) for i, e in enumerate(lst)]

    # e.g. [('A', 1), ('B', 4), ('C', 2), ('D', 3), ('E', 0)]
    lst_and_indices.sort(key=lambda pair: pair[0])

    # e.g. ['A' 'B' 'C' 'D' 'E']
    # and  [1, 4, 2, 3, 0]
    lst_sorted = [e for e, i in lst_and_indices]
    permutation = [i for e, i in lst_and_indices]
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
        if dtype == numpy.dtype("bool"):
            return numpy.dtype("uint8")

        # TileDB has no float16 -- cast up to float32
        if dtype == numpy.dtype("float16"):
            return numpy.dtype("float32")

        return dtype

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
            return x.astype("uint8")

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


# ----------------------------------------------------------------
def _X_and_ids_to_coo(
    Xdf: pd.DataFrame,
    row_dim_name: str,
    col_dim_name: str,
    attr_name: str,
    row_labels,
    col_labels,
) -> scipy.sparse.csr_matrix:
    """
    This is needed when we read a TileDB X.df[:]. Since TileDB X is sparse 2D string-dimensioned,
    the return value of which is a dict with three columns -- obs_id, var_id, and value. For
    conversion to anndata, we need make a sparse COO/IJV-format array where the indices are
    not strings but ints, matching the obs and var labels.
    """

    # Now we need to convert from TileDB's string indices to CSR integer indices.
    # Make a dict from string dimension values to integer indices.
    #
    # Example: suppose the sparse matrix looks like:
    #
    #     S T U V
    #   A 4 . . 3
    #   B: 5 . 6 .
    #   C . 1 . 2
    #   D 8 7 . .
    #
    # The return value from the X[:] query is (obs_id,var_id,value) triples like
    #
    #   A,S,4 A,V,3 B,S,5 B,U,6 C,V,2 C,T,1 D,S,8 D,T,7
    #
    # whereas scipy csr is going to want
    #
    #   0,0,4 0,3,3 1,0,5 1,2,6 2,3,2 2,1,1 3,0,8 3,1,7
    #
    # In order to accomplish this, we need to map ['A','B','C','D'] to [0,1,2,3] via {'A':0,
    # 'B':1, 'C':2, 'D':3} and similarly for the other dimension.
    row_labels_to_indices = dict(zip(row_labels, [i for i, e in enumerate(row_labels)]))
    col_labels_to_indices = dict(zip(col_labels, [i for i, e in enumerate(col_labels)]))

    # Apply the map.
    obs_indices = [row_labels_to_indices[row_label] for row_label in Xdf[row_dim_name]]
    var_indices = [col_labels_to_indices[col_label] for col_label in Xdf[col_dim_name]]

    return scipy.sparse.csr_matrix(
        (list(Xdf[attr_name]), (list(obs_indices), list(var_indices)))
    )


# ================================================================
class ETATracker:
    """
    Computes estimated time to completion for chunked writes.
    """

    percents: List[float]
    cumulative_seconds: List[float]

    def __init__(self):
        self.chunk_percents = []
        self.cumulative_seconds = []

    def ingest_and_predict(self, chunk_percent: float, chunk_seconds: float) -> str:
        """
        Updates from most recent chunk percent-done and chunk completion-seconds, then does a linear
        regression on all chunks done so far and estimates time to completion.
        :param chunk_percent: a percent done like 6.1 or 10.3.
        :param chunk_seconds: number of seconds it took to do the current chunk operation.
        """
        self._ingest(chunk_percent, chunk_seconds)
        eta_seconds = self._predict()
        eta_format = self._format_seconds(eta_seconds)
        return eta_format

    def _ingest(self, chunk_percent: float, chunk_seconds: float) -> None:
        """
        Takes the current percent done like 10.3 and current chunk seconds like 58.4 and grows an
        array of percent-dones and cumulative seconds. This means self.chunk_percents is a list of
        all the chunk_percent arguments from calling _ingest, while each self.chunk_seconds slot is
        the sum of all previous chunk_seconds arguments from calling _ingest.
        """
        if len(self.chunk_percents) == 0:
            self.chunk_percents = [chunk_percent]
            self.cumulative_seconds = [chunk_seconds]
        else:
            self.chunk_percents.append(chunk_percent)
            self.cumulative_seconds.append(self.cumulative_seconds[-1] + chunk_seconds)

    def _predict(self) -> float:
        """
        Does a linear regression on all chunks done so far and estimates time to completion.
        Returns ETA seconds as a number.
        """
        # Linear regression where x is cumulative seconds and y is percent done.
        x = numpy.array(self.cumulative_seconds)
        y = numpy.array(self.chunk_percents)
        A = numpy.vstack([x, numpy.ones(len(x))]).T
        m, b = numpy.linalg.lstsq(A, y, rcond=None)[0]
        # Solve for x where y == 100
        done_cumu_seconds = (100.0 - b) / m

        return done_cumu_seconds - self.cumulative_seconds[-1]

    def _format_seconds(self, seconds: float) -> str:
        """
        Formats the ETA seconds as a compact, human-readable string.
        """
        if seconds >= 86400:
            return "%.2f days" % (seconds / 86400)
        elif seconds >= 3600:
            return "%.2f hours" % (seconds / 3600)
        elif seconds >= 60:
            return "%.2f minutes" % (seconds / 60)
        else:
            return "%.2f seconds" % (seconds)

    def __str__(self):
        return str(self.chunk_percents) + " " + str(self.cumulative_seconds)

    def __repr__(self):
        return self.__str__()
