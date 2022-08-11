import time
from typing import TypeVar

import numpy as np
import pandas as pd
import scipy.sparse as sp

T = TypeVar("T", np.ndarray, pd.Series, pd.DataFrame, sp.spmatrix)

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


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


def is_local_path(path: str) -> bool:
    if path.startswith("file://"):
        return True
    if "://" in path:
        return False
    return True


def find_csr_chunk_size(
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


# This function is very similar to _find_csr_chunk_size. The code is largely repeated, and this is
# intentional.  Here we err on the side of increased readability, at the expense of line-count.
def find_csc_chunk_size(
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
