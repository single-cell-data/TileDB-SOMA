from typing import List, Union

import pandas as pd
import scipy.sparse as sp


def csr_from_tiledb_df(df: pd.DataFrame, num_rows: int, num_cols: int) -> sp.csr_matrix:
    """
    Given a tiledb dataframe, return a ``scipy.sparse.csr_matrx``.
    """
    return sp.csr_matrix(
        (df["soma_data"], (df["soma_dim_0"], df["soma_dim_1"])),
        shape=(num_rows, num_cols),
    )


def find_sparse_chunk_size(
    mat: Union[sp.csr_matrix, sp.csc_matrix],
    start_index: int,
    axis: int,
    goal_chunk_nnz: int,
) -> int:
    """
    Given a sparse matrix and a start index, return a step size, on the stride axis,
    which will achieve the cummulative nnz desired.

    :param mat: The input scipy.sparse matrix.
    :param start_index: the index at which to start a chunk.
    :param axis: the stride axis, across which to find a chunk.
    :param goal_chunk_nnz: Desired number of non-zero array entries for the chunk.
    """
    chunk_size = 1
    sum_nnz = 0
    coords: List[Union[slice, int]] = [slice(None), slice(None)]
    for index in range(start_index, mat.shape[axis]):
        coords[axis] = index
        sum_nnz += mat[tuple(coords)].nnz
        if sum_nnz > goal_chunk_nnz:
            break
        chunk_size += 1

    return chunk_size
