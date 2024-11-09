"""Test the fastercsx pybind11 C++ layer."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from scipy import sparse

import tiledbsoma.pytiledbsoma as clib
import tiledbsoma.pytiledbsoma.fastercsx as fastercsx


@pytest.fixture
def concurrency() -> int | None:
    return None


@pytest.fixture
def context(concurrency: int | None) -> clib.SOMAContext:
    if concurrency is None:
        return clib.SOMAContext()
    else:
        return clib.SOMAContext({"soma.compute_concurrency_level": f"{concurrency}"})


@pytest.mark.parametrize(
    "shape", [(0, 0), (1, 0), (0, 1), (10, 100), (100, 10), (9972, 1), (1, 9972)]
)
@pytest.mark.parametrize(
    "value_dtype",
    [
        np.float32,
        np.float64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.int8,
        np.int16,
        np.int32,
        np.int64,
    ],
)
@pytest.mark.parametrize(
    "csr_major_index_dtype", [np.uint16, np.int32, np.uint32, np.int64]
)
@pytest.mark.parametrize(
    "csr_minor_index_dtype", [np.uint16, np.int32, np.uint32, np.int64]
)
def test_construction(
    shape: tuple[int, int],
    csr_major_index_dtype: np.dtype[Any],
    csr_minor_index_dtype: np.dtype[Any],
    value_dtype: np.dtype[Any],
    context: clib.SOMAContext,
) -> None:

    rng = np.random.default_rng()
    sp = sparse.random(
        shape[0], shape[1], density=0.01, dtype=value_dtype, random_state=rng
    )

    if sp.nnz >= np.iinfo(csr_major_index_dtype).max:
        # only occur if we mess up the test params
        pytest.skip(reason="NNZ is too large for index type.")

    indptr = np.empty((sp.shape[0] + 1), dtype=csr_major_index_dtype)
    indices = np.empty((sp.nnz,), dtype=csr_minor_index_dtype)
    data = np.empty((sp.nnz,), dtype=sp.dtype)

    fastercsx.compress_coo(
        context,
        shape,
        (sp.row.astype(np.int64),),
        (sp.col.astype(np.int64),),
        (sp.data,),
        indptr,
        indices,
        data,
    )
    fastercsx.sort_indices(context, indptr, indices, data)

    # Verify equality with SciPy constructed CSR
    csr = sp.tocsr()
    assert np.array_equal(indptr, csr.indptr)
    assert np.array_equal(indices, csr.indices)
    assert np.array_equal(data, csr.data)

    # Verify no dups (which scipy.sparse will sum at create time)
    ncsr = sparse.csr_matrix((data, indices, indptr), shape=shape)
    assert ncsr.has_canonical_format
    assert (ncsr != csr).nnz == 0


@pytest.mark.parametrize(
    "shape,density",
    [
        ((100_000, 1_000), 0.1),
        ((1_000_001, 73), 0.1),
        ((463459, 181), 0.1),  # nnz == 8M-1
        ((8388608, 10), 0.1),  # nnz == 8M
        ((27962027, 3), 0.1),  # nnz == 8M+1
    ],
)
@pytest.mark.parametrize("concurrency", [1, 2, None])
def test_partitioning(
    shape: tuple[int, int],
    density: float,
    context: clib.SOMAContext,
) -> None:
    """
    The implementation (currently) starts to partition work in two different ways:
    1. row counting: when the nnz is > 8*1024**2
    2. reorg of minor/data: if row count > 1
    """
    rng = np.random.default_rng()
    sp = sparse.random(
        shape[0],
        shape[1],
        dtype=np.uint8,
        format="coo",
        density=density,
        random_state=rng,
    )

    indptr = np.empty((sp.shape[0] + 1), dtype=np.int32)
    indices = np.empty((sp.nnz,), dtype=np.int32)
    data = np.empty((sp.nnz,), dtype=sp.dtype)

    fastercsx.compress_coo(
        context,
        sp.shape,
        (sp.row,),
        (sp.col,),
        (sp.data,),
        indptr,
        indices,
        data,
    )
    fastercsx.sort_indices(context, indptr, indices, data)

    # Verify equality with SciPy constructed CSR
    csr = sp.tocsr()
    assert np.array_equal(indptr, csr.indptr)
    assert np.array_equal(indices, csr.indices)
    assert np.array_equal(data, csr.data)


@pytest.mark.parametrize(
    "shape,density",
    [
        ((10000, 100), 0.1),
        ((100_001, 83), 0.1),
    ],
)
@pytest.mark.parametrize("nchunks", [1, 2, 3, 7, 10, 13])
def test_multichunk(
    nchunks: int, shape: tuple[int, int], density: float, context: clib.SOMAContext
) -> None:
    """check that multi-chunk COO input functions correctly"""

    rng = np.random.default_rng()
    sp = sparse.random(
        shape[0],
        shape[1],
        dtype=np.uint8,
        format="coo",
        density=density,
        random_state=rng,
    )

    row_chunks = tuple(np.array_split(sp.row, nchunks))
    col_chunks = tuple(np.array_split(sp.col, nchunks))
    data_chunks = tuple(np.array_split(sp.data, nchunks))

    indptr = np.empty((sp.shape[0] + 1), dtype=np.int32)
    indices = np.empty((sp.nnz,), dtype=np.int32)
    data = np.empty((sp.nnz,), dtype=sp.dtype)

    fastercsx.compress_coo(
        context,
        sp.shape,
        row_chunks,
        col_chunks,
        data_chunks,
        indptr,
        indices,
        data,
    )
    fastercsx.sort_indices(context, indptr, indices, data)

    # Verify equality with SciPy constructed CSR
    csr = sp.tocsr()
    assert np.array_equal(indptr, csr.indptr)
    assert np.array_equal(indices, csr.indices)
    assert np.array_equal(data, csr.data)


"""
Other stuff to test:

1. all error checking
2. boundary conditions on index types, e.g.,
    - nnz is too large to fit into major/minor axis
    - bogus row/val in COO is negative or exceeds n_row/n_col
    - Bp is too small
3. compress w/ bad shape
4. copy_to_dense w/ and w/o slice
5. count rows
6. sort gets handed bogus shaped arrays
7. sort gets handed data with corrupted indptr
8. compress gets handled dups and/or coords outside shape

"""
