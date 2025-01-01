"""Test the fastercsx pybind11 C++ layer."""

from __future__ import annotations

import sys
from typing import Any

import numpy as np
import pytest
from scipy import sparse

import tiledbsoma.pytiledbsoma as clib
import tiledbsoma.pytiledbsoma.fastercsx as fastercsx

NON_NATIVE_BYTEORDER = ">" if sys.byteorder == "little" else "<"


@pytest.fixture
def concurrency() -> int | None:
    return None


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng()


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
        # only occurs if we mess up the test params
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
    assert fastercsx.sort_csx_indices(context, indptr, indices, data)

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
    assert fastercsx.sort_csx_indices(context, indptr, indices, data)

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
    """Check that multi-chunk COO input functions correctly"""

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
    assert fastercsx.sort_csx_indices(context, indptr, indices, data)

    # Verify equality with SciPy constructed CSR
    csr = sp.tocsr()
    assert np.array_equal(indptr, csr.indptr)
    assert np.array_equal(indices, csr.indices)
    assert np.array_equal(data, csr.data)


def test_sort_csx_indices_bad_args(
    rng: np.random.Generator, context: clib.SOMAContext
) -> None:

    sp = sparse.random(2830, 212, density=0.1, dtype=np.int8, random_state=rng).tocsr()
    p, j, d = sp.indptr, sp.indices, sp.data

    with pytest.raises(TypeError):
        fastercsx.sort_csx_indices(None, p, j, d)
    with pytest.raises(TypeError):
        fastercsx.sort_csx_indices(0, p, j, d)
    with pytest.raises(TypeError):
        fastercsx.sort_csx_indices(context, None, j, d)
    with pytest.raises(TypeError):
        fastercsx.sort_csx_indices(context, p, None, d)
    with pytest.raises(TypeError):
        fastercsx.sort_csx_indices(context, p, j, None)

    with pytest.raises(ValueError):
        fastercsx.sort_csx_indices(context, p[:-1], j, d)

    with pytest.raises(ValueError):
        fastercsx.sort_csx_indices(context, p, j[1:], d)

    with pytest.raises(ValueError):
        fastercsx.sort_csx_indices(context, p, j, d[1:])

    with pytest.raises(ValueError):
        fastercsx.sort_csx_indices(context, np.zeros_like(p), j, d)

    with pytest.raises(OverflowError):
        pbad = p.copy()
        pbad[1] = sp.nnz + 1000000
        fastercsx.sort_csx_indices(context, pbad, j, d)

    with pytest.raises(OverflowError):
        pbad = p.copy()
        pbad[1] = -1
        fastercsx.sort_csx_indices(context, pbad, j, d)

    # non-native byteorder should throw
    with pytest.raises(ValueError):
        fastercsx.sort_csx_indices(context, p.astype(f"{NON_NATIVE_BYTEORDER}i4"), j, d)
    with pytest.raises(ValueError):
        fastercsx.sort_csx_indices(context, p, j.astype(f"{NON_NATIVE_BYTEORDER}i4"), d)
    with pytest.raises(ValueError):
        fastercsx.sort_csx_indices(context, p, j, d.astype(f"{NON_NATIVE_BYTEORDER}i4"))


def test_compress_coo_bad_args(
    rng: np.random.Generator, context: clib.SOMAContext
) -> None:

    sp = sparse.random(2830, 212, density=0.1, dtype=np.int8, random_state=rng)
    i, j, d = sp.row, sp.col, sp.data

    indptr = np.empty((sp.shape[0] + 1), dtype=np.int32)
    indices = np.empty((sp.nnz,), dtype=np.int32)
    data = np.empty((sp.nnz,), dtype=sp.dtype)

    with pytest.raises(TypeError):
        fastercsx.compress_coo(None, sp.shape, (i,), (j,), (d,), indptr, indices, data)
    with pytest.raises(TypeError):
        fastercsx.compress_coo(context, None, (i,), (j,), (d,), indptr, indices, data)

    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context,
            (sp.shape[0] + 1, sp.shape[1]),
            (i,),
            (j,),
            (d,),
            indptr,
            indices,
            data,
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, (0, sp.shape[1]), (i,), (j,), (d,), indptr, indices, data
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, (-1, sp.shape[1]), (i,), (j,), (d,), indptr, indices, data
        )
    with pytest.raises(IndexError):
        fastercsx.compress_coo(
            context, (sp.shape[0], 0), (i,), (j,), (d,), indptr, indices, data
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, (sp.shape[0], -1), (i,), (j,), (d,), indptr, indices, data
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, sp.shape, (i,), (j,), (d,), indptr[:-1].copy(), indices, data
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, sp.shape, (i,), (j,), (d,), indptr, indices[:-1].copy(), data
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, sp.shape, (i,), (j,), (d,), indptr, indices, data[1:].copy()
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, sp.shape, (i[1:],), (j,), (d,), indptr, indices, data
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, sp.shape, (i,), (j[1:],), (d,), indptr, indices, data
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context, sp.shape, (i,), (j,), (d[1:],), indptr, indices, data
        )

    # non-native byteorder should throw
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context,
            sp.shape,
            (i.astype(f"{NON_NATIVE_BYTEORDER}i4"),),
            (j,),
            (d,),
            indptr,
            indices,
            data,
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context,
            sp.shape,
            (i,),
            (j.astype(f"{NON_NATIVE_BYTEORDER}i4"),),
            (d,),
            indptr,
            indices,
            data,
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context,
            sp.shape,
            (i,),
            (j,),
            (d.astype(f"{NON_NATIVE_BYTEORDER}i4"),),
            indptr,
            indices,
            data,
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context,
            sp.shape,
            (i,),
            (j,),
            (d,),
            indptr.astype(f"{NON_NATIVE_BYTEORDER}i4"),
            indices,
            data,
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context,
            sp.shape,
            (i,),
            (j,),
            (d,),
            indptr,
            indices.astype(f"{NON_NATIVE_BYTEORDER}i4"),
            data,
        )
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context,
            sp.shape,
            (i,),
            (j,),
            (d,),
            indptr,
            indices,
            data.astype(f"{NON_NATIVE_BYTEORDER}i4"),
        )


def test_ragged_chunk_error(
    rng: np.random.Generator, context: clib.SOMAContext
) -> None:
    # module assumes all chunks are regular across all input arrays. Check that this
    # is enforced.
    nnz = 3
    shape = (100, 20)
    indptr = np.empty(shape[0] + 1, dtype=np.int32)
    indices = np.empty(nnz, dtype=np.int32)
    data = np.empty(nnz, dtype=np.int8)
    with pytest.raises(ValueError):
        fastercsx.compress_coo(
            context,
            shape,
            (np.array([0]), np.array([1, 2])),
            (np.array([0]), np.array([1, 2])),
            (np.array([0, 1], dtype=np.int8), np.array([2], dtype=np.int8)),
            indptr,
            indices,
            data,
        )


def test_empty_first_frag(rng: np.random.Generator, context: clib.SOMAContext) -> None:
    fastercsx.compress_coo(
        context,
        shape := (56, 18),
        (
            np.array([], dtype=np.int32),
            i := np.array([9, 35, 20, 7, 46, 41, 16, 28, 8, 19], dtype=np.int32),
        ),
        (
            np.array([], dtype=np.int32),
            j := np.array([15, 14, 11, 9, 4, 0, 0, 5, 3, 1], dtype=np.int32),
        ),
        (
            np.array([], dtype=np.int8),
            d := np.array([85, 63, 51, 26, 30, 4, 7, 1, 17, 81], dtype=np.int8),
        ),
        indptr := np.zeros(57, dtype=np.int32),
        indices := np.empty((10,), dtype=np.int32),
        data := np.empty((10,), dtype=np.int8),
    )

    csr = sparse.csr_matrix((d, (i, j)), shape=shape)
    assert (csr != sparse.csr_matrix((data, indices, indptr), shape=shape)).nnz == 0
    assert np.array_equal(csr.indptr, indptr)
    assert np.array_equal(csr.indices, indices)
    assert np.array_equal(csr.data, data)
