"""Tests of Python wrappers."""

from __future__ import annotations

from typing import Any, Literal

import numpy as np
import pyarrow as pa
import pytest
import scipy.sparse as sparse
from typeguard import suppress_type_checks

import tiledbsoma as soma
import tiledbsoma._fastercsx as fastercsx

NP_VALUE_TYPES = [  # supported types - see fastercsx.h
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
]


@pytest.fixture
def concurrency() -> int | None:
    return None


@pytest.fixture
def rng() -> np.random.Generator:
    return np.random.default_rng()


@pytest.fixture
def context(concurrency: int | None) -> soma.SOMATileDBContext:
    if concurrency is None:
        return soma.SOMATileDBContext()
    else:
        return soma.SOMATileDBContext(
            tiledb_config={"soma.compute_concurrency_level": f"{concurrency}"}
        )


def assert_eq(sp: sparse.spmatrix, cm: fastercsx.CompressedMatrix) -> bool:

    sp_csx = sp.tocsr() if cm.format == "csr" else sp.tocsc()

    assert cm.data.shape == sp_csx.data.shape
    assert cm.indices.shape == sp_csx.indices.shape
    assert cm.indptr.shape == sp_csx.indptr.shape

    assert np.array_equal(cm.indptr, sp_csx.indptr)
    if cm.is_sorted:
        assert np.array_equal(cm.indices, sp_csx.indices)
        assert np.array_equal(cm.data, sp_csx.data)

    assert (sp_csx != cm.to_scipy()).nnz == 0

    assert sp_csx.dtype == cm.dtype
    assert sp_csx.nnz == cm.nnz
    assert sp_csx.shape == cm.shape

    return True


@pytest.mark.parametrize("format", ["csr", "csc"])
@pytest.mark.parametrize(
    "shape", [(0, 0), (1, 0), (0, 1), (10, 100), (100, 10), (9972, 1), (1, 9972)]
)
@pytest.mark.parametrize("value_dtype", NP_VALUE_TYPES)
def test_from_ijd(
    shape: tuple[int, int],
    value_dtype: np.dtype[Any],
    format: Literal["csc", "csr"],
    context: soma.SOMATileDBContext,
    rng: np.random.Generator,
) -> None:
    sp = sparse.random(
        shape[0], shape[1], density=0.01, dtype=value_dtype, random_state=rng
    )

    cm = fastercsx.CompressedMatrix.from_ijd(
        sp.row, sp.col, sp.data, sp.shape, format, make_sorted=True, context=context
    )
    assert_eq(sp, cm)
    assert cm.nbytes == (cm.indptr.nbytes + cm.indices.nbytes + cm.data.nbytes)
    assert cm.to_scipy().has_canonical_format


@pytest.mark.parametrize("format", ["csr", "csc"])
@pytest.mark.parametrize("n_tables", [1, 2, 3])
@pytest.mark.parametrize(
    "shape",
    [(0, 0), (1, 0), (0, 1), (10, 100), (100, 10), (9972, 1), (1, 9972), (10_001, 997)],
)
@pytest.mark.parametrize("value_dtype", NP_VALUE_TYPES)
def test_from_soma_array(
    n_tables: int,
    shape: tuple[int, int],
    value_dtype: np.dtype[Any],
    format: Literal["csc", "csr"],
    context: soma.SOMATileDBContext,
    rng: np.random.Generator,
) -> None:
    sp = sparse.random(
        shape[0],
        shape[1],
        density=0.1,
        dtype=value_dtype,
        random_state=rng,
        format="coo",
    )

    tables = []
    i = np.array_split(sp.row.astype(np.int64), n_tables)
    j = np.array_split(sp.col.astype(np.int64), n_tables)
    d = np.array_split(sp.data, n_tables)
    tables = [
        pa.Table.from_pydict(
            {
                "soma_dim_0": pa.array(i_),
                "soma_dim_1": pa.array(j_),
                "soma_data": pa.array(d_),
            }
        )
        for i_, j_, d_ in zip(i, j, d)
    ]
    cm = fastercsx.CompressedMatrix.from_soma(
        tables, sp.shape, format, make_sorted=True, context=context
    )
    assert_eq(sp, cm)
    assert cm.to_scipy().has_canonical_format


@pytest.mark.parametrize("format", ["csr", "csc"])
@pytest.mark.parametrize("n_tables", [1, 2, 3])
@pytest.mark.parametrize("n_chunks", [1, 2, 3, 7, 11])
@pytest.mark.parametrize(
    "shape",
    [(0, 0), (1, 0), (0, 1), (10, 100), (100, 10), (9972, 1), (1, 9972), (10_001, 997)],
)
@pytest.mark.parametrize("value_dtype", NP_VALUE_TYPES)
def test_from_soma_chunked_array(
    n_tables: int,
    n_chunks: int,
    shape: tuple[int, int],
    value_dtype: np.dtype[Any],
    format: Literal["csc", "csr"],
    context: soma.SOMATileDBContext,
    rng: np.random.Generator,
) -> None:
    sp = sparse.random(
        shape[0],
        shape[1],
        density=0.1,
        dtype=value_dtype,
        random_state=rng,
        format="coo",
    )

    i = np.array_split(sp.row.astype(np.int64), n_tables)
    j = np.array_split(sp.col.astype(np.int64), n_tables)
    d = np.array_split(sp.data, n_tables)
    tables = [
        pa.Table.from_pydict(
            {
                "soma_dim_0": pa.chunked_array(np.array_split(i_, n_chunks)),
                "soma_dim_1": pa.chunked_array(np.array_split(j_, n_chunks)),
                "soma_data": pa.chunked_array(np.array_split(d_, n_chunks)),
            }
        )
        for i_, j_, d_ in zip(i, j, d)
    ]

    assert len(tables) == n_tables
    assert tables[0]["soma_data"].num_chunks == n_chunks
    cm = fastercsx.CompressedMatrix.from_soma(
        tables, sp.shape, format, make_sorted=True, context=context
    )
    assert_eq(sp, cm)
    assert cm.to_scipy().has_canonical_format


@pytest.mark.parametrize("format", ["csr", "csc"])
@pytest.mark.parametrize("sorted", [True, False])
@pytest.mark.parametrize(
    "index", [None, slice(None), slice(0), slice(10), slice(0, 1), slice(1, -1)]
)
def test_to_scipy(
    format: Literal["csc", "csr"],
    sorted: bool,
    index: slice | None,
    context: soma.SOMATileDBContext,
    rng: np.random.Generator,
) -> None:
    sp = sparse.random(
        1000, 100, density=0.1, dtype=np.float32, random_state=rng, format="coo"
    )

    cm = fastercsx.CompressedMatrix.from_ijd(
        sp.row,
        sp.col,
        sp.data,
        sp.shape,
        format=format,
        make_sorted=sorted,
        context=context,
    )
    assert_eq(sp, cm)
    if format == "csr":
        assert (
            cm.to_scipy(index).tocsr() != sp.tocsr()[index or slice(None), :]
        ).nnz == 0
    else:
        assert (
            cm.to_scipy(index).tocsc() != sp.tocsc()[:, index or slice(None)]
        ).nnz == 0


@pytest.mark.parametrize("format", ["csr", "csc"])
@pytest.mark.parametrize("sorted", [True, False])
@pytest.mark.parametrize(
    "shape", [(0, 0), (1, 0), (0, 1), (10, 1000), (1000, 10), (3, 99), (99, 3)]
)
@pytest.mark.parametrize(
    "index", [None, slice(None), slice(0), slice(10), slice(0, 1), slice(1, -1)]
)
def test_to_numpy(
    format: Literal["csc", "csr"],
    sorted: bool,
    shape: tuple[int, int],
    index: slice | None,
    context: soma.SOMATileDBContext,
    rng: np.random.Generator,
) -> None:
    sp = sparse.random(
        shape[0],
        shape[1],
        density=0.1,
        dtype=np.float32,
        random_state=rng,
        format="coo",
    )

    cm = fastercsx.CompressedMatrix.from_ijd(
        sp.row,
        sp.col,
        sp.data,
        sp.shape,
        format=format,
        make_sorted=sorted,
        context=context,
    )
    assert_eq(sp, cm)
    assert np.array_equal(sp.toarray(), cm.to_numpy())
    if format == "csr":
        assert np.array_equal(
            sp.tocsr()[index or slice(None), :].toarray(), cm.to_numpy(index)
        )
    else:
        assert np.array_equal(
            sp.tocsc()[:, index or slice(None)].toarray(), cm.to_numpy(index)
        )


def test_bad_arguments(
    context: soma.SOMATileDBContext, rng: np.random.Generator
) -> None:
    """Test various bad argument types/values are caught."""
    sp = sparse.random(970, 31, density=0.01, dtype=np.float32, random_state=rng)

    with suppress_type_checks():  # w/o this, typeguard raises, which defeats the point of the test

        # context - bad type
        with pytest.raises(AttributeError):
            fastercsx.CompressedMatrix.from_ijd(
                sp.row, sp.col, sp.data, sp.shape, "csr", True, None
            )
        # context - missing
        with pytest.raises(TypeError):
            fastercsx.CompressedMatrix.from_ijd(
                sp.row, sp.col, sp.data, sp.shape, "csr", True
            )

        # unsupported data type
        unsup_data = np.full(sp.data.shape, "A", dtype="str")
        with pytest.raises(ValueError):
            fastercsx.CompressedMatrix.from_ijd(
                sp.row, sp.col, unsup_data, sp.shape, "csr", True, context
            )

        # unsupported index type
        with pytest.raises(ValueError):
            fastercsx.CompressedMatrix.from_ijd(
                sp.row.astype(np.float32),
                sp.col.astype(np.float32),
                sp.data,
                sp.shape,
                "csr",
                True,
                context,
            )

        # unsupported index type
        with pytest.raises(ValueError):
            fastercsx.CompressedMatrix.from_ijd(
                sp.row.astype(np.int16),
                sp.col.astype(np.int16),
                sp.data,
                sp.shape,
                "csr",
                True,
                context,
            )

        # mismatched index types
        with pytest.raises(TypeError):
            fastercsx.CompressedMatrix.from_ijd(
                sp.row.astype(np.int32),
                sp.col.astype(np.int64),
                sp.data,
                sp.shape,
                "csr",
                True,
                context,
            )


def test_bad_shapes(context: soma.SOMATileDBContext, rng: np.random.Generator) -> None:
    """Test bad/mismatched shape."""
    sp = sparse.identity(32, dtype=np.int8).tocoo()
    for shp in [
        (sp.shape[0] - 1, sp.shape[1]),
        (sp.shape[0], sp.shape[1] - 1),
    ]:
        with pytest.raises(IndexError):
            # this will also log an error - which is expected behavior for parallel_for/thread_pool
            fastercsx.CompressedMatrix.from_ijd(
                sp.row, sp.col, sp.data, shp, "csr", True, context
            )


@pytest.mark.parametrize("format", ["csr", "csc"])
@pytest.mark.parametrize("make_sorted", [True, False])
def test_duplicates(
    format: Literal["csc", "csr"],
    make_sorted: bool,
    context: soma.SOMATileDBContext,
    rng: np.random.Generator,
) -> None:
    shape = (10, 10)
    i = np.array([0, 0, 1, 1], dtype=np.int64)
    j = np.array([0, 0, 1, 1], dtype=np.int64)
    d = np.arange(len(i), dtype=np.int8)

    cm = fastercsx.CompressedMatrix.from_ijd(
        i, j, d, shape, format, make_sorted, context
    )
    assert (
        sparse.coo_matrix((d, (i, j)), shape=shape).asformat(format) != cm.to_scipy()
    ).nnz == 0
