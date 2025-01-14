"""Hypothesis tests for fastercsx module."""

from typing import Any, Literal, TypeAlias

import numpy as np
import numpy.typing as npt
import pyarrow as pa
import pytest
import scipy.sparse as sparse
from hypothesis import given, settings
from hypothesis import strategies as st

import tiledbsoma as soma
import tiledbsoma._fastercsx as fastercsx
import tiledbsoma.pytiledbsoma.fastercsx as clib_fastercsx

from tests.ht._ht_util import (
    arrow_array,
    arrow_array_fast,
    contiguous_slices,
    random_length_tuple,
    resolve_dtype,
    split_arrow_array,
    splitss,
)

# Supported types, i.e., these should work fine
CooIndexTypes = sorted(np.dtype(t) for t in (np.int32, np.int64))
CsxIndexTypes = sorted(np.dtype(t) for t in (np.int32, np.int64, np.uint16, np.uint32))
ValueTypes = sorted(
    np.dtype(t)
    for t in (
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.float32,
        np.float64,
    )
)

NDArrayIndex: TypeAlias = npt.NDArray[np.integer[Any]]
NDArrayNumber: TypeAlias = npt.NDArray[np.integer[Any] | np.floating[Any]]


def limit_value_range_element_strategy(
    dtype: np.dtype, divisor: float
) -> dict[str, Any] | None:
    if dtype.kind == "f":
        info = np.finfo(dtype)
        return {"min_value": -info.max / divisor, "max_value": info.max / divisor}
    if dtype.kind in ["i", "u"]:
        info = np.iinfo(dtype)
        return {"min_value": info.min // divisor, "max_value": info.max // divisor}
    return None


@st.composite
def coo_ijd(
    draw: st.DrawFn,
    dtype: npt.DTypeLike | pa.DataType | st.SearchStrategy[npt.DTypeLike | pa.DataType],
    shape: tuple[int, int] | st.SearchStrategy[tuple[int, int]],
    *,
    density: float | st.SearchStrategy[float] = 0.01,
    unique: bool = False,
) -> tuple[
    tuple[npt.NDArray[Any], ...],
    tuple[npt.NDArray[Any], ...],
    tuple[npt.NDArray[Any], ...],
]:

    dtype = resolve_dtype(draw, dtype)
    shape = draw(shape) if isinstance(shape, st.SearchStrategy) else shape
    assert isinstance(shape, tuple) and len(shape) == 2

    density = draw(density) if isinstance(density, st.SearchStrategy) else density
    assert isinstance(density, float) and 0 < density <= 1

    nnz = int(shape[0] * shape[1] * density)
    coord_dtype = draw(st.sampled_from(CooIndexTypes))

    """
    if not unique, we need to be cognizant of the potential to overflow
    when duplicates are summed (the default behavior for `to_scipy`). This
    can easily cause some types to overflow, and others to lose precision,
    which makes equality comparisons tricky.

    To avoid this, ONLY when `not unique`, constrain the range of generated
    values to a very limited range (currently 1/128th of the full range).
    This is extremely unlikely to overflow as it would require 128 identical
    coordinates to be drawn.

    In the case of `unique`, draw from the full range for the type.

    Currently, the only edge case that fails to do the right thing is timestamp
    generation (datetime64), as the underlying search strategy used does not
    obey min_value/max_value for that type.  TODO - FIXME.
    """
    if not unique:
        i = draw(
            arrow_array_fast(
                dtype=coord_dtype,
                shape=nnz,
                min_value=0,
                max_value=shape[0] - 1,
            )
        )
        j = draw(
            arrow_array_fast(
                dtype=coord_dtype,
                shape=nnz,
                min_value=0,
                max_value=shape[1] - 1,
            )
        )
        d = draw(
            arrow_array(
                dtype=dtype,
                shape=nnz,
                elements=limit_value_range_element_strategy(dtype, 128),
            )
        )

    else:
        # draw unique points, then split into I/J
        rng = np.random.default_rng(seed=draw(st.integers(min_value=0)))
        points = rng.choice(shape[0] * shape[1], size=nnz, replace=False)
        i, j = np.divmod(points, shape[1])
        i = pa.array(i, type=pa.from_numpy_dtype(coord_dtype))
        j = pa.array(j, type=pa.from_numpy_dtype(coord_dtype))
        d = draw(arrow_array(dtype=dtype, shape=nnz))

    if draw(st.booleans()):
        return (i.to_numpy(),), (j.to_numpy(),), (d.to_numpy(),)

    # else split into a chunked array
    n_splits = draw(st.integers(min_value=0, max_value=max(0, len(d) // 10)))
    split_points = draw(splitss(n_splits=n_splits, max_value=len(d)))
    return (
        tuple(c.to_numpy() for c in split_arrow_array(i, split_points).chunks),
        tuple(c.to_numpy() for c in split_arrow_array(j, split_points).chunks),
        tuple(c.to_numpy() for c in split_arrow_array(d, split_points).chunks),
    )


@given(
    do=st.data(),
    value_dtype=st.sampled_from(ValueTypes),
    unique=st.booleans(),
    shape=st.tuples(
        st.integers(min_value=0, max_value=1024),
        st.integers(min_value=0, max_value=1024),
    ),
    context=st.from_type(soma.SOMATileDBContext),
)
@settings(max_examples=500)
def test_fastercsx_clib_compress_coo(
    do: st.DataObject,
    value_dtype: np.dtype,
    unique: bool,
    shape: tuple[int, int],
    context: soma.SOMATileDBContext,
) -> None:
    i, j, d = do.draw(coo_ijd(dtype=value_dtype, shape=shape, unique=unique))
    nnz = sum(len(c) for c in i)
    assert nnz <= np.prod(shape)
    index_dtype = do.draw(
        st.sampled_from([t for t in CsxIndexTypes if np.iinfo(t).max >= nnz])
    )

    indptr = np.empty(shape[0] + 1, dtype=index_dtype)
    indices = np.empty(nnz, dtype=index_dtype)
    data = np.empty(nnz, dtype=value_dtype)
    clib_fastercsx.compress_coo(
        context.native_context, shape, i, j, d, indptr, indices, data
    )

    # check with the oracle. Be careful if dups allowed, as summing dups
    # in floats will be _approximately_ equal, not exactly equal
    csr = sparse.csr_matrix(
        (data, indices, indptr), shape=shape, dtype=value_dtype, copy=False
    )
    if not unique:
        csr.sum_duplicates()
    csr.sort_indices()

    scipy_csr = sparse.csr_matrix(
        (np.concatenate(d), (np.concatenate(i), np.concatenate(j))),
        shape=shape,
        dtype=value_dtype,
    )

    assert np.array_equal(csr.indptr, scipy_csr.indptr)
    assert np.array_equal(csr.indices, scipy_csr.indices)

    # XXX the non-unique case has several issues:
    #    - dups are added, which can overflow
    #    - dups are added, which for floats may not be exactly eq in some situations
    # XXX cleanup debug code
    if not (
        np.allclose(
            csr.data,
            scipy_csr.data,
            equal_nan=True if value_dtype.kind == "f" else False,
        )
        if not unique
        else np.array_equal(
            csr.data,
            scipy_csr.data,
            equal_nan=True if value_dtype.kind == "f" else False,
        )
    ):
        print(csr.data.dtype, scipy_csr.data.dtype)
        print(csr.data)
        print(scipy_csr.data)
        print(csr.data - scipy_csr.data)
    assert (
        np.allclose(
            csr.data,
            scipy_csr.data,
            equal_nan=True if value_dtype.kind == "f" else False,
        )
        if not unique
        else np.array_equal(
            csr.data,
            scipy_csr.data,
            equal_nan=True if value_dtype.kind == "f" else False,
        )
    )


@given(
    shape=random_length_tuple(elements=st.integers(), max_length=3),
    i=random_length_tuple(st.from_type(npt.NDArray[Any]), max_length=4),
    j=random_length_tuple(elements=st.from_type(npt.NDArray[Any]), max_length=4),
    d=random_length_tuple(elements=st.from_type(npt.NDArray[Any]), max_length=4),
    indptr=st.from_type(npt.NDArray[Any]),
    indices=st.from_type(npt.NDArray[Any]),
    data=st.from_type(npt.NDArray[Any]),
    context=st.from_type(soma.SOMATileDBContext),
)
def test_fuzz_fastercsx_clib_compress_coo(
    shape,
    i: tuple[npt.NDArray[Any], ...],
    j: tuple[npt.NDArray[Any], ...],
    d: tuple[npt.NDArray[Any], ...],
    indptr: npt.NDArray[Any],
    indices: npt.NDArray[Any],
    data: npt.NDArray[Any],
    context: soma.SOMATileDBContext,
) -> None:
    # TODO: exclude the rare case that would pass
    with pytest.raises(Exception):
        clib_fastercsx.compress_coo(
            context.native_context, shape, i, j, d, indptr, indices, data
        )


@given(
    indptr=st.from_type(npt.NDArray[Any]).filter(
        lambda a: a.dtype not in CsxIndexTypes
    ),
    indices=st.from_type(npt.NDArray[Any]).filter(
        lambda a: a.dtype not in CsxIndexTypes
    ),
    data=st.from_type(npt.NDArray[Any]).filter(lambda a: a.dtype not in ValueTypes),
    context=st.from_type(soma.SOMATileDBContext),
)
@settings(max_examples=250)
def test_fuzz_fastercsx_clib_sort_csx_indices(
    indptr: npt.NDArray[Any],
    indices: npt.NDArray[Any],
    data: npt.NDArray[Any],
    context: soma.SOMATileDBContext,
) -> None:
    # TODO: exclude the rare case that would pass
    with pytest.raises(Exception):
        clib_fastercsx.sort_csx_indices(context.native_context, indptr, indices, data)


@given(
    major_idx_start=st.integers(),
    major_idx_end=st.integers(),
    shape=random_length_tuple(elements=st.integers(), max_length=3),
    format=st.text(),
    indptr=st.from_type(npt.NDArray[Any]),
    indices=st.from_type(npt.NDArray[Any]),
    data=st.from_type(npt.NDArray[Any]),
    out=st.from_type(npt.NDArray[Any]),
    context=st.from_type(soma.SOMATileDBContext),
)
@settings(max_examples=250)
def test_fuzz_fastercsx_clib_copy_csx_to_dense(
    major_idx_start: int,
    major_idx_end: int,
    shape: tuple[int, int],
    format: str,
    indptr: npt.NDArray[Any],
    indices: npt.NDArray[Any],
    data: npt.NDArray[Any],
    out: npt.NDArray[Any],
    context: soma.SOMATileDBContext,
) -> None:
    # TODO: exclude the rare case that would pass
    with pytest.raises(Exception):
        clib_fastercsx.copy_csx_to_dense(
            context.native_context,
            major_idx_start,
            major_idx_end,
            shape,
            format,
            indptr,
            indices,
            data,
            out,
        )


@given(
    do=st.data(),
    value_dtype=st.sampled_from(ValueTypes),
    unique=st.booleans(),
    shape=st.tuples(
        st.integers(min_value=0, max_value=1024),
        st.integers(min_value=0, max_value=1024),
    ),
    make_sorted=st.booleans(),
    format=st.sampled_from(["csc", "csr"]),
    context=st.from_type(soma.SOMATileDBContext),
)
@settings(max_examples=500)
def test_fastercsx_from_ijd(
    do: st.DataObject,
    value_dtype: np.dtype,
    unique: bool,
    shape: tuple[int, int],
    format: Literal["csc", "csr"],
    make_sorted: bool,
    context: soma.SOMATileDBContext,
) -> None:
    i, j, d = do.draw(coo_ijd(dtype=value_dtype, shape=shape, unique=unique))
    assert all(a.dtype == value_dtype for a in d)

    cm = fastercsx.CompressedMatrix.from_ijd(
        i, j, d, shape, format, make_sorted, context
    ).to_scipy()
    assert cm.dtype == value_dtype

    # compare to oracle
    scipy_cm = sparse.coo_matrix(
        (np.concatenate(d), (np.concatenate(i), np.concatenate(j))),
        shape=shape,
        dtype=value_dtype,
    ).asformat(format)
    assert scipy_cm.has_canonical_format

    if not make_sorted or not unique:
        cm.sum_duplicates()

    assert np.array_equal(cm.indptr, scipy_cm.indptr)
    assert np.array_equal(cm.indices, scipy_cm.indices)

    # XXX cleanup debug code
    assert cm.data.dtype == scipy_cm.data.dtype
    if not (
        np.allclose(
            cm.data, scipy_cm.data, equal_nan=True if value_dtype.kind == "f" else False
        )
        if not unique
        else np.array_equal(
            cm.data, scipy_cm.data, equal_nan=True if value_dtype.kind == "f" else False
        )
    ):
        print(cm.data.dtype, scipy_cm.data.dtype)
        print(cm.data)
        print(scipy_cm.data)
    assert (
        np.allclose(
            cm.data, scipy_cm.data, equal_nan=True if value_dtype.kind == "f" else False
        )
        if not unique
        else np.array_equal(
            cm.data, scipy_cm.data, equal_nan=True if value_dtype.kind == "f" else False
        )
    )


@given(
    do=st.data(),
    value_dtype=st.sampled_from(ValueTypes),
    unique=st.booleans(),
    shape=st.tuples(
        st.integers(min_value=0, max_value=1024),
        st.integers(min_value=0, max_value=1024),
    ),
    make_sorted=st.booleans(),
    format=st.sampled_from(["csc", "csr"]),
    context=st.from_type(soma.SOMATileDBContext),
)
@settings(max_examples=500)
def test_fastercsx_to_scipy(
    do: st.DataObject,
    value_dtype: np.dtype,
    unique: bool,
    shape: tuple[int, int],
    format: Literal["csc", "csr"],
    make_sorted: bool,
    context: soma.SOMATileDBContext,
) -> None:
    i, j, d = do.draw(coo_ijd(dtype=value_dtype, shape=shape, unique=unique))
    cm = fastercsx.CompressedMatrix.from_ijd(
        i, j, d, shape, format, make_sorted, context
    )

    # compare to oracle
    scipy_cm = sparse.coo_matrix(
        (np.concatenate(d), (np.concatenate(i), np.concatenate(j))),
        shape=shape,
        dtype=value_dtype,
    ).asformat(format)
    assert scipy_cm.has_canonical_format

    major_index_slice = do.draw(contiguous_slices(shape[0]))

    cm_slc = cm.to_scipy(major_index_slice)
    if not make_sorted or not unique:
        cm_slc.sum_duplicates()
    assert cm_slc.has_canonical_format

    scipy_slc = (
        scipy_cm[major_index_slice]
        if format == "csr"
        else scipy_cm[:, major_index_slice]
    )

    assert np.array_equal(cm_slc.indptr, scipy_slc.indptr)
    assert np.array_equal(cm_slc.indices, scipy_slc.indices)
    assert (
        np.allclose(
            cm_slc.data,
            scipy_slc.data,
            equal_nan=True if value_dtype.kind == "f" else False,
        )
        if not unique
        else np.array_equal(
            cm_slc.data,
            scipy_slc.data,
            equal_nan=True if value_dtype.kind == "f" else False,
        )
    )
