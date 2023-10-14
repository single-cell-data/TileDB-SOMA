# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Read iterators.
"""
from __future__ import annotations

import abc
from typing import Iterator, Literal, Optional, Sequence, Tuple, TypeVar, Union

import numba
import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
import scipy.sparse as sparse
import somacore
from somacore import options
from somacore.query._eager_iter import EagerIterator

# This package's pybind11 code
import tiledbsoma.pytiledbsoma as clib

from . import _util
from ._exception import SOMAError
from ._types import NTuple


class TableReadIter(somacore.ReadIter[pa.Table]):
    """Iterator over `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_ elements"""

    def __init__(self, sr: clib.SOMAArray):
        self.sr = sr

    def __next__(self) -> pa.Table:
        arrow_table = self.sr.read_next()
        if arrow_table is None:
            raise StopIteration

        return arrow_table

    def concat(self) -> pa.Table:
        """Concatenate remainder of iterator, and return as a single `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_"""
        return pa.concat_tables(self)


RT = TypeVar("RT")


class SparseTensorReadIterBase(somacore.ReadIter[RT], metaclass=abc.ABCMeta):
    """Private implementation class"""

    def __init__(self, sr: clib.SOMAArray, shape: NTuple):
        self.sr = sr
        self.shape = shape

    @abc.abstractmethod
    def _from_table(self, arrow_table: pa.Table) -> RT:
        raise NotImplementedError()

    def __next__(self) -> RT:
        arrow_table = self.sr.read_next()
        if arrow_table is None:
            raise StopIteration

        return self._from_table(arrow_table)

    def concat(self) -> RT:
        """Returns all the requested data in a single operation.

        If some data has already been retrieved using ``next``, this will return
        the rest of the data after that is already returned.
        """
        arrow_tables = pa.concat_tables(TableReadIter(self.sr))
        return self._from_table(arrow_tables)


class SparseCOOTensorReadIter(SparseTensorReadIterBase[pa.SparseCOOTensor]):
    """Iterator over `Arrow SparseCOOTensor <https://arrow.apache.org/docs/cpp/api/tensor.html>`_ elements"""

    def _from_table(self, arrow_table: pa.Table) -> pa.SparseCOOTensor:
        coo_data = arrow_table.column("soma_data").to_numpy()
        coo_coords = np.array(
            [
                arrow_table.column(f"soma_dim_{n}").to_numpy()
                for n in range(len(self.shape))
            ]
        ).T
        return pa.SparseCOOTensor.from_numpy(coo_data, coo_coords, shape=self.shape)


# TODO: when we can access the schema via SOMAArray, remove this hard-wire
# (even though it is safe given the SOMA spec hard-wires the names)
_SparseNDArrayAxisNames = ("soma_dim_0", "soma_dim_1")


def scipy_sparse_iter(
    sr: clib.SOMAArray,
    coords: options.SparseDFCoords,
    axis: Literal[0, 1],
    step: int,
    compress: bool,
    reindex_sparse_axis: bool,
) -> Iterator[Union[sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix]]:
    # setup

    shape = sr.shape
    if len(shape) != 2 or len(coords) > 2:
        raise SOMAError("SciPy spmatrix iterator compatible only with 2D arrays")

    # TODO: need to add error checks for result_order and column_names compat.

    # normalize params tuple
    coords = tuple(coords[i] if i < len(coords) else None for i in range(2))

    # Step via the major axis
    major_axis = axis
    minor_axis = 1 - axis

    # Materialize all minor coordinates as we must yield them back
    minor_coords: npt.NDArray[np.int64] = np.concatenate(
        list(_coords_strider(coords[minor_axis], shape[minor_axis], shape[minor_axis]))
    )

    # Create reindexers. Only needed under specific conditions
    minor_axis_indexer: Optional[pd.Index] = None
    if reindex_sparse_axis:
        minor_axis_indexer = pd.Index(minor_coords)

    # TODO tighten up typing.  Could split this into two generator functions
    def _generate() -> Iterator[sparse.spmatrix]:
        for major_coords, coo_tbl in EagerIterator(
            _partitioned_arrow_table_reader(sr, coords, minor_coords, step, major_axis)
        ):
            # Arrow table multi-column sort
            st = coo_tbl.sort_by(
                [
                    (_SparseNDArrayAxisNames[major_axis], "ascending"),
                    (_SparseNDArrayAxisNames[minor_axis], "ascending"),
                ]
            )

            _sp_shape = list(shape)
            _sp_shape[major_axis] = len(major_coords)
            if reindex_sparse_axis:
                _sp_shape[minor_axis] = len(minor_coords)
            sp_shape = tuple(_sp_shape)

            if not compress:
                i = st.column("soma_dim_0").to_numpy()  # copies if multi-chunk
                j = st.column("soma_dim_1").to_numpy()  # copies if multi-chunk
                d = st.column("soma_data").to_numpy()  # copies if multi-chunk
                indices = [i, j]

                # reindex major axis
                indices[major_axis] = pd.Index(major_coords).get_indexer(indices[major_axis])  # type: ignore[no-untyped-call]

                # optionally reindex minor axis
                if minor_axis_indexer is not None:
                    indices[minor_axis] = minor_axis_indexer.get_indexer(  # type: ignore[no-untyped-call]
                        indices[minor_axis]
                    )

                i, j = indices
                sp = sparse.coo_matrix((d, (i, j)), shape=sp_shape)

                # SOMA disallows duplicates
                sp.has_canonical_format = True

            else:
                j = st.column(minor_axis).to_numpy()  # copies if multi-chunk
                d = st.column("soma_data").to_numpy()  # copies if multi-chunk
                indptr = _create_indptr(
                    major_coords,
                    tuple(c.to_numpy() for c in st.column(major_axis).iterchunks()),
                )

                # optional: reindex the minor dimension
                if minor_axis_indexer is not None:
                    j = minor_axis_indexer.get_indexer(j)  # type: ignore[no-untyped-call]

                # hack to bypass https://github.com/scipy/scipy/issues/11496
                cls = sparse.csr_matrix if axis == 0 else sparse.csc_matrix
                sp = cls((0, 0), dtype=d.dtype)
                sp.data = d
                sp.indices = j
                sp.indptr = indptr
                sp._shape = sp_shape

                # XXX debugging - move to unit tests
                sp.check_format()
                assert sp.has_canonical_format

            # debugging - move to unit test
            assert sp.has_canonical_format

            yield (major_coords, minor_coords), sp

    return _generate()


@numba.jit(nopython=True, nogil=True)  # type: ignore[misc]  # See https://github.com/numba/numba/issues/7424
def _create_indptr(
    index: npt.NDArray[np.int64], coords: Tuple[npt.NDArray[np.int64]]
) -> npt.NDArray[np.int64]:
    """
    Generate indptr. Assumes/knows that both index and coords are monotonically increasing.
    Note: numba does not recognize Arrow types, so use numpy (w/ attention to minimal copy).
    """
    indptr = np.empty((len(index) + 1,), dtype=np.int64)

    # count nnz for each row
    c: int = 0  # the coords item
    n: int = 0  # the element in the coord item
    for i in range(len(index)):
        cnt: int = 0  # the nnz per index
        while c < len(coords) and coords[c][n] == index[i]:
            # coord = coords[c]
            while n < len(coords[c]) and coords[c][n] == index[i]:
                cnt += 1
                n += 1
            if n == len(coords[c]):
                c += 1
                n = 0

        indptr[i] = cnt

    # accumulate nnz to create indptr
    cumsum = 0
    for i in range(len(index)):
        tmp = indptr[i]
        indptr[i] = cumsum
        cumsum += tmp
    indptr[len(index)] = cumsum

    return indptr


def _coords_strider(
    coords: options.SparseDFCoord, length: int, stride: int
) -> Iterator[npt.NDArray[np.int64]]:
    """
    Private.

    Iterate over major coordinates, in stride sized steps, materializing each step as an
    ndarray of coordinate values. Will be sorted in ascending order.
    """
    # convert to either a slice or ndarray
    if coords is None:
        coords = slice(None)
    elif isinstance(coords, Sequence):
        coords = np.array(coords).astype(np.int64)
    elif isinstance(coords, pa.Array) or isinstance(coords, pa.ChunkedArray):
        coords = coords.to_numpy()

    if isinstance(coords, slice):
        _util.validate_slice(coords)  # NB: this enforces step == 1, assumed below
        start, stop, _step = coords.indices(length)
        assert _step == 1
        # NB: SOMA slices are inclusive - add one to stop
        yield from (
            np.arange(i, min(i + stride, stop + 1), dtype=np.int64)
            for i in range(start, stop + 1, stride)
        )

    else:
        assert isinstance(coords, np.ndarray) and coords.dtype == np.int64
        sorted_coords = np.sort(coords)
        for i in range(0, len(sorted_coords), stride):
            yield sorted_coords[i : i + stride]


# reader - yields tables with complete reads on major axis
def _partitioned_arrow_table_reader(
    sr: clib.SOMAArray,
    coords: options.SparseDFCoords,
    minor_coords: npt.NDArray[np.int64],
    step: int,
    major_axis: int,
) -> Iterator[Tuple[npt.NDArray[np.int64], pa.Table]]:
    """
    Private. Yields step-sized tables with complete read over major axis.
    """
    minor_axis: int = 1 - major_axis

    def arrow_table_reader() -> Iterator[pa.Table]:
        tbl = sr.read_next()
        while tbl is not None:
            yield tbl
            tbl = sr.read_next()

    for coord_chunk in _coords_strider(coords[major_axis], sr.shape[major_axis], step):
        sr.reset()
        sr.set_dim_points_int64(_SparseNDArrayAxisNames[major_axis], coord_chunk)
        sr.set_dim_points_int64(_SparseNDArrayAxisNames[minor_axis], minor_coords)
        tbl = pa.concat_tables(arrow_table_reader())
        yield coord_chunk, tbl
