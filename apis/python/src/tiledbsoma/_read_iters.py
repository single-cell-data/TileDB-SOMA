# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Read iterators.
"""
from __future__ import annotations

import abc
import concurrent.futures
from typing import (
    TYPE_CHECKING,
    Dict,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
    Union,
    cast,
)

import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
import scipy.sparse as sparse
import somacore
from somacore import options
from somacore.query._eager_iter import EagerIterator
from typing_extensions import Literal

# This package's pybind11 code
import tiledbsoma.pytiledbsoma as clib

from . import _util
from ._exception import SOMAError
from ._types import NTuple

if TYPE_CHECKING:
    from . import SparseNDArray


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


IJDType = Tuple[
    Tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]], npt.NDArray[np.generic]
]


def scipy_sparse_iter(
    sr: clib.SOMAArray,
    coords: options.SparseNDCoords,
    axis: Literal[0, 1],
    step: int,
    compress: bool,
    reindex_sparse_axis: bool,
) -> Iterator[
    Tuple[
        Tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]],
        Union[sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix],
    ]
]:
    """
    Private.

    Iterator over SparseNDArray producing sequence of scipy sparse matrix.
    """
    if len(sr.shape) != 2 or len(coords) > 2 or axis not in [0, 1]:
        raise SOMAError("SciPy sparse matrix iterator compatible only with 2D arrays")

    shape = cast(Tuple[int, int], tuple(sr.shape))

    # TODO: need to add error checks for result_order once we can access it from the
    # SOMAArray reader.

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

    # our thread pool
    pool: Optional[concurrent.futures.Executor] = None

    # This is implemented as a set of stacked generators, to make the code a bit more
    # readable, and to enable injection of parallelism where it helps.

    def _reindexed_tbl_reader() -> Iterator[Tuple[npt.NDArray[np.int64], pa.Table]]:
        """Read stepped tables and reindex them"""
        for major_coords, coo_tbl in EagerIterator(
            _stepped_tbl_reader(sr, coords, minor_coords, step, major_axis), pool=pool
        ):
            ij = [
                coo_tbl.column(0).to_numpy(),  # copies if multi-chunk
                coo_tbl.column(1).to_numpy(),
            ]  # copies if multi-chunk
            d = coo_tbl.column(2).to_numpy()  # copies if multi-chunk

            # Reindex major axis
            ij[major_axis] = pd.Index(major_coords).get_indexer(ij[major_axis])  # type: ignore[no-untyped-call]

            # Optionally, reindex minor axis
            if minor_axis_indexer is not None:
                ij[minor_axis] = minor_axis_indexer.get_indexer(ij[minor_axis])  # type: ignore[no-untyped-call]

            yield major_coords, pa.Table.from_pydict(
                {"soma_dim_0": ij[0], "soma_dim_1": ij[1], "soma_data": d}
            )

    def _sorted_tbl_reader() -> Iterator[Tuple[npt.NDArray[np.int64], IJDType]]:
        """Read reindexed tables and sort them. Yield as ((i,j),d)"""
        for major_coords, coo_tbl in EagerIterator(_reindexed_tbl_reader(), pool=pool):
            coo_tbl = coo_tbl.sort_by(
                [
                    (f"soma_dim_{major_axis}", "ascending"),
                    (f"soma_dim_{minor_axis}", "ascending"),
                ]
            )
            yield major_coords, (
                (coo_tbl.column(0).to_numpy(), coo_tbl.column(1).to_numpy()),
                coo_tbl.column(2).to_numpy(),
            )

    def _mk_shape(
        major_coords: npt.NDArray[np.int64], minor_coords: npt.NDArray[np.int64]
    ) -> Tuple[int, int]:
        """Make shape of this iterator step"""
        assert len(shape) == 2
        _sp_shape: List[int] = list(shape)
        _sp_shape[major_axis] = len(major_coords)
        if reindex_sparse_axis:
            _sp_shape[minor_axis] = len(minor_coords)
        return cast(Tuple[int, int], tuple(_sp_shape))

    def _coo_reader() -> Iterator[sparse.coo_matrix]:
        """Uncompressed variants"""
        assert not compress
        for major_coords, ((i, j), d) in EagerIterator(_sorted_tbl_reader(), pool=pool):
            sp = sparse.coo_matrix(
                (d, (i, j)), shape=_mk_shape(major_coords, minor_coords)
            )

            # SOMA disallows duplicates
            sp.has_canonical_format = True

            _res = (major_coords, minor_coords)
            yield (_res[major_axis], _res[minor_axis]), sp

    def _cs_reader() -> Iterator[Union[sparse.csr_matrix, sparse.csc_matrix]]:
        """Compressed sparse variants"""
        assert compress
        for major_coords, ((i, j), d) in EagerIterator(_sorted_tbl_reader(), pool=pool):
            cls = sparse.csr_matrix if axis == 0 else sparse.csc_matrix
            sp = cls(
                sparse.coo_matrix(
                    (d, (i, j)), shape=_mk_shape(major_coords, minor_coords)
                )
            )
            _res = (major_coords, minor_coords)
            yield (_res[major_axis], _res[minor_axis]), sp

    with concurrent.futures.ThreadPoolExecutor() as _pool:
        pool = _pool
        yield from _cs_reader() if compress else _coo_reader()


def _arrow_table_reader(sr: clib.SOMAArray) -> Iterator[pa.Table]:
    tbl = sr.read_next()
    while tbl is not None:
        yield tbl
        tbl = sr.read_next()


def sparse_stepped_table_iter(
    array: "SparseNDArray",
    sr: clib.SOMAArray,
    coords: options.SparseNDCoords,
    axis: int,
    step: int,
    result_order: Optional[options.ResultOrderStr] = None,
) -> Iterator[pa.Table]:
    """Stepped iterator for Arrow Table"""
    ndim = len(sr.shape)
    if axis < 0 or axis >= ndim:
        raise ValueError("Axis argument must specify a dimension value.")

    # TODO: DRY this. Currently a copy/paste from the Array base class.
    kwargs: Dict[str, object] = {}
    if result_order:
        result_order_map = {
            "auto": clib.ResultOrder.automatic,
            "row-major": clib.ResultOrder.rowmajor,
            "column-major": clib.ResultOrder.colmajor,
        }
        result_order_enum = result_order_map[options.ResultOrder(result_order).value]
        kwargs["result_order"] = result_order_enum

    coords = tuple(coords[i] if i < len(coords) else None for i in range(ndim))
    for coord_chunk in _coords_strider(coords[axis], sr.shape[axis], step):
        sr.reset(**kwargs)
        step_coords = [*coords]
        step_coords[axis] = coord_chunk
        array._set_reader_coords(sr, step_coords)
        tbl = pa.concat_tables(_arrow_table_reader(sr))
        if not len(tbl):
            continue
        yield tbl


def _coords_strider(
    coords: options.SparseNDCoord, length: int, stride: int
) -> Iterator[npt.NDArray[np.int64]]:
    """
    Private.

    Iterate over major coordinates, in stride sized steps, materializing each step as an
    ndarray of coordinate values. Will be sorted in ascending order.

    NB: SOMA slices are _closed_ (i.e., inclusive of both range start and stop)
    """
    # normalize coord to either a slice or ndarray
    if coords is None or (isinstance(coords, slice) and coords == slice(None)):
        coords = slice(0, length - 1)
    elif isinstance(coords, slice):
        if coords == slice(None):
            coords = slice(0, length - 1)
    elif isinstance(coords, Sequence):
        coords = np.array(coords).astype(np.int64)
    elif isinstance(coords, (pa.Array, pa.ChunkedArray)):
        coords = coords.to_numpy()
    elif isinstance(coords, int):
        coords = np.array([coords], dtype=np.int64)
    elif not isinstance(coords, np.ndarray):
        raise TypeError("Unsupported slice coordinate type")

    if isinstance(coords, slice):
        _util.validate_slice(coords)  # NB: this enforces step == 1, assumed below
        start, stop, _step = coords.indices(length - 1)
        assert _step == 1
        yield from (
            np.arange(i, min(i + stride, stop + 1), dtype=np.int64)
            for i in range(start, stop + 1, stride)
        )

    else:
        assert isinstance(coords, np.ndarray) and coords.dtype == np.int64
        for i in range(0, len(coords), stride):
            yield cast(npt.NDArray[np.int64], coords[i : i + stride])


def _stepped_tbl_reader(
    sr: clib.SOMAArray,
    coords: options.SparseNDCoords,
    minor_coords: npt.NDArray[np.int64],
    step: int,
    major_axis: int,
) -> Iterator[Tuple[npt.NDArray[np.int64], pa.Table]]:
    """
    Private. Yields step-sized tables with complete read over major axis.
    """
    minor_axis: int = 1 - major_axis

    for coord_chunk in _coords_strider(coords[major_axis], sr.shape[major_axis], step):
        sr.reset()
        sr.set_dim_points_int64(f"soma_dim_{major_axis}", coord_chunk)
        sr.set_dim_points_int64(f"soma_dim_{minor_axis}", minor_coords)
        tbl = pa.concat_tables(_arrow_table_reader(sr))
        yield coord_chunk, tbl
