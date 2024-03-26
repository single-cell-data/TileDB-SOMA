# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Read iterators.
"""
from __future__ import annotations

import abc
from concurrent import futures
from concurrent.futures import ThreadPoolExecutor
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
import pyarrow as pa
import somacore
from scipy import sparse
from somacore import options
from somacore.query._eager_iter import EagerIterator

# This package's pybind11 code
import tiledbsoma.pytiledbsoma as clib

from . import _util
from ._exception import SOMAError
from ._indexer import IntIndexer
from ._types import NTuple
from .options import SOMATileDBContext

if TYPE_CHECKING:
    from . import SparseNDArray


# Convenience types
_RT = TypeVar("_RT")
BlockwiseTableReadIterResult = Tuple[pa.Table, Tuple[pa.Array, ...]]
BlockwiseSingleAxisTableIter = Iterator[BlockwiseTableReadIterResult]

BlockwiseScipyReadIterResult = Tuple[
    Union[sparse.csr_matrix, sparse.csc_matrix, sparse.coo_matrix],
    Tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]],
]

IndicesType = Tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]
IJDType = Tuple[
    Tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]], npt.NDArray[np.generic]
]


class TableReadIter(somacore.ReadIter[pa.Table]):
    """Iterator over `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_ elements"""

    def __init__(self, sr: clib.SOMAArray):
        self._reader = _arrow_table_reader(sr)

    def __next__(self) -> pa.Table:
        return next(self._reader)

    def concat(self) -> pa.Table:
        """Concatenate remainder of iterator, and return as a single `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_"""
        return pa.concat_tables(self)


_EagerRT = TypeVar("_EagerRT")


class BlockwiseReadIterBase(somacore.ReadIter[_RT], metaclass=abc.ABCMeta):
    """Private implementation class. Currently implemented as a single-axis blockwise iterator"""

    _reader: Iterator[_RT]

    def __init__(
        self,
        array: "SparseNDArray",
        sr: clib.SOMAArray,
        coords: options.SparseNDCoords,
        axis: Union[int, Sequence[int]],
        *,
        size: Optional[Union[int, Sequence[int]]] = None,
        reindex_disable_on_axis: Optional[Union[int, Sequence[int]]] = None,
        eager: bool = True,
        context: Optional[SOMATileDBContext] = None,
    ):
        super().__init__()

        self.ndim = len(sr.shape)
        self.array = array
        self.sr = sr
        self.eager = eager

        # Assign a thread pool from the context, or create a new one if no context
        # is available
        if context is not None:
            self._threadpool = context.threadpool
        else:
            self._threadpool = futures.ThreadPoolExecutor()

        self.context = context

        # raises on various error checks, AND normalizes args
        self.axis, self.size, self.reindex_disable_on_axis = self._validate_args(
            sr.shape, axis, size, reindex_disable_on_axis
        )

        self.major_axis = self.axis[0]
        self.coords = _pad_with_none(coords, self.ndim)

        # materialize all indexing info.
        self.joinids: List[pa.Array] = [
            pa.array(
                np.concatenate(
                    list(
                        _coords_strider(
                            self.coords[d], self.sr.shape[d], self.sr.shape[d]
                        )
                    )
                )
                if d != self.major_axis
                else np.array([], dtype=np.int64)
            )
            for d in range(self.ndim)
        ]

        # build indexers, as needed
        self.axes_to_reindex = set(range(self.ndim)) - set(self.reindex_disable_on_axis)
        assert context is not None
        self.minor_axes_indexer = {
            d: IntIndexer(self.joinids[d].to_numpy(), context=context)
            for d in (self.axes_to_reindex - set((self.major_axis,)))
        }

        # Ask subclass to create the type-specific reader/iterator
        self._reader = self._create_reader()

    @classmethod
    def _validate_args(
        cls,
        shape: Union[NTuple, Sequence[int]],
        axis: Union[int, Sequence[int]],
        size: Optional[Union[int, Sequence[int]]] = None,
        reindex_disable_on_axis: Optional[Union[int, Sequence[int]]] = None,
    ) -> Tuple[List[int], List[int], List[int]]:
        """
        Class method to validate and normalize common user-provided arguments axis, size and reindex_disable_on_axis.
        * normalize args to fully specify the arg per dimension, for convenience in later usage
        * error check and raise if a nonsense value
        * set defaults where supported.
        """
        ndim = len(shape)

        # convert to list
        axis = list(axis if isinstance(axis, Sequence) else [axis])

        if reindex_disable_on_axis is None:
            reindex_disable_on_axis = []
        elif isinstance(reindex_disable_on_axis, int):
            reindex_disable_on_axis = [reindex_disable_on_axis]
        elif isinstance(reindex_disable_on_axis, Sequence):
            reindex_disable_on_axis = list(reindex_disable_on_axis)
        else:
            raise TypeError(
                "reindex_disable_on_axis must be None, int or Sequence[int]"
            )

        # Currently, only support blockwise iteration on one dimension.
        if len(axis) != 1:
            raise NotImplementedError(
                "Multi-dimension blockwise iterators not implemented"
            )
        # all dim indices must be in acceptable range
        if not all(0 <= d < ndim for d in axis):
            raise ValueError("blockwise `axis` value must be in range [0, ndim)")
        if not all(0 <= d < ndim for d in reindex_disable_on_axis):
            raise ValueError(
                "blockwise `reindex_disable_on_axis` value must be in range [0, ndim)"
            )

        # if not specified, set default size for each axis. Default heuristic
        # assumes 2D array has many more rows than cols (i.e., n_obs>>n_vars).
        default_block_size = (2**16,) + (2**8,) * (ndim - 1)
        if size is None:
            size = [default_block_size[d] for d in axis]
        elif isinstance(size, int):
            size = [size] * len(axis)
        elif isinstance(size, Sequence):
            size = list(size) + [default_block_size[d] for d in axis[len(size) :]]
        else:
            raise TypeError(
                "blockwise iterator `size` must be None, int or Sequence[int]"
            )

        return axis, size, reindex_disable_on_axis

    @abc.abstractmethod
    def _create_reader(self) -> Iterator[_RT]:
        """Sub-class responsibility"""
        raise NotImplementedError()

    def __next__(self) -> _RT:
        return next(self._reader)

    def concat(self) -> _RT:
        """
        Unimplemented as there is little utility beyond that offered by a ragged
        read iterator concat, other than reindexing.
        """
        raise NotImplementedError("Blockwise iterators do not support concat operation")

    def _maybe_eager_iterator(
        self, x: Iterator[_EagerRT], _pool: Optional[ThreadPoolExecutor] = None
    ) -> Iterator[_EagerRT]:
        """Private"""
        return EagerIterator(x, pool=_pool) if self.eager else x

    def _table_reader(self) -> Iterator[BlockwiseTableReadIterResult]:
        """Private. Blockwise table reader. Helper function for sub-class use"""
        kwargs: Dict[str, object] = {"result_order": self.sr.result_order}
        for coord_chunk in _coords_strider(
            self.coords[self.major_axis],
            self.sr.shape[self.major_axis],
            self.size[0],
        ):
            self.sr.reset(**kwargs)
            step_coords = list(self.coords)
            step_coords[self.major_axis] = coord_chunk
            self.array._set_reader_coords(self.sr, step_coords)

            joinids = list(self.joinids)
            joinids[self.major_axis] = pa.array(coord_chunk)
            yield pa.concat_tables(_arrow_table_reader(self.sr)), tuple(joinids)

    def _reindexed_table_reader(
        self,
        _pool: Optional[ThreadPoolExecutor] = None,
    ) -> Iterator[BlockwiseTableReadIterResult]:
        """Private. Blockwise table reader w/ reindexing. Helper function for sub-class use"""
        for tbl, coords in self._maybe_eager_iterator(self._table_reader(), _pool):
            pytbl = {}
            for d in range(self.ndim):
                col = tbl.column(f"soma_dim_{d}")
                if d in self.axes_to_reindex:
                    if d == self.major_axis:
                        assert self.context is not None
                        col = IntIndexer(
                            coords[self.major_axis], context=self.context
                        ).get_indexer(
                            col.to_numpy(),
                        )
                    else:
                        col = self.minor_axes_indexer[d].get_indexer(col.to_numpy())
                pytbl[f"soma_dim_{d}"] = col
            pytbl["soma_data"] = tbl.column("soma_data")
            yield pa.Table.from_pydict(pytbl), coords


class BlockwiseTableReadIter(BlockwiseReadIterBase[BlockwiseTableReadIterResult]):
    """Blockwise iterator over `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_ elements"""

    def _create_reader(self) -> Iterator[BlockwiseTableReadIterResult]:
        """Private. Blockwise Arrow Table iterator, restricted to a single axis"""
        yield from (
            self._reindexed_table_reader(_pool=self._threadpool)
            if self.axes_to_reindex
            else self._table_reader()
        )


class BlockwiseScipyReadIter(BlockwiseReadIterBase[BlockwiseScipyReadIterResult]):
    """Blockwise iterator over `SciPy sparse matrix <https://docs.scipy.org/doc/scipy/reference/sparse.html>`_ elements"""

    def __init__(
        self,
        array: "SparseNDArray",
        sr: clib.SOMAArray,
        coords: options.SparseNDCoords,
        axis: Union[int, Sequence[int]],
        *,
        size: Optional[Union[int, Sequence[int]]] = None,
        reindex_disable_on_axis: Optional[Union[int, Sequence[int]]] = None,
        eager: bool = True,
        compress: bool = True,
        context: Optional[SOMATileDBContext] = None,
    ):
        self.compress = compress
        self.context = context
        super().__init__(
            array,
            sr,
            coords,
            axis,
            size=size,
            reindex_disable_on_axis=reindex_disable_on_axis,
            eager=eager,
            context=context,
        )

        if (
            len(self.sr.shape) != 2
            or len(self.coords) > 2
            or self.major_axis not in [0, 1]
        ):
            raise SOMAError(
                "SciPy sparse matrix iterator compatible only with 2D arrays"
            )

        if self.compress and self.major_axis in self.reindex_disable_on_axis:
            raise SOMAError(
                "Unable to disable reindexing of coordinates on CSC/CSR major axis"
            )

        # Sanity check: if CSC/CSR, we _must_ be reindexing
        assert not self.compress or self.axes_to_reindex

    @property
    def minor_axis(self) -> int:
        """For convienence's sake"""
        return 1 - self.major_axis

    def _create_reader(self) -> Iterator[BlockwiseScipyReadIterResult]:
        """
        Private. Iterator over SparseNDArray producing sequence of scipy sparse matrix.
        """
        yield from (
            self._cs_reader(_pool=self._threadpool)
            if self.compress
            else self._coo_reader(_pool=self._threadpool)
        )

    def _sorted_tbl_reader(
        self, _pool: Optional[ThreadPoolExecutor] = None
    ) -> Iterator[Tuple[IJDType, IndicesType]]:
        """Private. Read reindexed tables and sort them. Yield as ((i,j),d)"""
        for coo_tbl, indices in self._maybe_eager_iterator(
            self._reindexed_table_reader(_pool), _pool
        ):
            coo_tbl = coo_tbl.sort_by(
                [
                    (f"soma_dim_{self.major_axis}", "ascending"),
                    (f"soma_dim_{self.minor_axis}", "ascending"),
                ]
            )
            ijd = (
                (coo_tbl.column(0).to_numpy(), coo_tbl.column(1).to_numpy()),
                coo_tbl.column(2).to_numpy(),
            )
            yield ijd, (indices[0].to_numpy(), indices[1].to_numpy())

    def _mk_shape(
        self, major_coords: npt.NDArray[np.int64], minor_coords: npt.NDArray[np.int64]
    ) -> Tuple[int, int]:
        """Private. Make shape of this iterator step"""
        shape = cast(Tuple[int, int], tuple(self.sr.shape))
        assert len(shape) == 2
        _sp_shape: List[int] = list(shape)

        if self.major_axis not in self.reindex_disable_on_axis:
            _sp_shape[self.major_axis] = len(major_coords)
        if self.minor_axis not in self.reindex_disable_on_axis:
            _sp_shape[self.minor_axis] = len(minor_coords)

        return cast(Tuple[int, int], tuple(_sp_shape))

    def _coo_reader(
        self, _pool: Optional[ThreadPoolExecutor] = None
    ) -> Iterator[Tuple[sparse.coo_matrix, IndicesType]]:
        """Private. Uncompressed variants"""
        assert not self.compress
        for ((i, j), d), indices in self._maybe_eager_iterator(
            self._sorted_tbl_reader(_pool), _pool
        ):
            major_coords, minor_coords = (
                indices[self.major_axis],
                indices[self.minor_axis],
            )
            sp = sparse.coo_matrix(
                (d, (i, j)), shape=self._mk_shape(major_coords, minor_coords)
            )

            # SOMA disallows duplicates. Canonical implies sorted row-major, no dups
            if self.sr.result_order == clib.ResultOrder.rowmajor:
                sp.has_canonical_format = True

            yield sp, indices

    def _cs_reader(
        self, _pool: Optional[ThreadPoolExecutor] = None
    ) -> Iterator[Tuple[Union[sparse.csr_matrix, sparse.csc_matrix], IndicesType],]:
        """Private. Compressed sparse variants"""
        assert self.compress
        assert self.major_axis not in self.reindex_disable_on_axis
        for ((i, j), d), indices in self._maybe_eager_iterator(
            self._sorted_tbl_reader(_pool), _pool
        ):
            major_coords = indices[self.major_axis]
            minor_coords = indices[self.minor_axis]
            cls = sparse.csr_matrix if self.major_axis == 0 else sparse.csc_matrix
            sp = cls(
                sparse.coo_matrix(
                    (d, (i, j)), shape=self._mk_shape(major_coords, minor_coords)
                )
            )
            yield sp, indices


class SparseTensorReadIterBase(somacore.ReadIter[_RT], metaclass=abc.ABCMeta):
    """Private implementation class"""

    def __init__(self, sr: clib.SOMAArray, shape: NTuple):
        self.sr = sr
        self.shape = shape

    @abc.abstractmethod
    def _from_table(self, arrow_table: pa.Table) -> _RT:
        raise NotImplementedError()

    def __next__(self) -> _RT:
        arrow_table = self.sr.read_next()
        if arrow_table is None:
            raise StopIteration

        return self._from_table(arrow_table)

    def concat(self) -> _RT:
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


def _arrow_table_reader(sr: clib.SOMAArray) -> Iterator[pa.Table]:
    """Private. Simple Table iterator on any Array"""
    tbl = sr.read_next()
    while tbl is not None:
        yield tbl
        tbl = sr.read_next()


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

    # NB: type check on slice is to handle the case where coords is an NDArray,
    # and the equality check is broadcast to all elements of the array.
    if coords is None or (isinstance(coords, slice) and coords == slice(None)):
        coords = slice(0, length - 1)
    elif isinstance(coords, int):
        coords = np.array([coords], dtype=np.int64)
    elif isinstance(coords, Sequence):
        coords = np.array(coords, dtype=np.int64)
    elif isinstance(coords, (pa.Array, pa.ChunkedArray)):
        coords = coords.to_numpy().astype(np.int64, copy=False)
    elif isinstance(coords, np.ndarray):
        coords = coords.astype(np.int64)
    elif isinstance(coords, slice):
        pass
    else:
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


_ElemT = TypeVar("_ElemT")


def _pad_with_none(s: Sequence[_ElemT], to_length: int) -> Tuple[Optional[_ElemT], ...]:
    """Given a sequence, pad length to a user-specified length, with None values"""
    return tuple(s[i] if i < len(s) else None for i in range(to_length))
