# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
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
    Any,
    Iterator,
    List,
    Sequence,
    Tuple,
    TypeVar,
    Union,
    cast,
)

import attrs
import numpy as np
import numpy.typing as npt
import pyarrow as pa
import somacore
from scipy import sparse
from somacore import CoordinateSpace, options

# This package's pybind11 code
import tiledbsoma.pytiledbsoma as clib

from . import _util
from ._eager_iter import EagerIterator
from ._exception import SOMAError
from ._fastercsx import CompressedMatrix
from ._indexer import IntIndexer
from ._query_condition import QueryCondition
from ._types import NTuple
from .options import SOMATileDBContext

if TYPE_CHECKING:
    from . import SparseNDArray
    from ._soma_array import SOMAArray


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
    Tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]],
    npt.NDArray[Union[np.integer[Any], np.floating[Any]]],
]


class TableReadIter(somacore.ReadIter[pa.Table]):
    """Iterator over `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_ elements"""

    def __init__(
        self,
        array: SOMAArray,
        coords: Union[
            options.SparseDFCoords, options.SparseNDCoords, options.DenseNDCoords
        ],
        column_names: Sequence[str] | None,
        result_order: clib.ResultOrder,
        value_filter: str | None,
        platform_config: options.PlatformConfig | None,
        *,
        coord_space: CoordinateSpace | None = None,
    ):
        """Initalizes a new TableReadIter for SOMAArrays.

        Args:
            array (SOMAArray):
                The NDArray, DataFrame, or SpatialDataFrame being read.

            coords (Union[
                options.SparseDFCoords, options.SparseNDCoords, options.DenseNDCoords
            ]):
                for each index dimension, which rows to read.
                ``()`` means no constraint -- all IDs.

            column_names (Sequence[str] | None):
                The named columns to read and return.
                ``None`` means no constraint -- all column names.

            result_order (clib.ResultOrder):
                Order of read results.
                This can be one of automatic, rowmajor, or colmajor.

            value_filter (str | None):
                An optional [value filter] to apply to the results.

            platform_config (options.PlatformConfig | None):
                Pass in parameters for tuning reads.

        """
        self._reader = ArrowTableRead(
            array,
            coords,
            column_names,
            result_order,
            value_filter,
            platform_config,
            coord_space=coord_space,
        )

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
        array: SOMAArray,
        coords: options.SparseNDCoords,
        axis: Union[int, Sequence[int]],
        result_order: clib.ResultOrder,
        platform_config: options.PlatformConfig | None,
        *,
        size: int | Sequence[int] | None = None,
        reindex_disable_on_axis: int | Sequence[int] | None = None,
        eager: bool = True,
        context: SOMATileDBContext | None = None,
    ):
        super().__init__()

        self.array = array
        self.shape = array._handle._handle.shape
        self.ndim = len(self.shape)
        self.eager = eager
        self.result_order = result_order
        self.platform_config = platform_config

        # Assign a thread pool from the context, or create a new one if no context
        # is available
        if context is not None:
            self._threadpool = context.threadpool
        else:
            self._threadpool = futures.ThreadPoolExecutor()

        self.context = context

        # raises on various error checks, AND normalizes args
        self.axis, self.size, self.reindex_disable_on_axis = self._validate_args(
            self.shape, axis, size, reindex_disable_on_axis
        )

        self.major_axis = self.axis[0]
        self.coords = _pad_with_none(coords, self.ndim)

        # materialize all indexing info.
        self.joinids: List[pa.Array] = [
            pa.array(
                np.concatenate(
                    list(_coords_strider(self.coords[d], self.shape[d], self.shape[d]))
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
        size: int | Sequence[int] | None = None,
        reindex_disable_on_axis: int | Sequence[int] | None = None,
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
        self, x: Iterator[_EagerRT], _pool: ThreadPoolExecutor | None = None
    ) -> Iterator[_EagerRT]:
        """Private"""
        return EagerIterator(x, pool=_pool) if self.eager else x

    def _table_reader(self) -> Iterator[BlockwiseTableReadIterResult]:
        """Private. Blockwise table reader. Helper function for sub-class use"""
        for coord_chunk in _coords_strider(
            self.coords[self.major_axis],
            self.shape[self.major_axis],
            self.size[0],
        ):
            step_coords = list(self.coords)
            step_coords[self.major_axis] = coord_chunk

            joinids = list(self.joinids)
            joinids[self.major_axis] = pa.array(coord_chunk)
            yield pa.concat_tables(
                ArrowTableRead(
                    array=self.array,
                    coords=step_coords,
                    column_names=[],  # select all columns
                    result_order=self.result_order,
                    value_filter=None,
                    platform_config=self.platform_config,
                )
            ), tuple(joinids)

    def _reindexed_table_reader(
        self,
        _pool: ThreadPoolExecutor | None = None,
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
        array: SOMAArray,
        coords: options.SparseNDCoords,
        axis: Union[int, Sequence[int]],
        result_order: clib.ResultOrder,
        platform_config: options.PlatformConfig | None,
        *,
        size: int | Sequence[int] | None = None,
        reindex_disable_on_axis: int | Sequence[int] | None = None,
        eager: bool = True,
        compress: bool = True,
        context: SOMATileDBContext | None = None,
    ):
        self.compress = compress
        self.context = context
        super().__init__(
            array,
            coords,
            axis,
            result_order,
            platform_config,
            size=size,
            reindex_disable_on_axis=reindex_disable_on_axis,
            eager=eager,
            context=context,
        )

        if (
            len(self.shape) != 2
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
        self, _pool: ThreadPoolExecutor | None = None
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
        shape = cast(Tuple[int, int], tuple(self.shape))
        assert len(shape) == 2
        _sp_shape: List[int] = list(shape)

        if self.major_axis not in self.reindex_disable_on_axis:
            _sp_shape[self.major_axis] = len(major_coords)
        if self.minor_axis not in self.reindex_disable_on_axis:
            _sp_shape[self.minor_axis] = len(minor_coords)

        return cast(Tuple[int, int], tuple(_sp_shape))

    def _coo_reader(
        self, _pool: ThreadPoolExecutor | None = None
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
            if self.result_order == clib.ResultOrder.rowmajor:
                sp.has_canonical_format = True

            yield sp, indices

    def _cs_reader(
        self, _pool: ThreadPoolExecutor | None = None
    ) -> Iterator[Tuple[Union[sparse.csr_matrix, sparse.csc_matrix], IndicesType],]:
        """Private. Compressed sparse variants"""
        assert self.compress
        assert self.major_axis not in self.reindex_disable_on_axis
        for ((i, j), d), indices in self._maybe_eager_iterator(
            self._sorted_tbl_reader(_pool), _pool
        ):
            major_coords = indices[self.major_axis]
            minor_coords = indices[self.minor_axis]
            shape = self._mk_shape(major_coords, minor_coords)
            assert self.context is not None
            sp = CompressedMatrix.from_ijd(
                i,
                j,
                d,
                shape=shape,
                format="csr" if self.major_axis == 0 else "csc",
                make_sorted=True,
                context=self.context,
            ).to_scipy()
            yield sp, indices


@attrs.define(frozen=True)
class ManagedQuery:
    """Keep the lifetime of the SOMAArray tethered to ManagedQuery."""

    _array: SOMAArray
    _platform_config: options.PlatformConfig | None = None
    _handle: clib.ManagedQuery = attrs.field(init=False)

    def __attrs_post_init__(self) -> None:
        clib_handle = self._array._handle._handle

        if self._platform_config is not None:
            cfg = clib_handle.context().config()
            cfg.update(self._platform_config)
            ctx = clib.SOMAContext(cfg)
        else:
            ctx = clib_handle.context()

        object.__setattr__(self, "_handle", clib.ManagedQuery(clib_handle, ctx))


class SparseTensorReadIterBase(somacore.ReadIter[_RT], metaclass=abc.ABCMeta):
    """Private implementation class"""

    def __init__(
        self,
        array: SparseNDArray,
        coords: options.SparseDFCoords,
        shape: NTuple,
        result_order: clib.ResultOrder,
        platform_config: options.PlatformConfig | None,
    ):
        self.array = array
        self.coords = coords
        self.shape = shape
        self.result_order = result_order
        self.platform_config = platform_config

        self.mq = ManagedQuery(array, platform_config)

        self.mq._handle.set_layout(result_order)

        _util._set_coords(self.mq, coords)

    @abc.abstractmethod
    def _from_table(self, arrow_table: pa.Table) -> _RT:
        raise NotImplementedError()

    def __next__(self) -> _RT:
        return self._from_table(self.mq._handle.next())

    def concat(self) -> _RT:
        """Returns all the requested data in a single operation.

        If some data has already been retrieved using ``next``, this will return
        the rest of the data after that is already returned.
        """
        arrow_tables = pa.concat_tables(
            TableReadIter(
                array=self.array,
                coords=self.coords,
                column_names=[],  # select all columns
                result_order=self.result_order,
                value_filter=None,
                platform_config=self.platform_config,
            )
        )
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


class ArrowTableRead(Iterator[pa.Table]):
    def __init__(
        self,
        array: SOMAArray,
        coords: Union[
            options.SparseDFCoords, options.SparseNDCoords, options.DenseNDCoords
        ],
        column_names: Sequence[str] | None,
        result_order: clib.ResultOrder,
        value_filter: str | None,
        platform_config: options.PlatformConfig | None,
        *,
        coord_space: CoordinateSpace | None = None,
    ):
        clib_handle = array._handle._handle

        self.mq = ManagedQuery(array, platform_config)

        self.mq._handle.set_layout(result_order)

        column_names = column_names or array.schema.names
        for name in column_names:
            clib_handle.get_column(name).select_columns(self.mq._handle)

        if value_filter is not None:
            self.mq._handle.set_condition(
                QueryCondition(value_filter), clib_handle.schema
            )

        _util._set_coords(
            self.mq, coords, coord_space.axis_names if coord_space else None
        )

    def __next__(self) -> pa.Table:
        return self.mq._handle.next()


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


def _pad_with_none(s: Sequence[_ElemT], to_length: int) -> Tuple[_ElemT | None, ...]:
    """Given a sequence, pad length to a user-specified length, with None values"""
    return tuple(s[i] if i < len(s) else None for i in range(to_length))
