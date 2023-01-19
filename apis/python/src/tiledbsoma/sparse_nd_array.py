import collections.abc
from typing import Any, List, Optional, Union, cast

import numpy as np
import pyarrow as pa
import somacore
import tiledb
from somacore import options
from somacore.options import PlatformConfig

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from . import util, util_arrow
from .collection import CollectionBase
from .options import SOMATileDBContext, TileDBCreateOptions
from .tiledb_array import TileDBArray
from .types import NTuple
from .util_iter import (
    SparseCOOTensorReadIter,
    SparseCSCMatrixReadIter,
    SparseCSRMatrixReadIter,
    TableReadIter,
)

_UNBATCHED = options.BatchSize()


class SparseNDArray(TileDBArray, somacore.SparseNDArray):
    """
    Represents ``X`` and others.
    """

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[CollectionBase[Any]] = None,
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """

        super().__init__(uri=uri, parent=parent, context=context)

    # Inherited from somacore
    # soma_type: Final = "SOMASparseNDArray"

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
        platform_config: Optional[PlatformConfig] = None,
    ) -> "SparseNDArray":
        """
        Create a ``SparseNDArray`` named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the positive int64 range.

        :param platform_config: Platform-specific options used to create this Array, provided via "tiledb"->"create" nested keys
        """

        # check on shape
        if len(shape) == 0 or any(e <= 0 for e in shape):
            raise ValueError(
                "DenseNDArray shape must be non-zero length tuple of ints > 0"
            )

        if not pa.types.is_primitive(type):
            raise TypeError(
                "Unsupported type - DenseNDArray only supports primtive Arrow types"
            )

        tiledb_create_options = TileDBCreateOptions.from_platform_config(
            platform_config
        )

        dims = []
        for n, e in enumerate(shape):
            dim_name = f"soma_dim_{n}"
            dim = tiledb.Dim(
                name=dim_name,
                domain=(0, e - 1),
                tile=tiledb_create_options.dim_tile(dim_name, min(e, 2048)),
                dtype=np.int64,
                filters=tiledb_create_options.dim_filters(
                    dim_name,
                    [
                        dict(
                            _type="ZstdFilter",
                            level=tiledb_create_options.string_dim_zstd_level(),
                        )
                    ],
                ),
            )
            dims.append(dim)
        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = [
            tiledb.Attr(
                name="soma_data",
                dtype=util_arrow.tiledb_type_from_arrow_type(type),
                filters=tiledb_create_options.attr_filters("soma_data", ["ZstdFilter"]),
                ctx=self._ctx,
            )
        ]

        cell_order, tile_order = tiledb_create_options.cell_tile_orders()

        # TODO: code-dedupe w/ regard to DenseNDArray. The two creates are
        # almost identical & could share a common parent-class _create() method.
        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=True,
            offsets_filters=tiledb_create_options.offsets_filters(),
            capacity=tiledb_create_options.get("capacity", 100000),
            tile_order=tile_order,
            cell_order=cell_order,
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch)

        self._common_create(self.soma_type)  # object-type metadata etc

        return self

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndim``
        """
        with self._tiledb_open() as A:
            return cast(NTuple, A.schema.domain.shape)

    def reshape(self, shape: NTuple) -> None:
        """
        Unsupported operation for this object type.
        """
        raise NotImplementedError("reshape operation not implemented.")

    # Inherited from somacore
    # * ndim accessor
    # * is_sparse: Final = True

    @property
    def nnz(self) -> int:
        """
        Return the number of stored values in the array, including explicitly stored zeros.
        """
        return cast(
            int,
            clib.SOMAReader(
                self.uri,
                platform_config={} if self._ctx is None else self._ctx.config().dict(),
            ).nnz(),
        )

    def read(
        self,
        coords: Optional[options.SparseNDCoords] = None,
        *,
        result_order: options.StrOr[options.ResultOrder] = options.ResultOrder.AUTO,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[PlatformConfig] = None,
    ) -> "SparseNDArrayRead":
        """
        Read a user-defined slice of the SparseNDArray.

        Parameters
        ----------
        coords : Tuple[Union[int, slice, Tuple[int, ...], List[int], pa.IntegerArray], ...]
            Per-dimension tuple of scalar, slice, sequence of scalar or Arrow IntegerArray
            Arrow arrays currently unimplemented.

        Acceptable ways to index
        ------------------------
        * None
        * A sequence of coordinates is accepted, one per dimension.
        * Sequence length must be at least one and <= number of dimensions.
        * If the sequence contains missing coordinates (length less than number of dimensions),
          then `slice(None)` -- i.e. no constraint -- is assumed for the missing dimensions.
        * Per-dimension, explicitly specified coordinates can be one of: None, a value, a
          list/ndarray/paarray/etc of values, a slice, etc.
        * Slices are doubly inclusive: slice(2,4) means [2,3,4] not [2,3]. Slice steps can only be +1.
          Slices can be `slice(None)`, meaning select all in that dimension, but may not be half-specified:
          `slice(2,None)` and `slice(None,4)` are both unsupported.
        * Negative indexing is unsupported.

        Returns
        -------
        SparseNDArrayRead - which can be used to access an iterator of results in various formats.
        """
        del result_order, batch_size, partitions, platform_config  # Currently unused.

        if coords is None:
            coords = (slice(None),)

        with self._tiledb_open("r") as A:
            shape = A.shape
            sr = clib.SOMAReader(
                self._uri,
                name=self.__class__.__name__,
                schema=A.schema,
                platform_config={} if self._ctx is None else self._ctx.config().dict(),
            )

            if not isinstance(coords, (list, tuple)):
                raise TypeError(
                    f"coords type {type(coords)} unsupported; expected list or tuple"
                )
            if len(coords) < 1 or len(coords) > A.schema.domain.ndim:
                raise ValueError(
                    f"coords {coords} must have length between 1 and ndim ({A.schema.domain.ndim}); got {len(coords)}"
                )

            for i, coord in enumerate(coords):
                #                # Example: coords = [None, 3, slice(4,5)]
                #                # coord takes on values None, 3, and slice(4,5) in this loop body.
                dim_name = A.schema.domain.dim(i).name
                if coord is None:
                    pass  # No constraint; select all in this dimension
                elif isinstance(coord, int):
                    sr.set_dim_points(dim_name, [coord])
                elif isinstance(coord, np.ndarray):
                    if coord.ndim != 1:
                        raise ValueError(
                            f"only 1D numpy arrays may be used to index; got {coord.ndim}"
                        )
                    sr.set_dim_points(dim_name, coord)
                elif isinstance(coord, slice):
                    ned = A.nonempty_domain()  # None iff the array has no data
                    lo_hi = util.slice_to_range(coord, ned[i]) if ned else None
                    if lo_hi is not None:
                        lo, hi = lo_hi
                        if lo < 0 or hi < 0:
                            raise ValueError(
                                f"slice start and stop may not be negative; got ({lo}, {hi})"
                            )
                        if lo > hi:
                            raise ValueError(
                                f"slice start must be <= slice stop; got ({lo}, {hi})"
                            )
                        sr.set_dim_ranges(dim_name, [lo_hi])
                    # Else, no constraint in this slot. This is `slice(None)` which is like
                    # Python indexing syntax `[:]`.
                elif isinstance(
                    coord, (collections.abc.Sequence, pa.Array, pa.ChunkedArray)
                ):
                    sr.set_dim_points(dim_name, coord)
                else:
                    raise TypeError(f"coord type {type(coord)} at slot {i} unsupported")

        sr.submit()
        return SparseNDArrayRead(sr, shape)

    def write(
        self,
        values: Union[
            pa.SparseCOOTensor,
            pa.SparseCSRMatrix,
            pa.SparseCSCMatrix,
            pa.Table,
        ],
        *,
        platform_config: Optional[PlatformConfig] = None,
    ) -> None:
        """
        Write an Arrow object to the SparseNDArray.

        Arrow sparse tensor: the coordinates in the Arrow SparseTensor will be interpreted
        as the coordinates to write to. Supports the _experimental_ SparseCOOTensor,
        SparseCSRMatrix and SparseCSCMatrix. There is currently no support for Arrow
        SparseCSFTensor or dense Tensor.

        Arrow table: write a COO table, with columns named ``soma_dim_0``, ...,
        ``soma_dim_N`` and ``soma_data`` to the dense nD array.
        """
        del platform_config  # Currently unused.
        if isinstance(values, pa.SparseCOOTensor):
            data, coords = values.to_numpy()
            with self._tiledb_open("w") as A:
                A[tuple(c for c in coords.T)] = data
            return

        if isinstance(values, (pa.SparseCSCMatrix, pa.SparseCSRMatrix)):
            if self.ndim != 2:
                raise ValueError(
                    f"Unable to write 2D Arrow sparse matrix to {self.ndim}D SparseNDArray"
                )
            # TODO: the `to_scipy` function is not zero copy. Need to explore zero-copy options.
            sp = values.to_scipy().tocoo()
            with self._tiledb_open("w") as A:
                A[sp.row, sp.col] = sp.data
            return

        if isinstance(values, pa.Table):
            data = values.column("soma_data").to_numpy()
            coord_tbl = values.drop(["soma_data"])
            coords = tuple(
                coord_tbl.column(f"soma_dim_{n}").to_numpy()
                for n in range(coord_tbl.num_columns)
            )
            with self._tiledb_open("w") as A:
                A[coords] = data
            return

        raise TypeError("Unsupported Arrow type")


class SparseNDArrayRead(somacore.SparseRead):
    def __init__(self, sr: clib.SOMAReader, shape: NTuple):
        self.sr = sr
        self.shape = shape

    def coos(self) -> SparseCOOTensorReadIter:
        """Return an iterator of Arrow SparseCOOTensor"""
        return SparseCOOTensorReadIter(self.sr, self.shape)

    def cscs(self) -> SparseCSCMatrixReadIter:
        """Return an iterator of Arrow SparseCSCMatrix"""
        return SparseCSCMatrixReadIter(self.sr, self.shape)

    def csrs(self) -> SparseCSRMatrixReadIter:
        """Return an iterator of Arrow SparseCSRMatrix"""
        return SparseCSRMatrixReadIter(self.sr, self.shape)

    def dense_tensors(self) -> somacore.ReadIter[pa.Tensor]:
        """Return an iterator of Arrow Tensor"""
        raise NotImplementedError()

    def tables(self) -> TableReadIter:
        """Return an iterator of Arrow Table"""
        return TableReadIter(self.sr)
