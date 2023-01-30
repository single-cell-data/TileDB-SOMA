import collections.abc
from typing import Optional, Sequence, Union, cast

import numpy as np
import pyarrow as pa
import somacore
from somacore import options
from somacore.options import PlatformConfig

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from . import util
from .common_nd_array import build_tiledb_schema
from .options import SOMATileDBContext
from .options.tiledb_create_options import TileDBCreateOptions
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

    [lifecycle: experimental]
    """

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        shape: Sequence[int],
        platform_config: Optional[PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
    ) -> "SparseNDArray":
        """
        Create a ``SparseNDArray`` named with the URI.

        [lifecycle: experimental]

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the positive int64 range.

        :param platform_config: Platform-specific options used to create this Array, provided via "tiledb"->"create" nested keys
        """
        context = context or SOMATileDBContext()
        schema = build_tiledb_schema(
            type,
            shape,
            TileDBCreateOptions.from_platform_config(platform_config),
            context,
            is_sparse=True,
        )
        handle = cls._create_internal(uri, schema, context)
        return cls(handle, _this_is_internal_only="tiledbsoma-internal-code")

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndim``
        """
        return cast(NTuple, self._handle.reader.schema.domain.shape)

    def reshape(self, shape: NTuple) -> None:
        """
        Unsupported operation for this object type.

        [lifecycle: experimental]
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
        self._check_open_read()
        return cast(int, self._soma_reader().nnz())

    def read(
        self,
        coords: Optional[options.SparseNDCoords] = None,
        *,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[PlatformConfig] = None,
    ) -> "SparseNDArrayRead":
        """
        Read a user-defined slice of the SparseNDArray.

        [lifecycle: experimental]

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
        self._check_open_read()

        if coords is None:
            coords = (slice(None),)

        arr = self._handle.reader
        shape = arr.shape
        sr = self._soma_reader(schema=arr.schema)

        if not isinstance(coords, (list, tuple)):
            raise TypeError(
                f"coords type {type(coords)} unsupported; expected list or tuple"
            )
        if len(coords) < 1 or len(coords) > arr.schema.domain.ndim:
            raise ValueError(
                f"coords {coords} must have length between 1 and ndim ({arr.schema.domain.ndim}); got {len(coords)}"
            )

        for i, coord in enumerate(coords):
            #                # Example: coords = [None, 3, slice(4,5)]
            #                # coord takes on values None, 3, and slice(4,5) in this loop body.
            dim_name = arr.schema.domain.dim(i).name
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
                ned = arr.nonempty_domain()  # None iff the array has no data
                lo_hi = util.slice_to_range(coord, ned[i]) if ned else None
                if lo_hi is not None:
                    lo, hi = lo_hi
                    if lo < 0 or hi < 0:
                        raise ValueError(
                            f"slice start and stop may not be negative; got ({lo}, {hi})"
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

        [lifecycle: experimental]

        Arrow sparse tensor: the coordinates in the Arrow SparseTensor will be interpreted
        as the coordinates to write to. Supports the _experimental_ SparseCOOTensor,
        SparseCSRMatrix and SparseCSCMatrix. There is currently no support for Arrow
        SparseCSFTensor or dense Tensor.

        Arrow table: write a COO table, with columns named ``soma_dim_0``, ...,
        ``soma_dim_N`` and ``soma_data`` to the dense nD array.
        """
        del platform_config  # Currently unused.

        arr = self._handle.writer

        if isinstance(values, pa.SparseCOOTensor):
            data, coords = values.to_numpy()
            arr[tuple(c for c in coords.T)] = data
            return

        if isinstance(values, (pa.SparseCSCMatrix, pa.SparseCSRMatrix)):
            if self.ndim != 2:
                raise ValueError(
                    f"Unable to write 2D Arrow sparse matrix to {self.ndim}D SparseNDArray"
                )
            # TODO: the `to_scipy` function is not zero copy. Need to explore zero-copy options.
            sp = values.to_scipy().tocoo()
            arr[sp.row, sp.col] = sp.data
            return

        if isinstance(values, pa.Table):
            data = values.column("soma_data").to_numpy()
            coord_tbl = values.drop(["soma_data"])
            coords = tuple(
                coord_tbl.column(f"soma_dim_{n}").to_numpy()
                for n in range(coord_tbl.num_columns)
            )
            arr[coords] = data
            return

        raise TypeError(
            f"Unsupported Arrow type or non-arrow type for values argument: {type(values)}"
        )


class SparseNDArrayRead(somacore.SparseRead):
    """
    Intermediate type to choose result format when reading a sparse array.

    See docs for somacore.data.SparseRead

    [lifecycle: experimental]
    """

    def __init__(self, sr: clib.SOMAReader, shape: NTuple):
        """
        [lifecycle: experimental]
        """
        self.sr = sr
        self.shape = shape

    def coos(self) -> SparseCOOTensorReadIter:
        """
        Return an iterator of Arrow SparseCOOTensor

        [lifecycle: experimental]
        """
        return SparseCOOTensorReadIter(self.sr, self.shape)

    def cscs(self) -> SparseCSCMatrixReadIter:
        """
        Return an iterator of Arrow SparseCSCMatrix

        [lifecycle: experimental]
        """

        return SparseCSCMatrixReadIter(self.sr, self.shape)

    def csrs(self) -> SparseCSRMatrixReadIter:
        """
        Return an iterator of Arrow SparseCSRMatrix

        [lifecycle: experimental]
        """

        return SparseCSRMatrixReadIter(self.sr, self.shape)

    def dense_tensors(self) -> somacore.ReadIter[pa.Tensor]:
        """
        Return an iterator of Arrow Tensor

        [lifecycle: experimental]
        """

        raise NotImplementedError()

    def tables(self) -> TableReadIter:
        """
        Return an iterator of Arrow Table

        [lifecycle: experimental]
        """
        return TableReadIter(self.sr)
