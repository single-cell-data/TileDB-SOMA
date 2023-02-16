from typing import Optional, Sequence, Union, cast

import numpy as np
import pyarrow as pa
import somacore
import tiledb
from somacore import options
from somacore.options import PlatformConfig
from typing_extensions import Self

from . import _util

# This package's pybind11 code
from . import libtiledbsoma as clib
from ._common_nd_array import NDArray
from ._read_iters import (
    SparseCOOTensorReadIter,
    SparseCSCMatrixReadIter,
    SparseCSRMatrixReadIter,
    TableReadIter,
)
from ._types import NTuple

_UNBATCHED = options.BatchSize()


class SparseNDArray(NDArray, somacore.SparseNDArray):
    """
    Represents ``X`` and others.

    [lifecycle: experimental]
    """

    __slots__ = ()

    # Inherited from somacore
    # * ndim accessor
    # * is_sparse: Final = True

    @property
    def nnz(self) -> int:
        """
        The number of stored values in the array, including explicitly stored zeros.
        """
        self._check_open_read()
        return cast(int, self._soma_reader().nnz())

    def read(
        self,
        coords: options.SparseNDCoords = (),
        *,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[PlatformConfig] = None,
    ) -> "SparseNDArrayRead":
        """
        Read a user-defined slice of the SparseNDArray.

        [lifecycle: experimental]

        :param coords: A per-dimension tuple of scalar, slice, sequence of scalar
            or Arrow IntegerArray defining the region to read.
            (Arrow arrays currently unimplemented.)

        Acceptable ways to index
        ------------------------
        * A sequence of coordinates is accepted, one per dimension.
        * Sequence length must be <= number of dimensions.
        * If the sequence contains missing coordinates (length < number of dimensions),
          then ``slice(None)`` -- i.e. no constraint -- is assumed for the
          remaining dimensions.
        * Per-dimension, explicitly specified coordinates can be one of:
          None, a value, a list/ndarray/paarray/etc of values, a slice, etc.
        * Slices are doubly inclusive: slice(2,4) means [2,3,4] not [2,3].
          Slice steps can only be +1. Slices can be `slice(None)`, meaning
          select all in that dimension, but may not be half-specified:
          `slice(2,None)` and `slice(None,4)` are both unsupported.
        * Negative indexing is unsupported.

        :return: A SparseNDArrayRead to access result iterators in various formats.
        """
        del result_order, batch_size, platform_config  # Currently unused.
        self._check_open_read()
        _util.check_unpartitioned(partitions)

        schema = self._handle.schema
        sr = self._soma_reader(schema=schema)
        self._set_reader_coords(sr, coords)
        sr.submit()
        return SparseNDArrayRead(sr, schema.shape)

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
    ) -> Self:
        """
        Write an Arrow object to the SparseNDArray.

        [lifecycle: experimental]

        Arrow sparse tensor: the coordinates in the Arrow SparseTensor are interpreted
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
            return self

        if isinstance(values, (pa.SparseCSCMatrix, pa.SparseCSRMatrix)):
            if self.ndim != 2:
                raise ValueError(
                    f"Unable to write 2D Arrow sparse matrix to {self.ndim}D SparseNDArray"
                )
            # TODO: the `to_scipy` function is not zero copy. Need to explore zero-copy options.
            sp = values.to_scipy().tocoo()
            arr[sp.row, sp.col] = sp.data
            return self

        if isinstance(values, pa.Table):
            data = values.column("soma_data").to_numpy()
            coord_tbl = values.drop(["soma_data"])
            coords = tuple(
                coord_tbl.column(f"soma_dim_{n}").to_numpy()
                for n in range(coord_tbl.num_columns)
            )
            arr[coords] = data
            return self

        raise TypeError(
            f"Unsupported Arrow type or non-arrow type for values argument: {type(values)}"
        )

    def _set_reader_coord(
        self, sr: clib.SOMAReader, dim_idx: int, dim: tiledb.Dim, coord: object
    ) -> bool:
        if super()._set_reader_coord(sr, dim_idx, dim, coord):
            return True
        if isinstance(coord, (Sequence, pa.Array, pa.ChunkedArray, np.ndarray)):
            if isinstance(coord, np.ndarray) and coord.ndim != 1:
                raise ValueError(
                    f"only 1D numpy arrays may be used to index; got {coord.ndim}"
                )
            sr.set_dim_points(dim.name, coord)
            return True
        return False


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
