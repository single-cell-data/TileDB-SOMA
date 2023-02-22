from typing import Optional, Sequence, Tuple, Union, cast

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
    TableReadIter,
)
from ._types import NTuple
from .options._tiledb_create_options import TileDBCreateOptions

_UNBATCHED = options.BatchSize()


class SparseNDArray(NDArray, somacore.SparseNDArray):
    """
    ``SparseNDArray`` is a sparse, N-dimensional array, with offset (zero-based)
    integer indexing on each dimension. ``SparseNDArray`` has a user-defined
    schema, which includes:
    - the element type, expressed as an Arrow type, indicating the type of data
      contained within the array, and
    - the shape of the array, i.e., the number of dimensions and the length of
      each dimension

    All dimensions must have a positive, non-zero length, and there must be 1
    or more dimensions. Implicitly stored elements (i.e., those not explicitly
    stored in the array) are assumed to have a value of zero.

    Where explicitly referenced in the API, the dimensions are named
    ``soma_dim_N``, where ``N`` is the dimension number (e.g., ``soma_dim_0``),
    and elements are named ``soma_data``.

    [lifecycle: experimental]

    Examples:
    ---------
    >>> import tiledbsoma
    >>> import pyarrow as pa
    >>> import numpy as np
    >>> import scipy.sparse
    >>> with tiledbsoma.SparseNDArray.create(
    ...     "./test_sparse_ndarray", type=pa.float32(), shape=(1000, 100)
    ... ) as arr:
    ...     data = pa.SparseCOOTensor.from_scipy(
    ...         scipy.sparse.random(1000, 100, format="coo", dtype=np.float32)
    ...     )
    ...     arr.write(data)
    ... with tiledbsoma.SparseNDArray.open("./test_sparse_ndarray") as arr:
    ...     print(arr.schema)
    ...     print('---')
    ...     print(arr.read().coos().concat())
    ...
    soma_dim_0: int64
    soma_dim_1: int64
    soma_data: float
    ---
    <pyarrow.SparseCOOTensor>
    type: float
    shape: (1000, 100)
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
        if isinstance(coord, Sequence):
            if dim.dtype == np.int64:
                sr.set_dim_points_int64(dim.name, coord)
                return True
            elif dim.dtype == "str" or dim.dtype == "bytes":
                sr.set_dim_points_string_or_bytes(dim.name, coord)
                return True
            else:
                return False

        if isinstance(coord, np.ndarray):
            if isinstance(coord, np.ndarray) and coord.ndim != 1:
                raise ValueError(
                    f"only 1D numpy arrays may be used to index; got {coord.ndim}"
                )
            if dim.dtype == np.int64:
                sr.set_dim_points_int64(dim.name, coord)
                return True
            elif dim.dtype == "str" or dim.dtype == "bytes":
                sr.set_dim_points_string_or_bytes(dim.name, coord)
                return True

            return False
        if isinstance(coord, (pa.Array, pa.ChunkedArray)):
            sr.set_dim_points_arrow(dim.name, coord)
            return True
        return False

    @classmethod
    def _dim_capacity_and_extent(
        cls,
        dim_name: str,
        dim_shape: Optional[int],
        create_options: TileDBCreateOptions,
    ) -> Tuple[int, int]:
        """
        Given a user-specified shape (maybe ``None``) along a particular dimension,
        returns a tuple of the TileDB capacity and extent for that dimension, suitable
        for schema creation. If the user-specified shape is None, the largest possible
        int64 is returned for the capacity.
        """
        if dim_shape is None:
            dim_capacity = 2**63
            dim_extent = min(dim_capacity, create_options.dim_tile(dim_name, 2048))
            # TileDB requires that each signed-64-bit-int domain slot, rounded up to
            # a multiple of the tile extent in that slot, be representable as a
            # signed 64-bit int. So if the tile extent is 999, say, that would
            # exceed 2**63 - 1.
            dim_capacity -= dim_extent
        else:
            if dim_shape <= 0:
                raise ValueError(
                    "SOMASparseNDArray shape must be a non-zero-length tuple of positive ints or Nones"
                )
            dim_capacity = dim_shape
            dim_extent = min(dim_shape, create_options.dim_tile(dim_name, 2048))

        return (dim_capacity, dim_extent)


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
