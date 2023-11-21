# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Implementation of SOMA SparseNDArray.
"""
from __future__ import annotations

from typing import (
    Dict,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

import numpy as np
import pyarrow as pa
import pyarrow.compute as pacomp
import somacore
import tiledb
from somacore import options
from somacore.options import PlatformConfig
from typing_extensions import Self

from . import _util

# This package's pybind11 code
from . import pytiledbsoma as clib
from ._common_nd_array import NDArray
from ._exception import SOMAError
from ._read_iters import (
    BlockwiseScipyReadIter,
    BlockwiseTableReadIter,
    SparseCOOTensorReadIter,
    TableReadIter,
)
from ._types import NTuple
from .options._tiledb_create_options import TileDBCreateOptions

_UNBATCHED = options.BatchSize()


class SparseNDArray(NDArray, somacore.SparseNDArray):
    """:class:`SparseNDArray` is a sparse, N-dimensional array, with offset
    (zero-based) integer indexing on each dimension.
    :class:`SparseNDArray` has a user-defined schema, which includes:

    * The element type, expressed as an
      `Arrow type <https://arrow.apache.org/docs/python/api/datatypes.html>`_,
      indicating the type of data contained within the array.
    * The shape of the array, i.e., the number of dimensions and the length of
      each dimension.

    All dimensions must have a positive, non-zero length, and there must be 1
    or more dimensions. Implicitly stored elements (i.e., those not explicitly
    stored in the array) are assumed to have a value of zero.

    Where explicitly referenced in the API, the dimensions are named
    ``soma_dim_N``, where ``N`` is the dimension number (e.g., ``soma_dim_0``),
    and elements are named ``soma_data``.

    Lifecycle:
        Experimental.

    Examples:
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

        Lifecycle:
            Experimental.
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
        """Reads a user-defined slice of the :class:`SparseNDArray`.

        Args:
            coords:
                A per-dimension tuple of scalar, slice, sequence of scalar or
                `Arrow IntegerArray <https://arrow.apache.org/docs/python/generated/pyarrow.IntegerArray.html>`_
                defining the region to read.

        Returns:
            A :class:`SparseNDArrayRead` to access result iterators in various formats.

        Raises:
            SOMAError:
                If the object is not open for reading.

        Lifecycle:
            Experimental.

        Notes:
            Acceptable ways to index:

            * A sequence of coordinates is accepted, one per dimension.
            * Sequence length must be <= number of dimensions.
            * If the sequence contains missing coordinates (length < number of dimensions),
              then ``slice(None)`` -- i.e. no constraint -- is assumed for the
              remaining dimensions.
            * Per-dimension, explicitly specified coordinates can be one of:
              None, a value, a list/``numpy.ndarray``/``pyarrow.Array``/etc of values, a slice, etc.
            * Slices are doubly inclusive: ``slice(2,4)`` means [2,3,4] not [2,3].
              Slice steps can only be +1. Slices can be ``slice(None)``, meaning
              select all in that dimension, and may be half-specified, e.g.
              ``slice(2,None)`` or ``slice(None,4)``.
            * Negative indexing is unsupported.
        """
        del batch_size, platform_config  # Currently unused.
        self._check_open_read()
        _util.check_unpartitioned(partitions)

        sr = self._soma_reader(schema=self._handle.schema, result_order=result_order)
        return SparseNDArrayRead(sr, self, coords)

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
        Writes an Arrow object to the SparseNDArray.

        `Arrow SparseTensor <https://arrow.apache.org/docs/cpp/api/tensor.html>`_:
        the coordinates in the Arrow SparseTensor are interpreted as the
        coordinates to write to. Supports the _experimental_ SparseCOOTensor,
        SparseCSRMatrix and SparseCSCMatrix.  There is currently no support for
        Arrow SparseCSFTensor or dense Tensor.

        `Arrow table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_:
        write a COO table, with columns named ``soma_dim_0``, ...,
        ``soma_dim_N`` and ``soma_data`` to the dense nD array.

        Raises:
            TypeError:
                If the ``values`` parameter is an unsupported type.
            SOMAError:
                If the object is not open for writing.

        Lifecycle:
            Experimental.
        """

        arr = self._handle.writer
        tiledb_create_options = TileDBCreateOptions.from_platform_config(
            platform_config
        )

        if isinstance(values, pa.SparseCOOTensor):
            # Write bulk data
            data, coords = values.to_numpy()
            arr[tuple(c for c in coords.T)] = data

            # Write bounding-box metadata. Note COO can be N-dimensional.
            maxes = [e - 1 for e in values.shape]
            bounding_box = self._compute_bounding_box_metadata(maxes)
            self._set_bounding_box_metadata(bounding_box)

            if tiledb_create_options.consolidate_and_vacuum:
                # Consolidate non-bulk data
                self._consolidate_and_vacuum()
            return self

        if isinstance(values, (pa.SparseCSCMatrix, pa.SparseCSRMatrix)):
            if self.ndim != 2:
                raise ValueError(
                    f"Unable to write 2D Arrow sparse matrix to {self.ndim}D SparseNDArray"
                )
            # Write bulk data
            # TODO: the ``to_scipy`` function is not zero copy. Need to explore zero-copy options.
            sp = values.to_scipy().tocoo()
            arr[sp.row, sp.col] = sp.data

            # Write bounding-box metadata. Note CSR and CSC are necessarily 2-dimensional.
            nr, nc = values.shape
            bounding_box = self._compute_bounding_box_metadata([nr - 1, nc - 1])
            self._set_bounding_box_metadata(bounding_box)

            if tiledb_create_options.consolidate_and_vacuum:
                # Consolidate non-bulk data
                self._consolidate_and_vacuum()
            return self

        if isinstance(values, pa.Table):
            # Write bulk data
            data = values.column("soma_data").to_numpy()
            coord_tbl = values.drop(["soma_data"])
            coords = tuple(
                coord_tbl.column(f"soma_dim_{n}").to_numpy()
                for n in range(coord_tbl.num_columns)
            )
            arr[coords] = data

            # Write bounding-box metadata
            maxes = []
            for i in range(coord_tbl.num_columns):
                coords = values.column(f"soma_dim_{i}")
                if coords:
                    maxes.append(pacomp.max(coords).as_py())
                else:  # completely empty X
                    maxes.append(0)
            bounding_box = self._compute_bounding_box_metadata(maxes)
            self._set_bounding_box_metadata(bounding_box)

            if tiledb_create_options.consolidate_and_vacuum:
                # Consolidate non-bulk data
                self._consolidate_and_vacuum()
            return self

        raise TypeError(
            f"Unsupported Arrow type or non-arrow type for values argument: {type(values)}"
        )

    def _set_reader_coord(
        self, sr: clib.SOMAArray, dim_idx: int, dim: tiledb.Dim, coord: object
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
        """Given a user-specified shape (maybe ``None``) along a particular dimension,
        returns a tuple of the TileDB capacity and extent for that dimension, suitable
        for schema creation. If the user-specified shape is None, the largest possible
        int64 is returned for the capacity.
        """
        if dim_shape is None:
            dim_capacity = 2**31 - 2  # Make this friendly for reads by tiledbsoma-r
            dim_extent = min(dim_capacity, create_options.dim_tile(dim_name, 2048))
        else:
            if dim_shape <= 0:
                raise ValueError(
                    "SOMASparseNDArray shape must be a non-zero-length tuple of positive ints or Nones"
                )
            dim_capacity = dim_shape
            dim_extent = min(dim_shape, create_options.dim_tile(dim_name, 2048))

        return (dim_capacity, dim_extent)

    def used_shape(self) -> Tuple[Tuple[int, int], ...]:
        """
        Retrieve the range of indexes for a dimension that were explicitly written.
        Compare this to ``shape`` which returns the available/writable capacity.
        """
        retval = []
        with tiledb.open(self.uri, ctx=self.context.tiledb_ctx):
            for i in range(20):
                lower_key = f"soma_dim_{i}_domain_lower"
                lower_val = self.metadata.get(lower_key)
                upper_key = f"soma_dim_{i}_domain_upper"
                upper_val = self.metadata.get(upper_key)
                if lower_val is None or upper_val is None:
                    break
                retval.append((lower_val, upper_val))
        if not retval:
            raise SOMAError(
                f"Array {self.uri} was not written with bounding box support. "
                + "For an approximation, please use `non_empty_domain()` instead",
            )

        # In the unlikely event that a previous data update succeeded but the
        # subsequent metadata update did not, take the union of the core non-empty domain
        # (which is done as part of the data update) and the metadata bounding box.
        ned = self.non_empty_domain()
        for i, nedslot in enumerate(ned):
            ned_lower, ned_upper = nedslot
            bbox_lower, bbox_upper = retval[i]
            retval[i] = (min(ned_lower, bbox_lower), max(ned_upper, bbox_upper))
        return tuple(retval)

    def _compute_bounding_box_metadata(
        self,
        maxes: Sequence[int],
    ) -> Dict[str, int]:
        """
        This computes a bounding box for create or update. The former applies to initial ingest;
        the latter applies to append mode.
        """
        new_bounding_box = {}
        for i, slotmax in enumerate(maxes):
            lower_key = f"soma_dim_{i}_domain_lower"
            upper_key = f"soma_dim_{i}_domain_upper"
            old_lower = self.metadata.get(lower_key)
            old_upper = self.metadata.get(upper_key)

            if old_lower is None:
                new_lower = 0
            else:
                new_lower = min(0, old_lower)

            if old_upper is None:
                new_upper = slotmax
            else:
                new_upper = max(slotmax, old_upper)

            new_bounding_box[lower_key] = new_lower
            new_bounding_box[upper_key] = new_upper
        return new_bounding_box

    def _set_bounding_box_metadata(
        self,
        bounding_box: Dict[str, int],
    ) -> None:
        """Writes the bounding box to metadata storage."""
        self.metadata.update(bounding_box)


class _SparseNDArrayReadBase(somacore.SparseRead):
    """Base class for sparse reads"""

    def __init__(
        self,
        sr: clib.SOMAArray,
        array: SparseNDArray,
        coords: options.SparseNDCoords,
    ):
        """
        Lifecycle:
            Experimental.
        """
        self.sr = sr
        self.shape = tuple(sr.shape)
        self.array = array
        self.coords = coords


class SparseNDArrayRead(_SparseNDArrayReadBase):
    """Intermediate type to choose result format when reading a sparse array.

    Results returned by `coos` and `tables` iterate over COO coordinates in the user-specified result order,
    but with breaks between iterator steps at arbitrary coordinates (i.e., any given result may split a row or
    column across two separate steps of the iterator). See `blockwise` for iterators that will always yield
    complete "blocks" for any given user-specified dimension, eg., all coordinates in a given row in one
    iteration step. NB: `blockwise` iterators may utilize additional disk or network IO.

    See also:
        somacore.data.SparseRead

    Lifecycle:
        Experimental.
    """

    def coos(self, shape: Optional[NTuple] = None) -> SparseCOOTensorReadIter:
        """
        Returns an iterator of
        `Arrow SparseCOOTensor <https://arrow.apache.org/docs/cpp/api/tensor.html>`_.

        Args:
            shape:
                Optionally, a tuple that overrides the default capacity.

        Lifecycle:
            Experimental.
        """
        if shape is not None and (len(shape) != len(self.shape)):
            raise ValueError(f"shape must be a tuple of size {len(self.shape)}")
        self.array._set_reader_coords(self.sr, self.coords)
        return SparseCOOTensorReadIter(self.sr, shape or self.shape)

    def tables(self) -> TableReadIter:
        """
        Returns an iterator of
        `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_.

        Lifecycle:
            Experimental.
        """
        self.array._set_reader_coords(self.sr, self.coords)
        return TableReadIter(self.sr)

    def blockwise(
        self,
        axis: Union[int, Sequence[int]],
        *,
        size: Optional[Union[int, Sequence[int]]] = None,
        reindex_disable_on_axis: Optional[Union[int, Sequence[int]]] = None,
        eager: bool = True,
    ) -> SparseNDArrayBlockwiseRead:
        """
        Returns an intermediate type to choose a blockwise iterator of a specific format.

        Blockwise iterators yield results grouped by a user-specified axis. For example, a
        blockwise iterator with `axis=0` will yield results containing all coordinates for
        a given "row" in the array, regardless of the read `result_order` (i.e., the sort order).

        Blockwise iterators yield an array "block" in some user-specified format, as well as a
        list of coordinates contained in the individual block.

        All blockwise iterators will reindex coordinates (i.e., map them from soma_joinid to an integer
        in the range [0, N)), unless reindexing is specifically disabled for that axis, using the
        `reindex_disable_on_axis` argument. When reindexing:
        * the primary iterator axis coordinates, as indicated by the `axis` argument, will be reindexed into the range
          `[0, N)`, where `N` is the number of coordinates read for the block (controlled with the `size` argument).
        * all other axes will be reindexed to `[0, M)`, where `M` is the number of points read
          on that axis across all blocks.

        Args:
            axis:
                Required. The axis across which to yield blocks, indicated as the dimension number, e.g.,
                `axis=0` will step across `soma_dim_0` (the first dimension).
            size:
                Optional. Number of coordinates in each block yielded by the iterator. A reasonable default will
                be provided if the argument is omitted. Current defaults are 2^16 for dimension 0 and 2^8 for
                all other dimensions. Defaults are subject to change and will likely remain relatively small.
            reindex_disable_on_axis:
                Optional. Axis or sequence of axes which will _not_ be reindexed. Defaults to None, indicating
                all axes will be reindexed.
            eager:
                Optional. If `True`, the iterator will read ahead (using multi-threading) to improve overall
                performance when iterating over a large result. Setting this flag to `False` will reduce memory
                consumption, at the cost of additional processing time.

        Examples:

            A simple example iterating over the first 10000 elements of the first dimension, into
            blocks of SciPy sparse matrices:

            >>> import tiledbsoma
            >>> with tiledbsoma.open("a_sparse_nd_array") as X:
            ...     for (obs_coords, var_coords), matrix in X.read(
            ...         coords=(slice(9999),)
            ...     ).blockwise(
            ...         axis=0, size=4999
            ...     ).scipy():
            ...         print(repr(matrix))
            <4999x60664 sparse matrix of type '<class 'numpy.float32'>'
                    with 11509741 stored elements in Compressed Sparse Row format>
            <4999x60664 sparse matrix of type '<class 'numpy.float32'>'
                    with 13760197 stored elements in Compressed Sparse Row format>
            <2x60664 sparse matrix of type '<class 'numpy.float32'>'
                    with 3417 stored elements in Compressed Sparse Row format>

            To stride over the second dimension, returning a CSC matrix, specify `blockwise(axis=1)`.
            To iterate over COO matrices, on either axis, specify `scipy(compress=False)`.

        Lifecycle:
            Experimental.
        """
        return SparseNDArrayBlockwiseRead(
            self.sr,
            self.array,
            self.coords,
            axis,
            size=size,
            reindex_disable_on_axis=reindex_disable_on_axis,
            eager=eager,
        )


class SparseNDArrayBlockwiseRead(_SparseNDArrayReadBase):
    def __init__(
        self,
        sr: clib.SOMAArray,
        array: SparseNDArray,
        coords: options.SparseNDCoords,
        axis: Union[int, Sequence[int]],
        *,
        size: Optional[Union[int, Sequence[int]]],
        reindex_disable_on_axis: Optional[Union[int, Sequence[int]]],
        eager: bool = True,
    ):
        super().__init__(sr, array, coords)
        self.axis = axis
        self.size = size
        self.reindex_disable_on_axis = reindex_disable_on_axis
        self.eager = eager

    def tables(self) -> BlockwiseTableReadIter:
        """
        Returns a blockwise iterator of
        `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_.

        Yields:
            The iterator will yield a tuple of:
                Arrow Table, (dim_0_coords, ...),
            where the second element is a tuple containing the soma_joinid values for the queried array dimensions.

        Lifecycle:
            Experimental.
        """
        return BlockwiseTableReadIter(
            self.array,
            self.sr,
            self.coords,
            self.axis,
            size=self.size,
            reindex_disable_on_axis=self.reindex_disable_on_axis,
            eager=self.eager,
        )

    def coos(self) -> somacore.ReadIter[None]:
        """
        Unimplemented due to ARROW-17933, https://issues.apache.org/jira/browse/ARROW-17933,
        which causes failure on empty tensors (which are commonly yielded by blockwise
        iterators). Also tracked as https://github.com/single-cell-data/TileDB-SOMA/issues/668
        """
        raise NotImplementedError(
            "Blockwise SparseCOOTensor not implemented due to ARROW-17933."
        )

    def scipy(self, *, compress: bool = True) -> BlockwiseScipyReadIter:
        """
        Returns a blockwise iterator of
        `SciPy sparse matrix` <https://docs.scipy.org/doc/scipy/reference/sparse.html>
        over a 2D SparseNDArray.

        Args:
            compress:
                If True, a CSC or CSR matrix is returned, dependent on the value of the
                `axis` argument of the `blockwise()` method. If False, a COO matrix is returned.

                Note: implementation details of SciPy CSC and CSR compression effectively require
                reindexing of the major axis (columns and rows, respectively). Therefore, this method
                will throw an error if the major axis was included in the reindex_disable_on_axis argument to
                `blockwise()`. Reindexing can be disabled for the minor axis.

        Yields:
            The iterator will yield a tuple of:
                SciPy sparse matrix, (obs_coords, var_coords),
            where the first element is a tuple containing the soma_joinid values for the queried array dimensions.

        Lifecycle:
            Experimental.
        """
        return BlockwiseScipyReadIter(
            self.array,
            self.sr,
            self.coords,
            self.axis,
            size=self.size,
            compress=compress,
            reindex_disable_on_axis=self.reindex_disable_on_axis,
            eager=self.eager,
        )
