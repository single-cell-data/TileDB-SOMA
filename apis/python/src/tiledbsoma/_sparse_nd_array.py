# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""
Implementation of SOMA SparseNDArray.
"""

from __future__ import annotations

from typing import (
    Sequence,
    Tuple,
    Union,
    cast,
)

import numpy as np
import pyarrow as pa
import somacore
from somacore import options
from somacore.options import PlatformConfig
from typing_extensions import Self

from . import _util

# This package's pybind11 code
from . import pytiledbsoma as clib
from ._arrow_types import pyarrow_to_carrow_type
from ._common_nd_array import NDArray
from ._exception import SOMAError, map_exception_for_create
from ._read_iters import (
    BlockwiseScipyReadIter,
    BlockwiseTableReadIter,
    ManagedQuery,
    SparseCOOTensorReadIter,
    TableReadIter,
)
from ._tdb_handles import SparseNDArrayWrapper
from ._types import NTuple, OpenTimestamp
from .options._soma_tiledb_context import (
    SOMATileDBContext,
    _validate_soma_tiledb_context,
)
from .options._tiledb_create_write_options import (
    TileDBCreateOptions,
    TileDBWriteOptions,
)

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
        Maturing.

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

    _wrapper_type = SparseNDArrayWrapper

    # Inherited from somacore
    # * ndim accessor
    # * is_sparse: Final = True

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        type: pa.DataType,
        shape: Sequence[Union[int, None]],
        platform_config: options.PlatformConfig | None = None,
        context: SOMATileDBContext | None = None,
        tiledb_timestamp: OpenTimestamp | None = None,
    ) -> Self:
        context = _validate_soma_tiledb_context(context)

        # SOMA-to-core mappings:
        #
        # Before the current-domain feature was enabled (possible after core 2.25):
        #
        # * SOMA shape <-> core domain, AKA "max domain" which is a name we'll use for clarity
        #   o Specifically, (0, SOMA shape minus 1) = core domain
        # * core current domain did not exist
        #
        # After the current-domain feature was enabled:
        #
        # * SOMA maxshape <-> core domain, AKA "max domain" which is a name we'll use for clarity
        #   o Specifically, (0, SOMA maxshape minus 1) = core max domain
        # * SOMA shape <-> core current domain
        #   o Specifically, (0, SOMA shape minus 1) = core current domain

        # As far as the user is concerned, the SOMA shape is the _only_ thing they see and care
        # about. It's resizeable (up to max_domain anyway), reads and writes are bounds-checked
        # against it, etc.

        index_column_schema = []
        index_column_data = {}

        for dim_idx, dim_shape in enumerate(shape):
            dim_name = f"soma_dim_{dim_idx}"

            pa_field = pa.field(dim_name, pa.int64())
            index_column_schema.append(pa_field)

            # Here is our Arrow data API for communicating schema info between
            # Python/R and C++ libtiledbsoma:
            #
            # [0] core max domain lo
            # [1] core max domain hi
            # [2] core extent parameter
            # If present, these next two signal to use the current-domain feature:
            # [3] core current domain lo
            # [4] core current domain hi

            dim_capacity, dim_extent = cls._dim_capacity_and_extent(
                dim_name,
                # The user specifies current domain -- this is the max domain
                # which is taken from the max ranges for the dim datatype.
                # We pass None here to detect those.
                None,
                len(shape),
                TileDBCreateOptions.from_platform_config(platform_config),
            )

            if dim_shape == 0:
                raise ValueError("SparseNDArray shape slots must be at least 1")
            if dim_shape is None:
                # Core current-domain semantics are (lo, hi) with both
                # inclusive, with lo <= hi. This means smallest is (0, 0)
                # which is shape 1, not 0.
                dim_shape = 1

            index_column_data[pa_field.name] = [
                0,
                dim_capacity - 1,
                dim_extent,
                0,
                dim_shape - 1,
            ]

        index_column_info = pa.RecordBatch.from_pydict(
            index_column_data, schema=pa.schema(index_column_schema)
        )

        carrow_type = pyarrow_to_carrow_type(type)
        plt_cfg = _util.build_clib_platform_config(platform_config)
        timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
        try:
            clib.SOMASparseNDArray.create(
                uri,
                format=carrow_type,
                index_column_info=index_column_info,
                ctx=context.native_context,
                platform_config=plt_cfg,
                timestamp=(0, timestamp_ms),
            )
        except SOMAError as e:
            raise map_exception_for_create(e, uri) from None

        handle = cls._wrapper_type.open(uri, "w", context, tiledb_timestamp)
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    @property
    def nnz(self) -> int:
        """
        The number of stored values in the array, including explicitly stored zeros.

        Lifecycle:
            Maturing.
        """
        self._check_open_read()
        return cast(SparseNDArrayWrapper, self._handle).nnz

    def read(
        self,
        coords: options.SparseNDCoords = (),
        *,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        platform_config: PlatformConfig | None = None,
    ) -> "SparseNDArrayRead":
        """Reads a user-defined slice of the :class:`SparseNDArray`.

        Args:
            coords:
                A per-dimension ``Sequence`` of scalar, slice, sequence of scalar or
                `Arrow IntegerArray <https://arrow.apache.org/docs/python/generated/pyarrow.IntegerArray.html>` values
                defining the region to read.

        Returns:
            A :class:`SparseNDArrayRead` to access result iterators in various formats.

        Raises:
            SOMAError:
                If the object is not open for reading.

        Lifecycle:
            Maturing.

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
        del batch_size  # Currently unused.
        self._check_open_read()
        _util.check_unpartitioned(partitions)

        return SparseNDArrayRead(
            array=self,
            coords=coords,
            result_order=_util.to_clib_result_order(result_order),
            platform_config=platform_config,
        )

    def write(
        self,
        values: Union[
            pa.SparseCOOTensor,
            pa.SparseCSRMatrix,
            pa.SparseCSCMatrix,
            pa.Table,
        ],
        *,
        platform_config: PlatformConfig | None = None,
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
            Maturing.
        """

        write_options: Union[TileDBCreateOptions, TileDBWriteOptions]
        sort_coords = None
        if isinstance(platform_config, TileDBCreateOptions):
            raise ValueError(
                "As of TileDB-SOMA 1.13, the write method takes "
                "TileDBWriteOptions instead of TileDBCreateOptions"
            )
        write_options = TileDBWriteOptions.from_platform_config(platform_config)
        sort_coords = write_options.sort_coords

        clib_sparse_array = self._handle._handle

        if isinstance(values, pa.SparseCOOTensor):
            # Write bulk data
            data, coords = values.to_numpy()

            mq = ManagedQuery(self, platform_config)
            for i, c in enumerate(coords.T):
                mq._handle.set_column_data(
                    f"soma_dim_{i}",
                    np.array(
                        c,
                        dtype=self.schema.field(f"soma_dim_{i}").type.to_pandas_dtype(),
                    ),
                )
            mq._handle.set_column_data(
                "soma_data",
                np.array(
                    data, dtype=self.schema.field("soma_data").type.to_pandas_dtype()
                ),
            )
            mq._handle.submit_write(sort_coords or True)

            if write_options.consolidate_and_vacuum:
                # Consolidate non-bulk data
                clib_sparse_array.consolidate_and_vacuum()
            return self

        if isinstance(values, (pa.SparseCSCMatrix, pa.SparseCSRMatrix)):
            if self.ndim != 2:
                raise ValueError(
                    f"Unable to write 2D Arrow sparse matrix to {self.ndim}D SparseNDArray"
                )
            # Write bulk data
            # TODO: the ``to_scipy`` function is not zero copy. Need to explore zero-copy options.
            sp = values.to_scipy().tocoo()

            mq = ManagedQuery(self, platform_config)
            for i, c in enumerate([sp.row, sp.col]):
                mq._handle.set_column_data(
                    f"soma_dim_{i}",
                    np.array(
                        c,
                        dtype=self.schema.field(f"soma_dim_{i}").type.to_pandas_dtype(),
                    ),
                )
            mq._handle.set_column_data(
                "soma_data",
                np.array(
                    sp.data, dtype=self.schema.field("soma_data").type.to_pandas_dtype()
                ),
            )
            mq._handle.submit_write(sort_coords or True)

            if write_options.consolidate_and_vacuum:
                # Consolidate non-bulk data
                clib_sparse_array.consolidate_and_vacuum()
            return self

        if isinstance(values, pa.Table):
            # Write bulk data
            for batch in values.to_batches():
                # clib_sparse_array.write(batch, sort_coords or False)
                mq = ManagedQuery(self, None)
                mq._handle.set_array_data(batch)
                mq._handle.submit_write(sort_coords or False)

            if write_options.consolidate_and_vacuum:
                # Consolidate non-bulk data
                clib_sparse_array.consolidate_and_vacuum()
            return self

        raise TypeError(
            f"Unsupported Arrow type or non-arrow type for values argument: {type(values)}"
        )

    @classmethod
    def _dim_capacity_and_extent(
        cls,
        dim_name: str,
        dim_shape: int | None,
        ndim: int,  # not needed for sparse
        create_options: TileDBCreateOptions,
    ) -> Tuple[int, int]:
        """Given a user-specified shape (maybe ``None``) along a particular dimension,
        returns a tuple of the TileDB capacity and extent for that dimension, suitable
        for schema creation. If the user-specified shape is None, the largest possible
        int64 is returned for the capacity -- which is particularly suitable for
        maxdomain.
        """
        if dim_shape is None:
            dim_capacity = 2**63 - 1
            dim_extent = min(dim_capacity, create_options.dim_tile(dim_name, 2048))
            # For core: "domain max expanded to multiple of tile extent exceeds max value
            # representable by domain type. Reduce domain max by 1 tile extent to allow for
            # expansion."
            dim_capacity -= dim_extent
        else:
            if dim_shape <= 0:
                raise ValueError(
                    "SOMASparseNDArray shape must be a non-zero-length tuple of positive ints or Nones"
                )
            dim_capacity = dim_shape
            dim_extent = min(dim_shape, create_options.dim_tile(dim_name, 2048))

        return (dim_capacity, dim_extent)


class _SparseNDArrayReadBase(somacore.SparseRead):
    """Base class for sparse reads"""

    def __init__(
        self,
        array: SparseNDArray,
        coords: options.SparseNDCoords,
        result_order: clib.ResultOrder,
        platform_config: options.PlatformConfig | None,
    ):
        """
        Lifecycle:
            Maturing.
        """
        self.array = array
        self.coords = coords
        self.shape = tuple(array._handle._handle.shape)
        self.result_order = result_order
        self.platform_config = platform_config


class SparseNDArrayRead(_SparseNDArrayReadBase):
    """:class:`SparseNDArrayRead` is an intermediate type which supports multiple eventual result formats
     when reading a sparse array.

    Results returned by `coos` and `tables` iterate over COO coordinates in the user-specified result order,
    but with breaks between iterator steps at arbitrary coordinates (i.e., any given result may split a row or
    column across two separate steps of the iterator). See `blockwise` for iterators that will always yield
    complete "blocks" for any given user-specified dimension, eg., all coordinates in a given row in one
    iteration step. NB: `blockwise` iterators may utilize additional disk or network IO.

    See also:
        somacore.data.SparseRead

    Lifecycle:
        Maturing.
    """

    def coos(self, shape: NTuple | None = None) -> SparseCOOTensorReadIter:
        """
        Returns an iterator of
        `Arrow SparseCOOTensor <https://arrow.apache.org/docs/cpp/api/tensor.html>`_.

        Args:
            shape:
                Optionally, a tuple that overrides the default capacity.

        Lifecycle:
            Maturing.
        """
        if shape is not None and (len(shape) != len(self.shape)):
            raise ValueError(f"shape must be a tuple of size {len(self.shape)}")
        return SparseCOOTensorReadIter(
            self.array,
            self.coords,
            shape or self.shape,
            self.result_order,
            self.platform_config,
        )

    def tables(self) -> TableReadIter:
        """
        Returns an iterator of
        `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_.

        Lifecycle:
            Maturing.
        """
        return TableReadIter(
            array=self.array,
            coords=self.coords,
            column_names=[],
            result_order=self.result_order,
            value_filter=None,
            platform_config=self.platform_config,
        )

    def blockwise(
        self,
        axis: Union[int, Sequence[int]],
        *,
        size: int | Sequence[int] | None = None,
        reindex_disable_on_axis: int | Sequence[int] | None = None,
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
            Maturing.
        """
        return SparseNDArrayBlockwiseRead(
            self.array,
            self.coords,
            axis,
            self.result_order,
            self.platform_config,
            size=size,
            reindex_disable_on_axis=reindex_disable_on_axis,
            eager=eager,
        )


class SparseNDArrayBlockwiseRead(_SparseNDArrayReadBase):
    def __init__(
        self,
        array: SparseNDArray,
        coords: options.SparseNDCoords,
        axis: Union[int, Sequence[int]],
        result_order: clib.ResultOrder,
        platform_config: options.PlatformConfig | None,
        *,
        size: int | Sequence[int] | None,
        reindex_disable_on_axis: int | Sequence[int] | None,
        eager: bool = True,
    ):
        super().__init__(array, coords, result_order, platform_config)
        self.result_order = result_order
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
            Maturing.
        """
        return BlockwiseTableReadIter(
            self.array,
            self.coords,
            self.axis,
            self.result_order,
            self.platform_config,
            size=self.size,
            reindex_disable_on_axis=self.reindex_disable_on_axis,
            eager=self.eager,
            context=self.array.context,
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
            Maturing.
        """
        return BlockwiseScipyReadIter(
            self.array,
            self.coords,
            self.axis,
            self.result_order,
            self.platform_config,
            size=self.size,
            compress=compress,
            reindex_disable_on_axis=self.reindex_disable_on_axis,
            eager=self.eager,
            context=self.array.context,
        )
