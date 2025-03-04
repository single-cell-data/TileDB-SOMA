# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""
Implementation of SOMA DenseNDArray.
"""

from __future__ import annotations

import warnings
from typing import List, Sequence, Tuple, Union

import numpy as np
import pyarrow as pa
import somacore
from somacore import options
from typing_extensions import Self

from . import _util
from . import pytiledbsoma as clib
from ._arrow_types import pyarrow_to_carrow_type
from ._common_nd_array import NDArray
from ._exception import SOMAError, map_exception_for_create
from ._read_iters import ManagedQuery, TableReadIter
from ._tdb_handles import DenseNDArrayWrapper
from ._types import OpenTimestamp, Slice
from ._util import dense_indices_to_shape
from .options._soma_tiledb_context import (
    SOMATileDBContext,
    _validate_soma_tiledb_context,
)
from .options._tiledb_create_write_options import (
    TileDBCreateOptions,
    TileDBWriteOptions,
)


class DenseNDArray(NDArray, somacore.DenseNDArray):
    """:class:`DenseNDArray` is a dense, N-dimensional array, with offset (zero-based)
    integer indexing on each dimension. :class:`DenseNDArray` has a user-defined
    schema, which includes:

    * The element type, expressed as an
      `Arrow type <https://arrow.apache.org/docs/python/api/datatypes.html>`_
      indicating the type of data contained within the array.
    * The shape of the array, i.e., the number of dimensions and the length of
      each dimension.

    All dimensions must have a positive, non-zero length, and there must be 1
    or more dimensions.

    Where explicitly referenced in the API, the dimensions are named
    ``soma_dim_N``, where N is the dimension number (e.g., ``soma_dim_0``),
    and elements are named ``soma_data``.

    Lifecycle:
        Maturing.

    Examples:
        >>> import tiledbsoma
        >>> import pyarrow as pa
        >>> import numpy as np
        >>> with tiledbsoma.DenseNDArray.create(
        ...     "./test_dense_ndarray", type=pa.int32(), shape=(2, 3, 4)
        ... ) as arr:
        ...     data = pa.Tensor.from_numpy(
        ...         np.random.default_rng().integers(0, 10, 24).reshape(2, 3, 4)
        ...     )
        ...     arr.write((slice(None),), data)
        ... with tiledbsoma.open("./test_dense_ndarray") as arr:
        ...     print(arr.schema)
        ...     print("---")
        ...     print(arr.read())
        ...
        soma_dim_0: int64
        soma_dim_1: int64
        soma_dim_2: int64
        soma_data: int32
        ---
        <pyarrow.Tensor>
        type: int32
        shape: (2, 3, 4)
        strides: (48, 16, 4)
    """

    __slots__ = ()

    _wrapper_type = DenseNDArrayWrapper

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

        index_column_schema = []
        index_column_data = {}
        ndim = len(shape)

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

            if dim_shape is None:
                raise ValueError("DenseNDArray shape slots must be numeric")

            dim_capacity, dim_extent = cls._dim_capacity_and_extent(
                dim_name,
                # The user specifies current domain -- this is the max domain
                # which is taken from the max ranges for the dim datatype.
                # We pass None here to detect those.
                None,
                ndim,
                TileDBCreateOptions.from_platform_config(platform_config),
            )

            if dim_shape == 0:
                raise ValueError("DenseNDArray shape slots must be at least 1")

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
            clib.SOMADenseNDArray.create(
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

    def read(
        self,
        coords: options.DenseNDCoords = (),
        *,
        result_order: options.ResultOrderStr = somacore.ResultOrder.ROW_MAJOR,
        partitions: options.ReadPartitions | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> pa.Tensor:
        """Reads a user-defined dense slice of the array and return as an Arrow ``Tensor``.

        Coordinates must specify a contiguous subarray, and the number of coordinates
        must be less than or equal to the number of dimensions. For example,
        if the array is 10 by 20, then some acceptable values of ``coords`` include
        ``(3, 4)``, ``(slice(5, 10),)``, and ``(slice(5, 10), slice(6, 12))``.
        Slice indices are doubly inclusive.

        Args:
            coords:
                The coordinates for slicing the array.
            result_order:
                Order of read results.
                This can be one of 'row-major' (default) or 'column-major'
            partitions:
                An optional :class:`ReadPartitions` hint to indicate
                how results should be organized.

        Raises:
            SOMAError:
                If the object is not open for reading.

        Lifecycle:
            Maturing.
        """
        del partitions  # Currently unused.
        self._check_open_read()
        result_order = somacore.ResultOrder(result_order)

        # The dense_indices_to_shape includes, as one of its roles, how to handle default
        # coordinates -- e.g. `dnda.read()`. The default for a DenseNDArray should be "all the data"
        # -- but what is that? If the schema shape matches the non-empty domain -- e.g. at create,
        # shape was 100x200, and at write, 100x200 cells were written, those are both the same. But
        # if the array was written with room for growth -- e.g. created with shape
        # 1,000,000x1,000,000 but only 100x200 cells were written -- then we need the non-empty
        # domain.
        #
        # The non-empty domain is the correct choice in either case.
        #
        # The only exception is if the array has been created but no data have been written at
        # all, in which case the best we can do is use the schema shape.
        handle: clib.SOMADenseNDArray = self._handle._handle

        if result_order == options.ResultOrder.AUTO:
            warnings.warn(
                "The use of 'result_order=\"auto\"' is deprecated and will be "
                "removed in future versions. Please use 'row-order' (the default "
                "if no option is provided) or 'col-order' instead.",
                DeprecationWarning,
            )
            result_order = somacore.ResultOrder.ROW_MAJOR

        target_shape = dense_indices_to_shape(coords, tuple(handle.shape), result_order)

        arrow_table = TableReadIter(
            array=self,
            coords=coords,
            column_names=[],
            result_order=_util.to_clib_result_order(result_order),
            value_filter=None,
            platform_config=platform_config,
        ).concat()

        if arrow_table is None:
            raise SOMAError(
                "internal error: at least one table-piece should have been returned"
            )

        npval = arrow_table.column("soma_data").to_numpy()
        # TODO: as currently coded we're looking at the non-empty domain upper
        # bound but not its lower bound. That works fine if data are written at
        # the start: e.g. domain (0, 99) and data written at 0,1,2,3,4. It
        # doesn't work fine if data are written at say (40,41,42,43,44).
        #
        # This is tracked on https://github.com/single-cell-data/TileDB-SOMA/issues/3271
        reshaped = npval.reshape(target_shape)
        return pa.Tensor.from_numpy(reshaped)

    def write(
        self,
        coords: options.DenseNDCoords,
        values: pa.Tensor,
        *,
        platform_config: options.PlatformConfig | None = None,
    ) -> Self:
        """Writes a subarray, defined by ``coords`` and ``values``. Will overwrite existing
        values in the array.

        Args:
            coords:
                A per-dimension tuple of scalars or slices
                defining the bounds of the subarray to be written.
            values:
                The values to be written to the subarray.  Must have
                the same shape as ``coords``, and the type must match the DenseNDArray.
            platform_config:
                Optional platform-specific options to use
                in this write operation (currently unused).

        Raises:
            TypeError:
                If the ``values`` parameter is an unsupported type.
            SOMAError:
                If the object is not open for writing.

        Lifecycle:
            Maturing.
        """
        _util.check_type("values", values, (pa.Tensor,))

        clib_handle = self._handle._handle

        # Compute the coordinates for the dense array.
        new_coords: List[Union[int, Slice[int], None]] = []
        for c in coords:
            if isinstance(c, slice) and isinstance(c.stop, int):
                new_coords.append(slice(c.start, c.stop, c.step))
            else:
                new_coords.append(c)

        # Convert data to a numpy array.
        dtype = self.schema.field("soma_data").type.to_pandas_dtype()
        input = np.array(values, dtype=dtype)

        # Set the result order. If neither row nor col major, set to be row major.
        if input.flags.f_contiguous:
            order = clib.ResultOrder.colmajor
        else:
            if not input.flags.contiguous:
                input = np.ascontiguousarray(input)
            order = clib.ResultOrder.rowmajor

        mq = ManagedQuery(self, platform_config)
        mq._handle.set_layout(order)
        _util._set_coords(mq, new_coords)
        mq._handle.set_column_data("soma_data", input)
        mq._handle.submit_write()

        tiledb_write_options = TileDBWriteOptions.from_platform_config(platform_config)
        if tiledb_write_options.consolidate_and_vacuum:
            clib_handle.consolidate_and_vacuum()
        return self

    @classmethod
    def _dim_capacity_and_extent(
        cls,
        dim_name: str,
        dim_shape: int | None,
        ndim: int,
        create_options: TileDBCreateOptions,
    ) -> Tuple[int, int]:
        """Given a user-specified shape (maybe ``None``) along a particular dimension,
        returns a tuple of the TileDB capacity and extent for that dimension, suitable
        for schema creation. If the user-specified shape is None, the largest possible
        int64 is returned for the capacity -- which is particularly suitable for
        maxdomain.
        """

        # Old news: for dense n-dimensional arrays, the number of bytes for each
        # core cell is the product of the following:
        #
        # * element size (e.g. 4 for float32, 8 for float64)
        # * extent for dim 0
        # * extent for dim 1
        # * ...
        # * extent for dim n
        #
        # With float64 and extent 2048, this memory requirement for n = 1,2,3,4
        # multiplies out to:
        #
        # * n=1: ~10**4
        # * n=2: ~10**8
        # * n=3: ~10**11
        # * n=4: ~10**14
        #
        # New news: before core 2.27 we didn't have current-domain support for
        # dense arrays and so we made the core domain fit the shape -- and we of
        # course cap each dim slot's extent at the slotwise domain. So if the
        # array is 3D with shape = (10,20,30) we got tile extents 10, 20, 30.
        #
        # But as of core 2.27 we do have current-domain support for dense
        # arrays, and we make the core domain huge. (Core domain is immutable;
        # core current domain is upward-resizeable up to the limit which is the
        # core domain.) This means a default extent of 2048 blows up fast at
        # higher dimensions.
        if ndim == 1:
            default_extent = 2048
        elif ndim == 2:
            default_extent = 512
        elif ndim == 3:
            default_extent = 128
        elif ndim >= 4:
            default_extent = 4

        if dim_shape is None:
            dim_capacity = 2**63 - 1
            dim_extent = min(
                dim_capacity, create_options.dim_tile(dim_name, default_extent)
            )
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
            dim_extent = min(
                dim_shape, create_options.dim_tile(dim_name, default_extent)
            )

        return (dim_capacity, dim_extent)
