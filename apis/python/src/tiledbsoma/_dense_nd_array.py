# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Implementation of SOMA DenseNDArray.
"""

from typing import List, Optional, Sequence, Tuple, Union

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
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
    ) -> Self:
        context = _validate_soma_tiledb_context(context)

        index_column_schema = []
        index_column_data = {}
        for dim_idx, dim_shape in enumerate(shape):
            dim_name = f"soma_dim_{dim_idx}"
            pa_field = pa.field(dim_name, pa.int64())
            dim_capacity, dim_extent = cls._dim_capacity_and_extent(
                dim_name,
                dim_shape,
                TileDBCreateOptions.from_platform_config(platform_config),
            )
            index_column_schema.append(pa_field)
            # TODO: support current domain for dense arrays once we have core support.
            # https://github.com/single-cell-data/TileDB-SOMA/issues/2955
            index_column_data[pa_field.name] = [0, dim_capacity - 1, dim_extent]

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
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[options.PlatformConfig] = None,
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
                This can be one of 'row-major', 'col-major', or 'auto'.
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

        data_shape = handle.shape
        ned = self.non_empty_domain()
        if ned is not None:
            data_shape = tuple(slot[1] + 1 for slot in ned)
        target_shape = dense_indices_to_shape(coords, data_shape, result_order)

        context = handle.context()
        if platform_config is not None:
            config = context.tiledb_config.copy()
            config.update(platform_config)
            context = clib.SOMAContext(config)

        sr = clib.SOMADenseNDArray.open(
            uri=handle.uri,
            mode=clib.OpenMode.read,
            context=context,
            column_names=[],
            result_order=_util.to_clib_result_order(result_order),
            timestamp=handle.timestamp and (0, handle.timestamp),
        )

        self._set_reader_coords(sr, coords)

        arrow_tables = []
        while True:
            arrow_table_piece = sr.read_next()
            if not arrow_table_piece:
                break
            arrow_tables.append(arrow_table_piece)

        # For dense arrays there is no zero-output case: attempting to make a test case
        # to do that, say by indexing a 10x20 array by positions 888 and 999, results
        # in read-time errors of the form
        #
        # [TileDB::Subarray] Error: Cannot add range to dimension 'soma_dim_0'; Range [888, 888] is
        # out of domain bounds [0, 9]
        if not arrow_tables:
            raise SOMAError(
                "internal error: at least one table-piece should have been returned"
            )

        arrow_table = pa.concat_tables(arrow_tables)
        return pa.Tensor.from_numpy(
            arrow_table.column("soma_data").to_numpy().reshape(target_shape)
        )

    def write(
        self,
        coords: options.DenseNDCoords,
        values: pa.Tensor,
        *,
        platform_config: Optional[options.PlatformConfig] = None,
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

        clib_dense_array = self._handle._handle

        # Compute the coordinates for the dense array.
        new_coords: List[Union[int, Slice[int], None]] = []
        for c in coords:
            if isinstance(c, slice) and isinstance(c.stop, int):
                new_coords.append(slice(c.start, c.stop - 1, c.step))
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
        clib_dense_array.reset(result_order=order)
        self._set_reader_coords(clib_dense_array, new_coords)
        clib_dense_array.write(input)

        tiledb_write_options = TileDBWriteOptions.from_platform_config(platform_config)
        if tiledb_write_options.consolidate_and_vacuum:
            clib_dense_array.consolidate_and_vacuum()
        return self

    def resize(self, newshape: Sequence[Union[int, None]]) -> None:
        """Supported for ``SparseNDArray``; scheduled for implementation for
        ``DenseNDArray`` in TileDB-SOMA 1.15
        """
        # TODO: support current domain for dense arrays once we have core support.
        # https://github.com/single-cell-data/TileDB-SOMA/issues/2955
        raise NotImplementedError()

    @classmethod
    def _dim_capacity_and_extent(
        cls,
        dim_name: str,
        dim_shape: Optional[int],
        create_options: TileDBCreateOptions,
    ) -> Tuple[int, int]:
        """Given a user-specified shape along a particular dimension, returns a tuple of
        the TileDB capacity and extent for that dimension, suitable for schema creation.
        The user-specified shape cannot be ``None`` for :class:`DenseNDArray`.
        """
        if dim_shape is None or dim_shape <= 0:
            raise ValueError(
                "SOMADenseNDArray shape must be a non-zero-length tuple of positive ints"
            )

        dim_capacity = dim_shape
        dim_extent = min(dim_shape, create_options.dim_tile(dim_name, 2048))

        return (dim_capacity, dim_extent)
