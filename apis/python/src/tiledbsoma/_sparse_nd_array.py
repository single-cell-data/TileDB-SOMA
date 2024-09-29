# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Implementation of SOMA SparseNDArray.
"""
from __future__ import annotations

import itertools
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
from somacore import options
from somacore.options import PlatformConfig
from typing_extensions import Self

from tiledbsoma._flags import NEW_SHAPE_FEATURE_FLAG_ENABLED

from . import _util

# This package's pybind11 code
from . import pytiledbsoma as clib
from ._arrow_types import pyarrow_to_carrow_type
from ._common_nd_array import NDArray
from ._exception import SOMAError, map_exception_for_create
from ._read_iters import (
    BlockwiseScipyReadIter,
    BlockwiseTableReadIter,
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
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
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

            if NEW_SHAPE_FEATURE_FLAG_ENABLED:
                dim_capacity, dim_extent = cls._dim_capacity_and_extent(
                    dim_name,
                    # The user specifies current domain -- this is the max domain
                    # which is taken from the max ranges for the dim datatype.
                    # We pass None here to detect those.
                    None,
                    TileDBCreateOptions.from_platform_config(platform_config),
                )

                if dim_shape == 0:
                    raise ValueError("SparseNDArray shape slots must be at least 1")
                if dim_shape is None:
                    dim_shape = dim_capacity

                index_column_data[pa_field.name] = [
                    0,
                    dim_capacity - 1,
                    dim_extent,
                    0,
                    dim_shape - 1,
                ]

            else:
                dim_capacity, dim_extent = cls._dim_capacity_and_extent(
                    dim_name,
                    dim_shape,
                    TileDBCreateOptions.from_platform_config(platform_config),
                )
                index_column_data[pa_field.name] = [0, dim_capacity - 1, dim_extent]

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
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[PlatformConfig] = None,
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
        handle: clib.SOMASparseNDArray = self._handle._handle

        self._check_open_read()
        _util.check_unpartitioned(partitions)

        context = handle.context()
        if platform_config is not None:
            config = context.tiledb_config.copy()
            config.update(platform_config)
            context = clib.SOMAContext(config)

        sr = clib.SOMASparseNDArray.open(
            uri=handle.uri,
            mode=clib.OpenMode.read,
            context=context,
            column_names=[],
            result_order=_util.to_clib_result_order(result_order),
            timestamp=handle.timestamp and (0, handle.timestamp),
        )

        return SparseNDArrayRead(sr, self, coords)

    def resize(self, newshape: Sequence[Union[int, None]]) -> None:
        """Increases the shape of the array as specfied. Raises an error if the new
        shape is less than the current shape in any dimension. Raises an error if
        the new shape exceeds maxshape in any dimension. Raises an error if the
        array doesn't already have a shape: in that case please call
        tiledbsoma_upgrade_shape.
        """
        self._handle.resize(newshape)

    def tiledbsoma_upgrade_shape(self, newshape: Sequence[Union[int, None]]) -> None:
        """Allows the array to have a resizeable shape as described in the TileDB-SOMA
        1.15 release notes.  Raises an error if the new shape exceeds maxshape in
        any dimension. Raises an error if the array already has a shape.
        """
        self._handle.tiledbsoma_upgrade_shape(newshape)

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
            clib_sparse_array.write_coords(
                [
                    np.array(
                        c,
                        dtype=self.schema.field(f"soma_dim_{i}").type.to_pandas_dtype(),
                    )
                    for i, c in enumerate(coords.T)
                ],
                np.array(
                    data, dtype=self.schema.field("soma_data").type.to_pandas_dtype()
                ),
                sort_coords or True,
            )

            # Write bounding-box metadata. Note COO can be N-dimensional.
            maxes = [e - 1 for e in values.shape]
            bounding_box = self._compute_bounding_box_metadata(maxes)
            self._set_bounding_box_metadata(bounding_box)

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
            clib_sparse_array.write_coords(
                [
                    np.array(
                        c,
                        dtype=self.schema.field(f"soma_dim_{i}").type.to_pandas_dtype(),
                    )
                    for i, c in enumerate([sp.row, sp.col])
                ],
                np.array(
                    sp.data, dtype=self.schema.field("soma_data").type.to_pandas_dtype()
                ),
                sort_coords or True,
            )

            # Write bounding-box metadata. Note CSR and CSC are necessarily 2-dimensional.
            nr, nc = values.shape
            bounding_box = self._compute_bounding_box_metadata([nr - 1, nc - 1])
            self._set_bounding_box_metadata(bounding_box)

            if write_options.consolidate_and_vacuum:
                # Consolidate non-bulk data
                clib_sparse_array.consolidate_and_vacuum()
            return self

        if isinstance(values, pa.Table):
            # Write bulk data
            values = _util.cast_values_to_target_schema(values, self.schema)
            for batch in values.to_batches():
                clib_sparse_array.write(batch, sort_coords or False)

            # Write bounding-box metadata
            maxes = []
            coord_tbl = values.drop(["soma_data"])
            for i in range(coord_tbl.num_columns):
                coords = values.column(f"soma_dim_{i}")
                if coords:
                    maxes.append(pacomp.max(coords).as_py())
                else:  # completely empty X
                    maxes.append(0)
            bounding_box = self._compute_bounding_box_metadata(maxes)
            self._set_bounding_box_metadata(bounding_box)

            if write_options.consolidate_and_vacuum:
                # Consolidate non-bulk data
                clib_sparse_array.consolidate_and_vacuum()
            return self

        raise TypeError(
            f"Unsupported Arrow type or non-arrow type for values argument: {type(values)}"
        )

    def _set_reader_coord(
        self, sr: clib.SOMAArray, dim_idx: int, dim: pa.Field, coord: object
    ) -> bool:
        if super()._set_reader_coord(sr, dim_idx, dim, coord):
            return True
        if isinstance(coord, Sequence):
            if pa.types.is_int64(dim.type):
                sr.set_dim_points_int64(dim.name, coord)
                return True
            elif _util.pa_types_is_string_or_bytes(dim.type):
                sr.set_dim_points_string_or_bytes(dim.name, coord)
                return True
            else:
                return False

        if isinstance(coord, np.ndarray):
            if isinstance(coord, np.ndarray) and coord.ndim != 1:
                raise ValueError(
                    f"only 1D numpy arrays may be used to index; got {coord.ndim}"
                )
            if pa.types.is_int64(dim.type):
                sr.set_dim_points_int64(dim.name, coord)
                return True
            elif _util.pa_types_is_string_or_bytes(dim.type):
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

    def used_shape(self) -> Tuple[Tuple[int, int], ...]:
        """
        Retrieve the range of indexes for a dimension that were explicitly written.
        Compare this to ``shape`` which returns the available/writable capacity.

        This method is deprecated as of TileDB-SOMA 1.13, and will be removed in TileDB-SOMA 1.15.
        """
        retval = []
        for i in itertools.count():
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
                "For an approximation, please use `non_empty_domain()` instead"
            )

        # In the unlikely event that a previous data update succeeded but the
        # subsequent metadata update did not, take the union of the core non-empty domain
        # (which is done as part of the data update) and the metadata bounding box.
        ned = self.non_empty_domain() or ()
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
            Maturing.
        """
        self.sr = sr
        self.shape = tuple(sr.shape)
        self.array = array
        self.coords = coords


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

    def coos(self, shape: Optional[NTuple] = None) -> SparseCOOTensorReadIter:
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
        self.array._set_reader_coords(self.sr, self.coords)
        return SparseCOOTensorReadIter(self.sr, shape or self.shape)

    def tables(self) -> TableReadIter:
        """
        Returns an iterator of
        `Arrow Table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_.

        Lifecycle:
            Maturing.
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
            Maturing.
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
            Maturing.
        """
        return BlockwiseTableReadIter(
            self.array,
            self.sr,
            self.coords,
            self.axis,
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
            self.sr,
            self.coords,
            self.axis,
            size=self.size,
            compress=compress,
            reindex_disable_on_axis=self.reindex_disable_on_axis,
            eager=self.eager,
            context=self.array.context,
        )
