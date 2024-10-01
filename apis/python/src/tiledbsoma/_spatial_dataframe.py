# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Implementation of a base class shared between GeometryDataFrame and PointCloudDataFrame
"""
from typing import Any, Optional, Sequence, Tuple, Type, Union

import numpy as np
import pyarrow as pa
import somacore
from somacore import CoordinateSpace, CoordinateTransform, options
from typing_extensions import Self

from . import _util
from . import pytiledbsoma as clib
from ._read_iters import TableReadIter
from ._soma_array import SOMAArray
from ._types import Slice, is_slice_of

_UNBATCHED = options.BatchSize()


class SpatialDataFrame(SOMAArray):

    __slots__ = ()

    def keys(self) -> Tuple[str, ...]:
        """Returns the names of the columns when read back as a spatial dataframe.

        Examples:
            >>> with tiledbsoma.open("a_dataframe") as soma_df:
            ...     k = soma_df.keys()
            ...
            >>> k
            ('soma_joinid', 'col1')

        Lifecycle:
            Experimental.
        """
        return self._tiledb_array_keys()

    @property
    def index_column_names(self) -> Tuple[str, ...]:
        """Returns index (dimension) column names.

        Lifecycle:
            Experimental.
        """
        return self._tiledb_dim_names()

    @property
    def axis_names(self) -> Tuple[str, ...]:
        """The names of the axes of the coordinate space the data is defined on.

        Lifecycle: experimental
        """
        raise NotImplementedError("Must be implemented by the child class")

    @property
    def domain(self) -> Tuple[Tuple[Any, Any], ...]:
        """Returns a tuple of minimum and maximum values, inclusive, storable
        on each index column of the dataframe.

        Lifecycle:
            Experimental.
        """
        return self._domain()

    def read(
        self,
        coords: options.SparseDFCoords = (),
        column_names: Optional[Sequence[str]] = None,
        *,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        value_filter: Optional[str] = None,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: Optional[options.ReadPartitions] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> TableReadIter:
        """Reads a user-defined slice of data into Arrow tables.

        Args:
            coords: for each index dimension, which rows to read.
                Defaults to ``()``, meaning no constraint -- all IDs.
            column_names: the named columns to read and return.
                Defaults to ``None``, meaning no constraint -- all column names.
            partitions: If present, specifies that this is part of
                a partitioned read, and which part of the data to include.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            value_filter: an optional value filter to apply to the results.
                The default of ``None`` represents no filter. Value filter
                syntax is implementation-defined; see the documentation
                for the particular SOMA implementation for details.
        Returns:
            A :class:`ReadIter` of :class:`pa.Table`s.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError("must be implemented by child class")

    def read_spatial_region(
        self,
        region: Optional[options.SpatialRegion] = None,
        column_names: Optional[Sequence[str]] = None,
        *,
        region_transform: Optional[CoordinateTransform] = None,
        region_coord_space: Optional[CoordinateSpace] = None,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: Optional[options.ReadPartitions] = None,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        value_filter: Optional[str] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> somacore.SpatialRead[somacore.ReadIter[pa.Table]]:
        """Reads data intersecting an user-defined region of space into a
        :class:`SpatialRead` with data in Arrow tables.


        Args:
            region: The region to query. May be a box in the form
                [x_min, y_min, x_max, y_max] (for 2D images), a box in the form
                [x_min, y_min, z_min, x_max, y_max, z_max] (for 3D images), or
                a shapely Geometry.
            column_names: The named columns to read and return.
                Defaults to ``None``, meaning no constraint -- all column names.
            region_transform: An optional coordinate transform from the read region to
                the coordinate system of the spatial dataframe.
                Defaults to ``None``, meaning an identity transform.
            region_coord_space: An optional coordinate space for the region being read.
                Defaults to ``None``, coordinate space will be inferred from transform.
            batch_size: The size of batched reads.
                Defaults to `unbatched`.
            partitions: If present, specifies that this is part of a partitioned read,
                and which part of the data to include.
            result_order: the order to return results, specified as a
                :class:`~options.ResultOrder` or its string value.
            value_filter: an optional value filter to apply to the results.
                The default of ``None`` represents no filter. Value filter
                syntax is implementation-defined; see the documentation
                for the particular SOMA implementation for details.

        Returns:
            A :class:`SpatialRead` with :class:`ReadIter` of :class:`pa.Table`s data.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError("must be implemented by child class")

    def write(
        self,
        values: Union[pa.RecordBatch, pa.Table],
        *,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> Self:
        """Writes the data from an Arrow table to the persistent object.

        As duplicate index values are not allowed, index values already present
        in the object are overwritten and new index values are added.

        Args:
            values: An Arrow table containing all columns, including
                the index columns. The schema for the values must match
                the schema for the ``DataFrame``.

        Returns: ``self``, to enable method chaining.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError("must be implemented by child class")

    def _set_reader_coord(
        self,
        sr: clib.SOMAArray,
        dim_idx: int,
        dim: pa.Field,
        coord: object,
    ) -> bool:
        if coord is None:
            return True  # No constraint; select all in this dimension

        if isinstance(coord, (str, bytes)):
            sr.set_dim_points_string_or_bytes(dim.name, [coord])
            return True

        if isinstance(coord, (pa.Array, pa.ChunkedArray)):
            # sr.set_dim_points_arrow does type disambiguation based on array schema -- so we do
            # not.
            sr.set_dim_points_arrow(dim.name, coord)
            return True

        if isinstance(coord, (Sequence, np.ndarray)):
            if self._set_reader_coord_by_py_seq_or_np_array(sr, dim_idx, dim, coord):
                return True

        if isinstance(coord, slice):
            _util.validate_slice(coord)
            if coord.start is None and coord.stop is None:
                return True

        if isinstance(coord, slice):
            _util.validate_slice(coord)
            if self._set_reader_coord_by_numeric_slice(sr, dim_idx, dim, coord):
                return True

        domain = self.domain[dim_idx]

        # Note: slice(None, None) matches the is_slice_of part, unless we also check the dim-type
        # part.
        if (
            is_slice_of(coord, str) or is_slice_of(coord, bytes)
        ) and _util.pa_types_is_string_or_bytes(dim.type):
            _util.validate_slice(coord)
            # Figure out which one.
            dim_type: Union[Type[str], Type[bytes]] = type(domain[0])
            # A ``None`` or empty start is always equivalent to empty str/bytes.
            start = coord.start or dim_type()
            if coord.stop is None:
                # There's no way to specify "to infinity" for strings.
                # We have to get the nonempty domain and use that as the end.
                ned = self._handle.non_empty_domain()
                _, stop = ned[dim_idx]
            else:
                stop = coord.stop
            sr.set_dim_ranges_string_or_bytes(dim.name, [(start, stop)])
            return True

        # Note: slice(None, None) matches the is_slice_of part, unless we also check the dim-type
        # part.
        if is_slice_of(coord, np.datetime64) and pa.types.is_timestamp(dim.type):
            _util.validate_slice(coord)
            # These timestamp types are stored in Arrow as well as TileDB as 64-bit integers (with
            # distinguishing metadata of course). For purposes of the query logic they're just
            # int64.
            istart = coord.start or domain[0]
            istart = int(istart.astype("int64"))
            istop = coord.stop or domain[1]
            istop = int(istop.astype("int64"))
            sr.set_dim_ranges_int64(dim.name, [(istart, istop)])
            return True

        if super()._set_reader_coord(sr, dim_idx, dim, coord):
            return True

        return False

    def _set_reader_coord_by_py_seq_or_np_array(
        self,
        sr: clib.SOMAArray,
        dim_idx: int,
        dim: pa.Field,
        coord: object,
    ) -> bool:
        if isinstance(coord, np.ndarray):
            if coord.ndim != 1:
                raise ValueError(
                    f"only 1D numpy arrays may be used to index; got {coord.ndim}"
                )

        try:
            set_dim_points = getattr(sr, f"set_dim_points_{dim.type}")
        except AttributeError:
            # We have to handle this type specially below
            pass
        else:
            set_dim_points(dim.name, coord)
            return True

        if _util.pa_types_is_string_or_bytes(dim.type):
            sr.set_dim_points_string_or_bytes(dim.name, coord)
            return True
        elif pa.types.is_timestamp(dim.type):
            if not isinstance(coord, (tuple, list, np.ndarray)):
                raise ValueError(
                    f"unhandled coord type {type(coord)} for index column named {dim.name}"
                )
            icoord = [
                int(e.astype("int64")) if isinstance(e, np.datetime64) else e
                for e in coord
            ]
            sr.set_dim_points_int64(dim.name, icoord)
            return True

        # TODO: bool

        raise ValueError(f"unhandled type {dim.type} for index column named {dim.name}")

    def _set_reader_coord_by_numeric_slice(
        self, sr: clib.SOMAArray, dim_idx: int, dim: pa.Field, coord: Slice[Any]
    ) -> bool:
        try:
            lo_hi = _util.slice_to_numeric_range(coord, self.domain[dim_idx])
        except _util.NonNumericDimensionError:
            return False  # We only handle numeric dimensions here.

        if not lo_hi:
            return True

        try:
            set_dim_range = getattr(sr, f"set_dim_ranges_{dim.type}")
            set_dim_range(dim.name, [lo_hi])
            return True
        except AttributeError:
            return False
