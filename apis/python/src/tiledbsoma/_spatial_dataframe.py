# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Implementation of a base class shared between GeometryDataFrame and PointCloudDataFrame
"""
from typing import Any, Optional, Sequence, Tuple, Union

import pyarrow as pa
import somacore
from somacore import CoordinateSpace, CoordinateTransform, options
from typing_extensions import Self

from ._read_iters import TableReadIter
from ._soma_array import SOMAArray

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
