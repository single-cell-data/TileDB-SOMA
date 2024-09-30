# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.
"""
Implementation of a SOMA Geometry DataFrame
"""

import warnings
from typing import Any, Optional, Sequence, Tuple, Union

import pyarrow as pa
import somacore
from somacore import CoordinateSpace, CoordinateTransform, options
from typing_extensions import Self

from ._constants import SOMA_GEOMETRY, SOMA_JOINID, SPATIAL_DISCLAIMER
from ._dataframe import Domain
from ._read_iters import TableReadIter
from ._types import OpenTimestamp
from .options import SOMATileDBContext

_UNBATCHED = options.BatchSize()


class GeometryDataFrame(somacore.GeometryDataFrame):
    """A specialized SOMA object for storing complex geometries with spatial indexing.

    The ``GeometryDataFrame`` class is designed to store and manage geometric shapes
    such as polygons, lines, and multipoints, along with additional columns for
    associated attributes.

    Lifecycle:
        Experimental.
    """

    __slots__ = ()

    # Lifecycle

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID, SOMA_GEOMETRY),
        axis_names: Sequence[str] = ("x", "y"),
        domain: Optional[Domain] = None,
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
    ) -> Self:
        """Creates a new ``GeometryDataFrame`` at the given URI.

        The schema of the created geometry dataframe will include a column named
        ``soma_joinid`` of type ``pyarrow.int64``, with negative values
        disallowed, and a column named ``soma_geometry of type ``pyarrow.binary`` or
        ``pyarrow.large_binary``.  If a ``soma_joinid`` column or ``soma_geometry``
        are present in the provided schema, they must be of the correct type.  If
        either the ``soma_joinid`` column or ``soma_geometry`` column are not provided,
        one will be added. The ``soma_joinid`` may be an index column. The
        ``soma_geometry`` column must be an index column.

        Args:
            uri: The URI where the dataframe will be created.
            schema: Arrow schema defining the per-column schema. This schema
                must define all columns, including columns to be named as index
                columns.  If the schema includes types unsupported by the SOMA
                implementation, a ValueError will be raised.
            index_column_names: A list of column names to use as user-defined
                index columns (e.g., ``['cell_type', 'tissue_type']``).
                All named columns must exist in the schema, and at least one
                index column name is required.
            axis_names: An ordered list of axis column names that correspond to the
                names of the axes of the coordinate space the geometries are defined
                on.
            domain: An optional sequence of tuples specifying the domain of each
                index column. Two tuples must be provided for the ``soma_geometry``
                column which store the width followed by the height. Each tuple should
                be a pair consisting of the minimum and maximum values storable in the
                index column. If omitted entirely, or if ``None`` in a given dimension,
                the corresponding index-column domain will use the minimum and maximum
                possible values for the column's datatype.  This makes a dataframe
                growable.

        Returns:
            The newly created geometry dataframe, opened for writing.

        Lifecycle:
            Experimental.
        """
        warnings.warn(SPATIAL_DISCLAIMER)
        raise NotImplementedError()

    # Data operations

    def read(
        self,
        coords: options.SparseDFCoords = (),
        column_names: Optional[Sequence[str]] = None,
        *,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: Optional[options.ReadPartitions] = None,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        value_filter: Optional[str] = None,
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
        raise NotImplementedError()

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
            region_transform: An optional coordinate transform from the read region to the
                coordinate system of the spatial dataframe.
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
        raise NotImplementedError()

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
        raise NotImplementedError()

    # Metadata operations

    @property
    def schema(self) -> pa.Schema:
        """The schema of the data in this dataframe.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def index_column_names(self) -> Tuple[str, ...]:
        """The names of the index (dimension) columns.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def axis_names(self) -> Tuple[str, ...]:
        """The names of the axes of the coordinate space the data is defined on.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def coordinate_space(self) -> CoordinateSpace:
        """Coordinate space for this geometry dataframe.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        """Coordinate space for this geometry dataframe.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()

    @property
    def domain(self) -> Tuple[Tuple[Any, Any], ...]:
        """The allowable range of values in each index column.

        Returns: a tuple of minimum and maximum values, inclusive,
            storable on each index column of the dataframe.

        Lifecycle:
            Experimental.
        """
        raise NotImplementedError()
