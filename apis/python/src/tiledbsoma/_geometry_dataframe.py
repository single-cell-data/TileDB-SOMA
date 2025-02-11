# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""
Implementation of a SOMA Geometry DataFrame
"""

from __future__ import annotations

import warnings
from typing import Any, Sequence, Tuple, Union

import pyarrow as pa
import somacore
from somacore import CoordinateSpace, CoordinateTransform, options
from typing_extensions import Self

from tiledbsoma._tdb_handles import GeometryDataFrameWrapper
from tiledbsoma.options._soma_tiledb_context import _validate_soma_tiledb_context
from tiledbsoma.options._tiledb_create_write_options import TileDBCreateOptions

from . import _arrow_types, _util
from . import pytiledbsoma as clib
from ._constants import (
    SOMA_COORDINATE_SPACE_METADATA_KEY,
    SOMA_GEOMETRY,
    SOMA_JOINID,
    SPATIAL_DISCLAIMER,
)
from ._dataframe import (
    Domain,
    _canonicalize_schema,
    _fill_out_slot_soma_domain,
    _find_extent_for_domain,
    _revise_domain_for_extent,
)
from ._exception import SOMAError, map_exception_for_create
from ._read_iters import TableReadIter
from ._spatial_dataframe import SpatialDataFrame
from ._spatial_util import (
    coordinate_space_from_json,
    coordinate_space_to_json,
)
from ._types import OpenTimestamp
from .options import SOMATileDBContext

_UNBATCHED = options.BatchSize()


class GeometryDataFrame(SpatialDataFrame, somacore.GeometryDataFrame):
    """A specialized SOMA object for storing complex geometries with spatial indexing.

    The ``GeometryDataFrame`` class is designed to store and manage geometric shapes
    such as polygons, lines, and multipoints, along with additional columns for
    associated attributes.

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_coord_space",)
    _wrapper_type = GeometryDataFrameWrapper

    # Lifecycle

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        coordinate_space: Union[Sequence[str], CoordinateSpace] = ("x", "y"),
        domain: Domain | None = None,
        platform_config: options.PlatformConfig | None = None,
        context: SOMATileDBContext | None = None,
        tiledb_timestamp: OpenTimestamp | None = None,
    ) -> Self:
        """Creates a new ``GeometryDataFrame`` at the given URI.

        The schema of the created geometry dataframe will include a column named
        ``soma_joinid`` of type ``pyarrow.int64``, with negative values
        disallowed, and a column named ``soma_geometry of type ``pyarrow.binary`` or
        ``pyarrow.large_binary``.  If a ``soma_joinid`` column or ``soma_geometry``
        are present in the provided schema, they must be of the correct type.  If
        either the ``soma_joinid`` column or ``soma_geometry`` column are not provided,
        one will be added.

        Args:
            uri: The URI where the dataframe will be created.
            schema: Arrow schema defining the per-column schema. This schema
                must define all columns, including columns to be named as index
                columns.  If the schema includes types unsupported by the SOMA
                implementation, a ValueError will be raised.
            coordinate_space: Either the coordinate space or the axis names for the
                coordinate space the point cloud is defined on.
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

        # Get coordinate space axis data.
        if isinstance(coordinate_space, CoordinateSpace):
            axis_names = tuple(axis.name for axis in coordinate_space)
            axis_units = tuple(axis.unit for axis in coordinate_space)
        else:
            axis_names = tuple(coordinate_space)
            axis_units = tuple(len(axis_names) * [None])

        index_column_names = (
            SOMA_GEOMETRY,
            SOMA_JOINID,
        )

        context = _validate_soma_tiledb_context(context)
        schema = _canonicalize_schema(
            schema, index_column_names, [SOMA_JOINID, SOMA_GEOMETRY]
        )

        # SOMA-to-core mappings:
        #
        # Before the current-domain feature was enabled (possible after core 2.25):
        #
        # * SOMA domain <-> core domain, AKA "max domain" which is a name we'll use for clarity
        # * core current domain did not exist
        #
        # After the current-domain feature was enabled:
        #
        # * SOMA max_domain <-> core domain
        # * SOMA domain <-> core current domain
        #
        # As far as the user is concerned, the SOMA-level domain is the only
        # thing they see and care about. Before 2.25 support, it was immutable
        # (since it was implemented by core domain). After 2.25 support, it is
        # mutable/up-resizeable (since it is implemented by core current domain).

        # At this point shift from API terminology "domain" to specifying a soma_ or core_
        # prefix for these variables. This is crucial to avoid developer confusion.
        soma_domain = domain
        domain = None

        # Check if domain has the right size (number of axis + 1 for SOMA_JOINID)

        if soma_domain is None:
            soma_domain = tuple(None for _ in index_column_names)
        else:
            ndom = len(soma_domain)
            nidx = len(index_column_names)
            if ndom != nidx:
                raise ValueError(
                    f"if domain is specified, it must have the same length as "
                    f"index_column_names; got {ndom} != {nidx}"
                )

        mutable_soma_domain = list(soma_domain)
        soma_geometry_domain = mutable_soma_domain[
            index_column_names.index(SOMA_GEOMETRY)
        ]

        if soma_geometry_domain is None:
            soma_geometry_domain = [None for _ in axis_names]
        elif not isinstance(soma_geometry_domain, list):
            raise ValueError(
                f"'{SOMA_GEOMETRY}' domain should be a list of tuple[float, float]"
            )
        elif len(soma_geometry_domain) != len(axis_names):
            raise ValueError(
                f"Dimension mishmatch between '{SOMA_GEOMETRY}' domain and coordinate system"
            )

        mutable_soma_domain[index_column_names.index(SOMA_GEOMETRY)] = (
            soma_geometry_domain
        )
        soma_domain = tuple(mutable_soma_domain)

        index_column_schema = []
        index_column_data = {}

        for index_column_name, slot_soma_domain in zip(index_column_names, soma_domain):
            pa_field = schema.field(index_column_name)
            dtype = _arrow_types.tiledb_type_from_arrow_type(
                pa_field.type, is_indexed_column=True
            )

            (slot_core_current_domain, saturated_cd) = _fill_out_slot_soma_domain(
                slot_soma_domain, False, index_column_name, pa_field.type, dtype
            )

            if index_column_name == SOMA_GEOMETRY:
                (slot_core_max_domain, saturated_md) = _fill_out_slot_soma_domain(
                    [None for _ in axis_names],
                    True,
                    index_column_name,
                    pa_field.type,
                    dtype,
                )
            else:
                (slot_core_max_domain, saturated_md) = _fill_out_slot_soma_domain(
                    None, True, index_column_name, pa_field.type, dtype
                )

            extent = _find_extent_for_domain(
                index_column_name,
                TileDBCreateOptions.from_platform_config(platform_config),
                dtype,
                slot_core_current_domain,
            )

            # Necessary to avoid core array-creation error "Reduce domain max by
            # 1 tile extent to allow for expansion."
            slot_core_current_domain = _revise_domain_for_extent(
                slot_core_current_domain, extent, saturated_cd
            )
            slot_core_max_domain = _revise_domain_for_extent(
                slot_core_max_domain, extent, saturated_md
            )

            # Here is our Arrow data API for communicating schema info between
            # Python/R and C++ libtiledbsoma:
            #
            # [0] core max domain lo
            # [1] core max domain hi
            # [2] core extent parameter
            # If present, these next two signal to use the current-domain feature:
            # [3] core current domain lo
            # [4] core current domain hi
            if index_column_name == SOMA_GEOMETRY:
                # SOMA_GEOMETRY has a specific schema
                index_column_schema.append(
                    pa.field(
                        pa_field.name,
                        pa.struct({axis: pa.float64() for axis in axis_names}),
                    )
                )
                index_column_data[pa_field.name] = [
                    [
                        (axis, slot_core_max_domain[0][idx])
                        for idx, axis in enumerate(axis_names)
                    ],
                    [
                        (axis, slot_core_max_domain[1][idx])
                        for idx, axis in enumerate(axis_names)
                    ],
                    [(axis, extent) for axis in axis_names],
                    [
                        (axis, slot_core_current_domain[0][idx])
                        for idx, axis in enumerate(axis_names)
                    ],
                    [
                        (axis, slot_core_current_domain[1][idx])
                        for idx, axis in enumerate(axis_names)
                    ],
                ]
            else:
                index_column_schema.append(pa_field)
                index_column_data[pa_field.name] = [
                    *slot_core_max_domain,
                    extent,
                    *slot_core_current_domain,
                ]

        index_column_info = pa.RecordBatch.from_pydict(
            index_column_data, schema=pa.schema(index_column_schema)
        )

        plt_cfg = _util.build_clib_platform_config(platform_config)
        timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
        try:
            clib.SOMAGeometryDataFrame.create(
                uri,
                schema=schema,
                index_column_info=index_column_info,
                axis_names=axis_names,
                axis_units=axis_units,
                ctx=context.native_context,
                platform_config=plt_cfg,
                timestamp=(0, timestamp_ms),
            )
        except SOMAError as e:
            raise map_exception_for_create(e, uri) from None

        return cls(
            cls._wrapper_type.open(uri, "w", context, tiledb_timestamp),
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    def __init__(
        self,
        handle: GeometryDataFrameWrapper,
        **kwargs: Any,
    ):
        super().__init__(handle, **kwargs)

        # Get and validate coordinate space.
        try:
            coord_space = self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY]
        except KeyError as ke:
            raise SOMAError("Missing coordinate space metadata") from ke
        self._coord_space = coordinate_space_from_json(coord_space)

    # Data operations

    def read(
        self,
        coords: options.SparseDFCoords = (),
        column_names: Sequence[str] | None = None,
        *,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
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
        region: options.SpatialRegion | None = None,
        column_names: Sequence[str] | None = None,
        *,
        region_transform: CoordinateTransform | None = None,
        region_coord_space: CoordinateSpace | None = None,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        value_filter: str | None = None,
        platform_config: options.PlatformConfig | None = None,
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
        platform_config: options.PlatformConfig | None = None,
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
    def axis_names(self) -> Tuple[str, ...]:
        """The names of the axes of the coordinate space the data is defined on.

        Lifecycle:
            Experimental.
        """
        return self._coord_space.axis_names

    @property
    def coordinate_space(self) -> CoordinateSpace:
        """Coordinate space for this geometry dataframe.

        Lifecycle:
            Experimental.
        """
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        """Coordinate space for this geometry dataframe.

        Lifecycle:
            Experimental.
        """
        if self._coord_space is not None:
            if value.axis_names != self._coord_space.axis_names:
                raise ValueError(
                    f"Cannot change axis names of a geometry dataframe. Existing "
                    f"axis names are {self._coord_space.axis_names}. New coordinate "
                    f"space has axis names {value.axis_names}."
                )
        self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY] = coordinate_space_to_json(
            value
        )
        self._coord_space = value
