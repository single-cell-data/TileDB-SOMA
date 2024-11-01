# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.
"""
Implementation of a SOMA Geometry DataFrame
"""

import sys
import warnings
from typing import Any, Mapping, Optional, Sequence, Tuple, Union

import pyarrow as pa
import somacore
from somacore import CoordinateSpace, CoordinateTransform, options
from typing_extensions import Self

from tiledbsoma import _arrow_types, _util
from tiledbsoma._dataframe import (
    _fill_out_slot_soma_domain,
    _find_extent_for_domain,
    _revise_domain_for_extent,
)
from tiledbsoma._exception import SOMAError, map_exception_for_create
from tiledbsoma._flags import NEW_SHAPE_FEATURE_FLAG_ENABLED
from tiledbsoma._query_condition import QueryCondition
from tiledbsoma._spatial_dataframe import SpatialDataFrame
from tiledbsoma._spatial_util import (
    coordinate_space_from_json,
    coordinate_space_to_json,
    process_spatial_geometry_df_region,
)
from tiledbsoma._tdb_handles import GeometryDataFrameWrapper
from tiledbsoma.options._soma_tiledb_context import _validate_soma_tiledb_context
from tiledbsoma.options._tiledb_create_write_options import TileDBCreateOptions

from . import pytiledbsoma as clib
from ._constants import (
    SOMA_COORDINATE_SPACE_METADATA_KEY,
    SOMA_GEOMETRY,
    SOMA_JOINID,
    SPATIAL_DISCLAIMER,
)
from ._dataframe import AxisDomain, Domain
from ._read_iters import TableReadIter
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
        index_column_names: Optional[Sequence[str]] = (SOMA_JOINID, SOMA_GEOMETRY),
        coordinate_space: Union[Sequence[str], CoordinateSpace] = ("x", "y"),
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
        if not isinstance(coordinate_space, CoordinateSpace):
            coordinate_space = CoordinateSpace.from_axis_names(coordinate_space)
        if index_column_names is None:
            index_column_names = (SOMA_JOINID, SOMA_GEOMETRY)

        context = _validate_soma_tiledb_context(context)
        schema = _canonicalize_schema(schema, index_column_names)

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

        if soma_domain is None:
            soma_domain = tuple(None for _ in index_column_names)
        else:
            ndom = len(soma_domain)
            nidx = len(index_column_names)
            if ndom != nidx:
                raise ValueError(
                    f"if domain is specified, it must have the same length as index_column_names; got {ndom} != {nidx}"
                )

        index_column_schema = []
        index_column_data = {}

        spatial_column_schema: Sequence[pa.Field] = []
        spatial_column_data: Mapping[str, Any] = {}

        for index_column_name, slot_soma_domain in zip(index_column_names, soma_domain):
            if index_column_name == SOMA_GEOMETRY:
                spatial_column_schema, spatial_column_data = _spatial_domain(
                    [axis.name for axis in coordinate_space],
                    slot_soma_domain,
                    platform_config,
                )

            pa_field = schema.field(index_column_name)
            dtype = _arrow_types.tiledb_type_from_arrow_type(
                pa_field.type, is_indexed_column=True
            )

            (slot_core_current_domain, saturated_cd) = _fill_out_slot_soma_domain(
                slot_soma_domain, False, index_column_name, pa_field.type, dtype
            )
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

            index_column_schema.append(pa_field)
            if NEW_SHAPE_FEATURE_FLAG_ENABLED:

                index_column_data[pa_field.name] = [
                    *slot_core_max_domain,
                    extent,
                    *slot_core_current_domain,
                ]

            else:
                index_column_data[pa_field.name] = [*slot_core_current_domain, extent]

        index_column_info = pa.RecordBatch.from_pydict(
            index_column_data, schema=pa.schema(index_column_schema)
        )

        spatial_column_info = pa.RecordBatch.from_pydict(
            spatial_column_data, schema=pa.schema(spatial_column_schema)
        )

        print(spatial_column_info)

        plt_cfg = _util.build_clib_platform_config(platform_config)
        timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
        try:
            clib.SOMAGeometryDataFrame.create(
                uri,
                schema=schema,
                index_column_info=index_column_info,
                spatial_column_info=spatial_column_info,
                ctx=context.native_context,
                platform_config=plt_cfg,
                timestamp=(0, timestamp_ms),
            )
        except SOMAError as e:
            raise map_exception_for_create(e, uri) from None

        handle = cls._wrapper_type.open(uri, "w", context, tiledb_timestamp)
        handle.meta[SOMA_COORDINATE_SPACE_METADATA_KEY] = coordinate_space_to_json(
            coordinate_space
        )
        return cls(
            handle,
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
        for name in self._coord_space.axis_names:
            if name not in self.axis_names:
                raise SOMAError(
                    f"Geometry dataframe axis '{name}' does not match any of the "
                    f"spatial column names."
                )

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
        del batch_size  # Currently unused.
        _util.check_unpartitioned(partitions)
        self._check_open_read()

        handle = self._handle._handle

        context = handle.context()
        if platform_config is not None:
            config = context.tiledb_config.copy()
            config.update(platform_config)
            context = clib.SOMAContext(config)

        sr = clib.SOMAGeometryDataFrame.open(
            uri=handle.uri,
            mode=clib.OpenMode.read,
            context=context,
            column_names=column_names or [],
            result_order=_util.to_clib_result_order(result_order),
            timestamp=handle.timestamp and (0, handle.timestamp),
        )

        if value_filter is not None:
            sr.set_condition(QueryCondition(value_filter), handle.schema)

        self._set_reader_coords(sr, coords)

        # # TODO: batch_size
        return TableReadIter(sr)

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
        # Set/check transform and region coordinate space.
        if region_transform is None:
            region_transform = somacore.IdentityTransform(
                self.axis_names, self.axis_names
            )
            if region_coord_space is not None:
                raise ValueError(
                    "Cannot specify the output coordinate space when region transform i"
                    "is ``None``."
                )
            region_coord_space = self._coord_space
        else:
            if region_coord_space is None:
                region_coord_space = CoordinateSpace.from_axis_names(
                    region_transform.input_axes
                )
            elif region_transform.input_axes != region_coord_space.axis_names:
                raise ValueError(
                    f"The input axes '{region_transform.input_axes}' of the region "
                    f"transform must match the axes '{region_coord_space.axis_names}' "
                    f"of the coordinate space the requested region is defined in."
                )
            if region_transform.output_axes != self._coord_space.axis_names:
                raise ValueError(
                    f"The output axes of '{region_transform.output_axes}' of the "
                    f"transform must match the axes '{self._coord_space.axis_names}' "
                    f"of the coordinate space of this point cloud dataframe."
                )

        # Process the user provided region.
        coords, data_region, inv_transform = process_spatial_geometry_df_region(
            region,
            region_transform,
            dict(),  # Move index value_filters into this dict to optimize queries
            self.index_column_names,
            self._coord_space.axis_names,
            self._handle.schema,
        )

        return somacore.SpatialRead(
            self.read(
                coords,
                column_names,
                result_order=result_order,
                value_filter=value_filter,
                batch_size=batch_size,
                partitions=partitions,
                platform_config=platform_config,
            ),
            self.coordinate_space,
            region_coord_space,
            inv_transform,
        )

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
    def index_column_names(self) -> Tuple[str, ...]:
        """The names of the index (dimension) columns.

        Lifecycle:
            Experimental.
        """
        return tuple(self._handle._handle.index_column_names)

    @property
    def axis_names(self) -> Tuple[str, ...]:
        """The names of the axes of the coordinate space the data is defined on.

        Lifecycle:
            Experimental.
        """
        return tuple(self._handle._handle.spatial_column_names)

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

    def _set_reader_coord(
        self, sr: clib.SOMAGeometryDataFrame, dim_idx: int, dim: pa.Field, coord: object
    ) -> bool:
        """Parses a single coordinate entry.

        The base implementation parses the most fundamental types shared by all
        TileDB Array types; subclasses can implement their own readers that
        handle types not recognized here.

        Returns:
            True if successful, False if unrecognized.
        """

        if coord is None:
            return True  # No constraint; select all in this dimension

        if dim.name == SOMA_GEOMETRY:
            if isinstance(coord, list):
                if len(coord) > len(self.axis_names):
                    # Got more ranges than the number of spatial axes
                    return False

                for axis, axis_coord in zip(self.axis_names, coord):
                    if axis_coord is None:
                        continue
                    if isinstance(axis_coord, float) or isinstance(axis_coord, int):
                        sr.set_spatial_dim_ranges(axis, [coord, coord])
                        continue
                    if isinstance(axis_coord, slice):
                        _util.validate_slice(axis_coord)
                        try:
                            # Domain sanitization is handled for each spatial axis internally
                            spatial_lo_hi = _util.slice_to_numeric_range(
                                axis_coord, (-sys.float_info.max, sys.float_info.max)
                            )
                        except _util.NonNumericDimensionError:
                            return False  # We only handle numeric dimensions here.
                        if spatial_lo_hi:
                            sr.set_spatial_dim_ranges(axis, [spatial_lo_hi])
                        # If `None`, coord was `slice(None)` and there is no constraint.
                        continue
                    return False
                return True
            return False

        if isinstance(coord, int):
            sr.set_dim_points_int64(dim.name, [coord])
            return True
        if isinstance(coord, slice):
            _util.validate_slice(coord)
            try:
                dom = self._handle.domain[dim_idx]
                lo_hi = _util.slice_to_numeric_range(coord, dom)
            except _util.NonNumericDimensionError:
                return False  # We only handle numeric dimensions here.
            if lo_hi:
                sr.set_dim_ranges_int64(dim.name, [lo_hi])
            # If `None`, coord was `slice(None)` and there is no constraint.
            return True
        return False


def _canonicalize_schema(
    schema: pa.Schema, index_column_names: Sequence[str]
) -> pa.Schema:
    """Turns an Arrow schema into the canonical version and checks for errors.

    Returns a schema, which may be modified by the addition of required columns
    (e.g. ``soma_joinid``, ``soma_geometry``).
    """
    _util.check_type("schema", schema, (pa.Schema,))
    if not index_column_names:
        raise ValueError("DataFrame requires one or more index columns")

    if SOMA_JOINID in schema.names:
        joinid_type = schema.field(SOMA_JOINID).type
        if joinid_type != pa.int64():
            raise ValueError(
                f"{SOMA_JOINID} field must be of type Arrow int64 but is {joinid_type}"
            )
    else:
        # add SOMA_JOINID
        schema = schema.append(pa.field(SOMA_JOINID, pa.int64()))

    if SOMA_GEOMETRY in schema.names:
        geometry_type = schema.field(SOMA_GEOMETRY).type
        if geometry_type != pa.binary():
            raise ValueError(
                f"{SOMA_GEOMETRY} field must be of type Arrow binary but is {geometry_type}"
            )
    else:
        # add SOMA_GEOMETRY
        schema = schema.append(
            pa.field(SOMA_GEOMETRY, pa.binary(), metadata={"dtype": "WKB"})
        )

    # verify no illegal use of soma_ prefix
    for field_name in schema.names:
        if (
            field_name.startswith("soma_")
            and field_name != SOMA_JOINID
            and field_name != SOMA_GEOMETRY
        ):
            raise ValueError(
                f"DataFrame schema may not contain fields with name prefix ``soma_``: got ``{field_name}``"
            )

    # verify that all index_column_names are present in the schema
    schema_names_set = set(schema.names)
    for index_column_name in index_column_names:
        if (
            index_column_name.startswith("soma_")
            and index_column_name != SOMA_JOINID
            and index_column_name != SOMA_GEOMETRY
        ):
            raise ValueError(
                f'index_column_name other than "soma_joinid" must not begin with "soma_"; got "{index_column_name}"'
            )
        if index_column_name not in schema_names_set:
            schema_names_string = "{}".format(list(schema_names_set))
            raise ValueError(
                f"All index names must be defined in the dataframe schema: '{index_column_name}' not in {schema_names_string}"
            )
        dtype = schema.field(index_column_name).type
        if not pa.types.is_dictionary(dtype) and dtype not in [
            pa.int8(),
            pa.uint8(),
            pa.int16(),
            pa.uint16(),
            pa.int32(),
            pa.uint32(),
            pa.int64(),
            pa.uint64(),
            pa.float32(),
            pa.float64(),
            pa.binary(),
            pa.large_binary(),
            pa.string(),
            pa.large_string(),
            pa.timestamp("s"),
            pa.timestamp("ms"),
            pa.timestamp("us"),
            pa.timestamp("ns"),
        ]:
            raise TypeError(
                f"Unsupported index type {schema.field(index_column_name).type}"
            )

    return schema


def _spatial_domain(
    spatial_column_names: Sequence[str],
    domain: AxisDomain,
    platform_config: options.PlatformConfig,
) -> Tuple[Sequence[pa.Field], Mapping[str, Any]]:
    spatial_column_schema = []
    spatial_column_data = {}

    if domain is None:
        spatial_domain: Domain = tuple(None for _ in spatial_column_names)
    else:
        spatial_domain = domain

    for spatial_column_name, slot_spatial_domain in zip(
        spatial_column_names, spatial_domain
    ):
        pa_field = pa.field(spatial_column_name, pa.float64())
        dtype = _arrow_types.tiledb_type_from_arrow_type(
            pa_field.type, is_indexed_column=True
        )

        (slot_core_current_domain, saturated_cd) = _fill_out_slot_soma_domain(
            slot_spatial_domain, False, spatial_column_name, pa_field.type, dtype
        )
        (slot_core_max_domain, saturated_md) = _fill_out_slot_soma_domain(
            None, True, spatial_column_name, pa_field.type, dtype
        )

        extent = _find_extent_for_domain(
            spatial_column_name,
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

        spatial_column_schema.append(pa_field)
        if NEW_SHAPE_FEATURE_FLAG_ENABLED:
            spatial_column_data[pa_field.name] = [
                *slot_core_max_domain,
                extent,
                *slot_core_current_domain,
            ]

        else:
            spatial_column_data[pa_field.name] = [*slot_core_current_domain, extent]

    return spatial_column_schema, spatial_column_data
