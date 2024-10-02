# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.
"""
Implementation of a SOMA Point Cloud DataFrame
"""

import warnings
from typing import Any, Optional, Sequence, Tuple, Union, cast

import pyarrow as pa
import somacore
from somacore import Axis, CoordinateSpace, CoordinateTransform, options
from typing_extensions import Self

from . import _arrow_types, _util
from . import pytiledbsoma as clib
from ._constants import (
    SOMA_COORDINATE_SPACE_METADATA_KEY,
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
from ._flags import NEW_SHAPE_FEATURE_FLAG_ENABLED
from ._query_condition import QueryCondition
from ._read_iters import TableReadIter
from ._spatial_dataframe import SpatialDataFrame
from ._spatial_util import (
    coordinate_space_from_json,
    coordinate_space_to_json,
    process_spatial_df_region,
)
from ._tdb_handles import PointCloudDataFrameWrapper
from ._types import OpenTimestamp
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context
from .options._tiledb_create_write_options import (
    TileDBCreateOptions,
    TileDBWriteOptions,
)

_UNBATCHED = options.BatchSize()


class PointCloudDataFrame(SpatialDataFrame, somacore.PointCloudDataFrame):
    """A specialized SOMA DataFrame for storing collections of points in
    multi-dimensional space.

    The ``PointCloudDataFrame`` class is designed to efficiently store and query point data,
    where each point is represented by coordinates in one or more spatial dimensions
    (e.g., x, y, z) and may have additional columns for associated attributes.

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_coord_space",)
    _wrapper_type = PointCloudDataFrameWrapper

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID, "x", "y"),
        axis_names: Sequence[str] = ("x", "y"),
        domain: Optional[Domain] = None,
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
    ) -> Self:
        """Creates a new ``PointCloudDataFrame`` at the given URI.

        The schema of the created point cloud dataframe will include a column named
        ``soma_joinid`` of type ``pyarrow.int64``, with negative values disallowed, and
        at least one axis with numeric type.  If a ``soma_joinid`` column is
        present in the provided schema, it must be of the correct type.  If the
        ``soma_joinid`` column is not provided, one will be added. The ``soma_joinid``
        may be an index column. The axis columns must be index columns.

        Args:
            uri: The URI where the dataframe will be created.
            schema: Arrow schema defining the per-column schema. This schema
                must define all columns, including columns to be named as index
                columns.  If the schema includes types unsupported by the SOMA
                implementation, a ValueError will be raised.
            index_column_names: A list of column names to use as user-defined index
                columns (e.g., ``['x', 'y']``). All named columns must exist in the
                schema, and at least one index column name is required.
                Default is ``("soma_joinid", "x", "y")``.
            axis_names: An ordered list of axis column names that correspond to the
                names of axes of the the coordinate space the points are defined on.
                Must be the name of index columns. Default is ``("x", "y")``.
            domain: An optional sequence of tuples specifying the domain of each
                index column. Each tuple should be a pair consisting of the minimum
                and maximum values storable in the index column. If omitted entirely,
                or if ``None`` in a given dimension, the corresponding index-column
                domain will use the minimum and maximum possible values for the
                column's datatype.  This makes a point cloud dataframe growable.

        Returns:
            The newly created point cloud, opened for writing.

        Lifecycle:
            Experimental.
        """
        warnings.warn(SPATIAL_DISCLAIMER)

        axis_dtype: Optional[pa.DataType] = None
        for column_name in axis_names:
            if column_name not in index_column_names:
                raise ValueError(f"Spatial column '{column_name}' must an index column")
            # Check axis column type is valid and all axis columns have the same type.
            if axis_dtype is None:
                axis_dtype = schema.field(column_name).type
                if not (
                    pa.types.is_integer(axis_dtype) or pa.types.is_floating(axis_dtype)
                ):
                    raise ValueError(
                        f"Spatial column '{column_name}' must have an integer or "
                        f"floating-point type. Column type is {axis_dtype!r}"
                    )
            else:
                column_dtype = schema.field(column_name).type
                if column_dtype != axis_dtype:
                    raise ValueError("All spatial axes must have the same datatype.")

        # mypy false positive https://github.com/python/mypy/issues/5313
        coord_space = CoordinateSpace(
            tuple(Axis(axis_name) for axis_name in axis_names)  # type: ignore
        )
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

        for index_column_name, slot_soma_domain in zip(index_column_names, soma_domain):
            pa_field = schema.field(index_column_name)
            dtype = _arrow_types.tiledb_type_from_arrow_type(
                pa_field.type, is_indexed_column=True
            )

            (slot_core_current_domain, saturated_cd) = _fill_out_slot_soma_domain(
                slot_soma_domain, index_column_name, pa_field.type, dtype
            )
            (slot_core_max_domain, saturated_md) = _fill_out_slot_soma_domain(
                None, index_column_name, pa_field.type, dtype
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

        plt_cfg = _util.build_clib_platform_config(platform_config)
        timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
        try:
            clib.SOMAPointCloudDataFrame.create(
                uri,
                schema=schema,
                index_column_info=index_column_info,
                ctx=context.native_context,
                platform_config=plt_cfg,
                timestamp=(0, timestamp_ms),
            )
        except SOMAError as e:
            raise map_exception_for_create(e, uri) from None

        handle = cls._wrapper_type.open(uri, "w", context, tiledb_timestamp)
        handle.meta[SOMA_COORDINATE_SPACE_METADATA_KEY] = coordinate_space_to_json(
            coord_space
        )
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    def __init__(
        self,
        handle: PointCloudDataFrameWrapper,
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
            if name not in self.index_column_names:
                raise SOMAError(
                    f"Point cloud dataframe axis '{name}' does not match any of the "
                    f"index column names."
                )

    # Data operations

    def __len__(self) -> int:
        """Returns the number of rows in the point cound dataframe."""
        return self.count

    @property
    def count(self) -> int:
        """Returns the number of rows in the dataframe."""
        self._check_open_read()
        # if is it in read open mode, then it is a PointCloudDataFrameWrapper
        return cast(PointCloudDataFrameWrapper, self._handle).count

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

        sr = clib.SOMAPointCloudDataFrame.open(
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
                # mypy false positive https://github.com/python/mypy/issues/5313
                region_coord_space = CoordinateSpace(
                    tuple(Axis(axis_name) for axis_name in region_transform.input_axes)  # type: ignore
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
        coords, data_region, inv_transform = process_spatial_df_region(
            region,
            region_transform,
            dict(),  #  Move index value_filters into this dict to optimize queries
            self._tiledb_dim_names(),
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
        _util.check_type("values", values, (pa.Table,))

        write_options: Union[TileDBCreateOptions, TileDBWriteOptions]
        sort_coords = None
        if isinstance(platform_config, TileDBCreateOptions):
            raise ValueError(
                "As of TileDB-SOMA 1.13, the write method takes "
                "TileDBWriteOptions instead of TileDBCreateOptions"
            )
        write_options = TileDBWriteOptions.from_platform_config(platform_config)
        sort_coords = write_options.sort_coords

        clib_dataframe = self._handle._handle

        for batch in values.to_batches():
            clib_dataframe.write(batch, sort_coords or False)

        if write_options.consolidate_and_vacuum:
            clib_dataframe.consolidate_and_vacuum()

        return self

    # Metadata operations

    @property
    def coordinate_space(self) -> CoordinateSpace:
        """Coordinate space for this point cloud dataframe.

        Lifecycle:
            Experimental.
        """
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        """Coordinate space for this point cloud dataframe.

        Lifecycle:
            Experimental.
        """
        if self._coord_space is not None:
            if value.axis_names != self._coord_space.axis_names:
                raise ValueError(
                    f"Cannot change axis names of a point cloud dataframe. Existing "
                    f"axis names are {self._coord_space.axis_names}. New coordinate "
                    f"space has axis names {value.axis_names}."
                )
        self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY] = coordinate_space_to_json(
            value
        )
        self._coord_space = value

    @property
    def axis_names(self) -> Tuple[str, ...]:
        """The names of the axes of the coordinate space the data is defined on.

        Lifecycle:
            Experimental.
        """
        return self._coord_space.axis_names
