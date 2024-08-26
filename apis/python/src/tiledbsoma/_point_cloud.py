# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Implementation of a SOMA Point Cloud
"""
from typing import Any, Optional, Sequence, Tuple, Union, cast

import pyarrow as pa
import somacore
from somacore import options
from typing_extensions import Self

from . import _arrow_types, _util
from . import pytiledbsoma as clib
from ._constants import SOMA_COORDINATE_SPACE_METADATA_KEY, SOMA_JOINID
from ._coordinates import Axis, CoordinateSpace
from ._dataframe import (
    Domain,
    _canonicalize_schema,
    _fill_out_slot_domain,
    _find_extent_for_domain,
)
from ._exception import SOMAError, map_exception_for_create
from ._query_condition import QueryCondition
from ._read_iters import TableReadIter
from ._spatial_dataframe import SpatialDataFrame
from ._tdb_handles import PointCloudWrapper
from ._types import OpenTimestamp
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context
from .options._tiledb_create_write_options import (
    TileDBCreateOptions,
    TileDBWriteOptions,
)

_UNBATCHED = options.BatchSize()


class PointCloud(SpatialDataFrame, somacore.PointCloud):

    __slots__ = ("_axis_names", "_coord_space")
    _wrapper_type = PointCloudWrapper

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
        for column_name in axis_names:
            if column_name not in index_column_names:
                raise ValueError(f"Spatial column '{column_name}' must an index column")
            column_dtype = schema.field(column_name).type
            if not (
                pa.types.is_integer(column_dtype) or pa.types.is_floating(column_dtype)
            ):
                raise ValueError(
                    f"Spatial column '{column_name}' must have an integer or "
                    f"floating-point type. Column type is {column_dtype!r}"
                )
        coord_space = CoordinateSpace(
            tuple(Axis(axis_name) for axis_name in axis_names)
        )
        context = _validate_soma_tiledb_context(context)
        schema = _canonicalize_schema(schema, index_column_names)
        if domain is None:
            domain = tuple(None for _ in index_column_names)
        else:
            ndom = len(domain)
            nidx = len(index_column_names)
            if ndom != nidx:
                raise ValueError(
                    f"if domain is specified, it must have the same length as index_column_names; got {ndom} != {nidx}"
                )

        index_column_schema = []
        index_column_data = {}

        for index_column_name, slot_domain in zip(index_column_names, domain):
            pa_field = schema.field(index_column_name)
            dtype = _arrow_types.tiledb_type_from_arrow_type(
                pa_field.type, is_indexed_column=True
            )

            slot_domain = _fill_out_slot_domain(
                slot_domain, index_column_name, pa_field.type, dtype
            )

            extent = _find_extent_for_domain(
                index_column_name,
                TileDBCreateOptions.from_platform_config(platform_config),
                dtype,
                slot_domain,
            )

            index_column_schema.append(pa_field)
            index_column_data[pa_field.name] = [*slot_domain, extent]

        index_column_info = pa.RecordBatch.from_pydict(
            index_column_data, schema=pa.schema(index_column_schema)
        )

        plt_cfg = _util.build_clib_platform_config(platform_config)
        timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
        try:
            clib.SOMAPointCloud.create(
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
        handle.meta[SOMA_COORDINATE_SPACE_METADATA_KEY] = coord_space.to_json()
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    def __init__(
        self,
        handle: PointCloudWrapper,
        **kwargs: Any,
    ):
        super().__init__(handle, **kwargs)

        # Get and validate coordinate space.
        try:
            coord_space = self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY]
        except KeyError as ke:
            raise SOMAError("Missing axis name metadata") from ke
        self._coord_space = CoordinateSpace.from_json(coord_space)
        for name in self._coord_space.axis_names:
            if name not in self.index_column_names:
                raise SOMAError(
                    f"Point cloud axis '{name}' does not match any of the index column"
                    f" names."
                )

    def __len__(self) -> int:
        """Returns the number of rows in the dataframe. Same as ``df.count``."""
        return self.count

    @property
    def axis_names(self) -> Tuple[str, ...]:
        return self._coord_space.axis_names

    @property
    def coordinate_space(self) -> CoordinateSpace:
        """Coordinate system for this scene."""
        return self._coord_space

    @coordinate_space.setter
    def coordinate_space(self, value: CoordinateSpace) -> None:
        if self._coord_space is not None:
            if value.axis_names != self._coord_space.axis_names:
                raise ValueError(
                    f"Cannot change axis names of a point cloud. Existing axis names "
                    f"are {self._coord_space.axis_names}. New coordinate space has "
                    f"axis names {self._coord_space.axis_names}."
                )
        self.metadata[SOMA_COORDINATE_SPACE_METADATA_KEY] = value.to_json()
        self._coord_space = value

    @property
    def count(self) -> int:
        """Returns the number of rows in the dataframe. Same as ``len(df)``.

        Lifecycle:
            Experimental.
        """
        self._check_open_read()
        # if is it in read open mode, then it is a DataFrameWrapper
        return cast(PointCloudWrapper, self._handle).count

    # Note: This will be updated to support spatial reads.
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
        """TODO: Add docs

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

        sr = clib.SOMAPointCloud.open(
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

    def write(
        self, values: pa.Table, platform_config: Optional[options.PlatformConfig] = None
    ) -> Self:
        """TODO: Add docs

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
