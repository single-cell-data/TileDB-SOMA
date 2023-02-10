import collections.abc
from typing import Any, Mapping, Optional, Sequence, Tuple, TypeVar, cast

import numpy as np
import pyarrow as pa
import somacore
import tiledb
from somacore import options
from typing_extensions import Self

from . import arrow_types, util
from .constants import SOMA_JOINID
from .options import SOMATileDBContext
from .options.tiledb_create_options import TileDBCreateOptions
from .query_condition import QueryCondition
from .read_iters import TableReadIter
from .tiledb_array import TileDBArray
from .types import NPFloating, NPInteger

Slice = TypeVar("Slice", bound=Sequence[int])
_UNBATCHED = options.BatchSize()


class DataFrame(TileDBArray, somacore.DataFrame):
    """
    Represents ``obs``, ``var``, and others.

    Every ``DataFrame`` must contain a column called ``soma_joinid``, of type ``int64``.
    The ``soma_joinid`` column contains a unique value for each row in the dataframe,
    and acts as a joint key for other objects, such as ``SparseNDArray``.

    [lifecycle: experimental]
    """

    __slots__ = ()

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID,),
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
    ) -> "DataFrame":
        """
        Creates the data structure on disk/S3/cloud.

        [lifecycle: experimental]

        :param schema: Arrow schema defining the per-column schema. This schema
            must define all columns, including columns to be named as index columns.
            If the schema includes types unsupported by the SOMA implementation,
            an error will be raised.

        :param index_column_names: A list of column names to use as user-defined
            index columns (e.g., ``['cell_type', 'tissue_type']``).
            All named columns must exist in the schema, and at least one
            index column name is required.

        :param platform_config: Platform-specific options used to create this DataFrame,
            provided via ``{"tiledb": {"create": ...}}`` nested keys.
        """
        context = context or SOMATileDBContext()
        schema = _canonicalize_schema(schema, index_column_names)
        tdb_schema = _build_tiledb_schema(
            schema,
            index_column_names,
            TileDBCreateOptions.from_platform_config(platform_config),
            context,
        )
        handle = cls._create_internal(uri, tdb_schema, context)
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    def keys(self) -> Tuple[str, ...]:
        """
        Returns the names of the columns when read back as a dataframe.

        [lifecycle: experimental]
        """
        return self._tiledb_array_keys()

    @property
    def index_column_names(self) -> Tuple[str, ...]:
        """
        Return index (dimension) column names.
        """
        return self._tiledb_dim_names()

    @property
    def count(self) -> int:
        """
        Return the number of rows in the dataframe. Same as `len(df)`.
        """
        # A.domain.shape at the tiledb level gives us the 0..2^63 range which is not what we want
        self._check_open_read()
        return cast(int, self._soma_reader().nnz())

    def __len__(self) -> int:
        """
        Return the number of rows in the dataframe. Same as `df.count`.
        """
        return self.count

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
        """
        Read a user-defined subset of data, addressed by the dataframe indexing columns,
        optionally filtered, and return results as one or more Arrow tables.

        [lifecycle: experimental]

        :param coords: for each index dimension, which rows to read.
            Defaults to ``None``, meaning no constraint -- all IDs.
        :param column_names: the named columns to read and return.
            Defaults to ``None``, meaning no constraint -- all column names.
        :param partitions: an optional ``ReadPartitions`` hint to indicate
            how results should be organized.
        :param result_order: order of read results.
            This can be one of 'row-major', 'col-major', or 'auto'.
        :param value_filter: an optional [value filter] to apply to the results.
            Defaults to no filter.

        **Indexing**: the ``coords`` parameter will support, per dimension: a list of values of the type of the indexed column.

        Acceptable ways to index
        ------------------------

        * A sequence of coordinates is accepted, one per dimension.
        * Sequence length must be <= number of dimensions.
        * If the sequence contains missing coordinates (length less than number of dimensions),
          then `slice(None)` -- i.e. no constraint -- is assumed for the missing dimensions.
        * Per-dimension, explicitly specified coordinates can be one of: None, a value, a
          list/ndarray/paarray/etc of values, a slice, etc.
        * Slices are doubly inclusive: slice(2,4) means [2,3,4] not [2,3].
          Slice steps are not supported.
          Slices can be `slice(None)`, meaning select all in that dimension, but may not be half-specified:
          `slice(2,None)` and `slice(None,4)` are both unsupported.
        * Negative indexing is unsupported.
        """
        del batch_size, platform_config  # Currently unused.
        util.check_unpartitioned(partitions)
        self._check_open_read()
        result_order = options.ResultOrder(result_order)

        schema = self._handle.schema
        query_condition = None
        if value_filter is not None:
            query_condition = QueryCondition(value_filter)

        sr = self._soma_reader(
            schema=schema,  # query_condition needs this
            column_names=column_names,
            query_condition=query_condition,
            result_order=result_order.value,
        )

        if not isinstance(coords, (list, tuple)):
            raise TypeError(
                f"coords type {type(coords)} unsupported; expected list or tuple"
            )
        if schema.domain.ndim < len(coords):
            raise ValueError(
                f"coords {coords} must have length between 1 and ndim ({schema.domain.ndim}); got {len(coords)}"
            )

        for i, dim_coords in enumerate(coords):
            # Example: coords = [None, 3, slice(4,5)]
            # dim_coords takes on values None, 3, and slice(4,5) in this loop body.
            dim_name = schema.domain.dim(i).name
            if dim_coords is None:
                pass  # No constraint; select all in this dimension
            elif isinstance(dim_coords, (int, str, bytes)):
                sr.set_dim_points(dim_name, [dim_coords])
            elif isinstance(dim_coords, np.ndarray):
                if dim_coords.ndim != 1:
                    raise ValueError(
                        f"only 1D numpy arrays may be used to index; got {dim_coords.ndim}"
                    )
                sr.set_dim_points(dim_name, dim_coords)
            elif isinstance(dim_coords, slice):
                ned = self._handle.reader.nonempty_domain()
                # ned is None iff the array has no data
                lo_hi = util.slice_to_range(dim_coords, ned[i]) if ned else None
                if lo_hi is not None:
                    lo, hi = lo_hi
                    if lo < 0 or hi < 0:
                        raise ValueError(
                            f"slice start and stop may not be negative; got ({lo}, {hi})"
                        )
                    if lo > hi:
                        raise ValueError(
                            f"coordinate at slot {i} must have lo <= hi; got {lo} > {hi}"
                        )
                    sr.set_dim_ranges(dim_name, [lo_hi])
                # Else, no constraint in this slot. This is `slice(None)` which is like
                # Python indexing syntax `[:]`.
            elif isinstance(
                dim_coords,
                (collections.abc.Sequence, pa.Array, pa.ChunkedArray),
            ):
                sr.set_dim_points(dim_name, dim_coords)
            else:
                raise TypeError(f"coords[{i}] type {type(dim_coords)} is unsupported")

        # TODO: platform_config
        # TODO: batch_size

        sr.submit()
        return TableReadIter(sr)

    def write(
        self, values: pa.Table, platform_config: Optional[Mapping[str, Any]] = None
    ) -> Self:
        """
        Write an Arrow table to the persistent object. As duplicate index values
        are not allowed, index values already present in the object are overwritten
        and new index values are added.

        [lifecycle: experimental]

        :param values: An Arrow table containing all columns, including
            the index columns. The schema for the values must match
            the schema for the ``DataFrame``.
        """
        util.check_type("values", values, (pa.Table,))

        del platform_config  # unused
        dim_cols_list = []
        attr_cols_map = {}
        dim_names_set = self.index_column_names
        n = None

        for name in values.schema.names:
            n = len(values.column(name))
            if name in dim_names_set:
                dim_cols_list.append(values.column(name).to_pandas())
            else:
                attr_cols_map[name] = values.column(name).to_pandas()
        if n is None:
            raise ValueError(f"did not find any column names in {values.schema.names}")

        dim_cols_tuple = tuple(dim_cols_list)
        self._handle.writer[dim_cols_tuple] = attr_cols_map

        return self


def _canonicalize_schema(
    schema: pa.Schema, index_column_names: Sequence[str]
) -> pa.Schema:
    """Turn a Arrow schema into the canonical version and check for errors.

    Returns a schema, which may be modified by the addition of required columns
    (e.g. ``soma_joinid``).
    """
    util.check_type("schema", schema, (pa.Schema,))
    if not index_column_names:
        raise ValueError("DataFrame requires one or more index columns")

    if SOMA_JOINID in schema.names:
        if schema.field(SOMA_JOINID).type != pa.int64():
            raise ValueError(f"{SOMA_JOINID} field must be of type Arrow int64")
    else:
        # add SOMA_JOINID
        schema = schema.append(pa.field(SOMA_JOINID, pa.int64()))

    # verify no illegal use of soma_ prefix
    for field_name in schema.names:
        if field_name.startswith("soma_") and field_name != SOMA_JOINID:
            raise ValueError(
                f"DataFrame schema may not contain fields with name prefix `soma_`: got `{field_name}`"
            )

    # verify that all index_column_names are present in the schema
    schema_names_set = set(schema.names)
    for index_column_name in index_column_names:
        if index_column_name.startswith("soma_") and index_column_name != SOMA_JOINID:
            raise ValueError(
                f'index_column_name other than "soma_joinid" must not begin with "soma_"; got "{index_column_name}"'
            )
        if index_column_name not in schema_names_set:
            schema_names_string = "{}".format(list(schema_names_set))
            raise ValueError(
                f"All index names must be dataframe schema: '{index_column_name}' not in {schema_names_string}"
            )
        if schema.field(index_column_name).type not in [
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
        ]:
            raise TypeError(
                f"Unsupported index type {schema.field(index_column_name).type}"
            )

    return schema


def _build_tiledb_schema(
    schema: pa.Schema,
    index_column_names: Sequence[str],
    tiledb_create_options: TileDBCreateOptions,
    context: SOMATileDBContext,
) -> tiledb.ArraySchema:
    """Converts an Arrow schema into a TileDB ArraySchema for creation."""
    dims = []
    for index_column_name in index_column_names:
        pa_type = schema.field(index_column_name).type
        dtype = arrow_types.tiledb_type_from_arrow_type(pa_type, is_indexed_column=True)
        domain: Tuple[Any, Any]
        if isinstance(dtype, str):
            domain = None, None
        elif np.issubdtype(dtype, NPInteger):
            iinfo = np.iinfo(cast(NPInteger, dtype))
            domain = iinfo.min, iinfo.max - 1
        elif np.issubdtype(dtype, NPFloating):
            finfo = np.finfo(cast(NPFloating, dtype))
            domain = finfo.min, finfo.max
        else:
            raise TypeError(f"Unsupported dtype {dtype}")

        # Default 2048 mods to 0 for 8-bit types and 0 is an invalid extent
        extent = tiledb_create_options.dim_tile(index_column_name)
        if isinstance(dtype, np.dtype) and dtype.itemsize == 1:
            extent = 64

        dim = tiledb.Dim(
            name=index_column_name,
            domain=domain,
            tile=extent,
            dtype=dtype,
            filters=tiledb_create_options.dim_filters(
                index_column_name,
                [
                    dict(
                        _type="ZstdFilter",
                        level=tiledb_create_options.dataframe_dim_zstd_level(),
                    )
                ],
            ),
        )
        dims.append(dim)

    dom = tiledb.Domain(dims, ctx=context.tiledb_ctx)

    attrs = []
    for attr_name in schema.names:
        if attr_name in index_column_names:
            continue
        attr = tiledb.Attr(
            name=attr_name,
            dtype=arrow_types.tiledb_type_from_arrow_type(schema.field(attr_name).type),
            filters=tiledb_create_options.attr_filters(attr_name, ["ZstdFilter"]),
            ctx=context.tiledb_ctx,
        )
        attrs.append(attr)

    cell_order, tile_order = tiledb_create_options.cell_tile_orders()

    return tiledb.ArraySchema(
        domain=dom,
        attrs=attrs,
        sparse=True,
        allows_duplicates=tiledb_create_options.get("allows_duplicates", False),
        offsets_filters=tiledb_create_options.offsets_filters(),
        validity_filters=tiledb_create_options.validity_filters(),
        capacity=tiledb_create_options.get("capacity", 100000),
        cell_order=cell_order,
        # As of TileDB core 2.8.2, we cannot consolidate string-indexed sparse arrays with
        # col-major tile order: so we write ``X`` with row-major tile order.
        tile_order=tile_order,
        ctx=context.tiledb_ctx,
    )
