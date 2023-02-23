from typing import Any, Mapping, Optional, Sequence, Tuple, Type, Union, cast

import numpy as np
import pyarrow as pa
import somacore
import tiledb
from somacore import options
from typing_extensions import Self

from . import _arrow_types, _util
from . import libtiledbsoma as clib
from ._constants import SOMA_JOINID
from ._query_condition import QueryCondition
from ._read_iters import TableReadIter
from ._tiledb_array import TileDBArray
from ._types import NPFloating, NPInteger, Slice, is_slice_of
from .options import SOMATileDBContext
from .options._tiledb_create_options import TileDBCreateOptions

_UNBATCHED = options.BatchSize()


class DataFrame(TileDBArray, somacore.DataFrame):
    """
    ``DataFrame`` is a multi-column table with a user-defined schema. The
    schema is expressed as an Arrow Schema, and defines the column names
    and value types.

    Every ``DataFrame`` must contain a column called ``soma_joinid``, of type
    ``int64``, with negative values explicitly disallowed. The ``soma_joinid``
    column contains a unique value for each row in the dataframe, and in some
    cases (e.g., as part of an ``Experiment``), acts as a join key for other
    objects, such as ``SparseNDArray``.

    [lifecycle: experimental]

    Examples:
    ---------
    >>> import pyarrow as pa
    >>> import tiledbsoma
    >>> schema = pa.schema(
    ...     [
    ...         ("soma_joinid", pa.int64()),
    ...         ("A", pa.float32()),
    ...         ("B", pa.large_string()),
    ...     ]
    ... )
    ... with tiledbsoma.DataFrame.create("./test_dataframe", schema=schema) as df:
    ...     data = pa.Table.from_pydict(
    ...         {
    ...             "soma_joinid": [0, 1, 2],
    ...             "A": [1.0, 2.7182, 3.1214],
    ...             "B": ["one", "e", "pi"],
    ...         }
    ...     )
    ...     df.write(data)
    ... with tiledbsoma.DataFrame.open("./test_dataframe") as df:
    ...     print(df.schema)
    ...     print("---")
    ...     print(df.read().concat().to_pandas())
    ...
    soma_joinid: int64
    A: float
    B: large_string
    ---
       soma_joinid       A    B
    0            0  1.0000  one
    1            1  2.7182    e
    2            2  3.1214   pi
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
        _util.check_unpartitioned(partitions)
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

        self._set_reader_coords(sr, coords)

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
        _util.check_type("values", values, (pa.Table,))

        del platform_config  # unused
        dim_cols_map = {}
        attr_cols_map = {}
        dim_names_set = self.index_column_names
        n = None

        for name in values.schema.names:
            n = len(values.column(name))
            if name in dim_names_set:
                dim_cols_map[name] = values.column(name).to_pandas()
            else:
                attr_cols_map[name] = values.column(name).to_pandas()
        if n is None:
            raise ValueError(f"did not find any column names in {values.schema.names}")

        # We need to produce the dim cols in the same order as they're present in the TileDB schema
        # (tracked by self.index_column_names). This is important in the multi-index case.  Suppose
        # the Arrow schema has two index columns in the order "burger" and "meister", and suppose
        # the user set index_column_names = ["meister", "burger"] when creating the TileDB schema.
        # Then the above for-loop over the Arrow schema will find the former ordering, but for the
        # `writer[dims] = attrs` below we must have dims with the latter ordering.
        dim_cols_list = [dim_cols_map[name] for name in self.index_column_names]
        dim_cols_tuple = tuple(dim_cols_list)
        self._handle.writer[dim_cols_tuple] = attr_cols_map

        return self

    def _set_reader_coord(
        self, sr: clib.SOMAReader, dim_idx: int, dim: tiledb.Dim, coord: object
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

        # Note: slice(None, None) matches the is_slice_of part, unless we also check the dim-type
        # part.
        if (is_slice_of(coord, str) or is_slice_of(coord, bytes)) and (
            dim.dtype == "str" or dim.dtype == "bytes"
        ):
            _util.validate_slice(coord)
            # Figure out which one.
            dim_type: Union[Type[str], Type[bytes]] = type(dim.domain[0])
            # A `None` or empty start is always equivalent to empty str/bytes.
            start = coord.start or dim_type()
            if coord.stop is None:
                # There's no way to specify "to infinity" for strings.
                # We have to get the nonempty domain and use that as the end.
                _, stop = self._handle.reader.nonempty_domain()[dim_idx]
            else:
                stop = coord.stop
            sr.set_dim_ranges_string_or_bytes(dim.name, [(start, stop)])
            return True

        # Note: slice(None, None) matches the is_slice_of part, unless we also check the dim-type
        # part.
        if is_slice_of(coord, np.datetime64) and dim.dtype.name.startswith(
            "datetime64"
        ):
            _util.validate_slice(coord)
            # These timestamp types are stored in Arrow as well as TileDB as 64-bit integers (with
            # distinguishing metadata of course). For purposes of the query logic they're just
            # int64.
            istart = coord.start or dim.domain[0]
            istart = int(istart.astype("int64"))
            istop = coord.stop or dim.domain[1]
            istop = int(istop.astype("int64"))
            sr.set_dim_ranges_int64(dim.name, [(istart, istop)])
            return True

        if super()._set_reader_coord(sr, dim_idx, dim, coord):
            return True

        return False

    def _set_reader_coord_by_py_seq_or_np_array(
        self,
        sr: clib.SOMAReader,
        dim_idx: int,
        dim: tiledb.Dim,
        coord: object,
    ) -> bool:
        if isinstance(coord, np.ndarray):
            if coord.ndim != 1:
                raise ValueError(
                    f"only 1D numpy arrays may be used to index; got {coord.ndim}"
                )

        # See libtiledbsoma.cc for more context on why we need the
        # explicit type-check here.

        if dim.dtype == np.int64:
            sr.set_dim_points_int64(dim.name, coord)
        elif dim.dtype == np.int32:
            sr.set_dim_points_int32(dim.name, coord)
        elif dim.dtype == np.int16:
            sr.set_dim_points_int16(dim.name, coord)
        elif dim.dtype == np.int8:
            sr.set_dim_points_int8(dim.name, coord)

        elif dim.dtype == np.uint64:
            sr.set_dim_points_uint64(dim.name, coord)
        elif dim.dtype == np.uint32:
            sr.set_dim_points_uint32(dim.name, coord)
        elif dim.dtype == np.uint16:
            sr.set_dim_points_uint16(dim.name, coord)
        elif dim.dtype == np.uint8:
            sr.set_dim_points_uint8(dim.name, coord)

        elif dim.dtype == np.float64:
            sr.set_dim_points_float64(dim.name, coord)
        elif dim.dtype == np.float32:
            sr.set_dim_points_float32(dim.name, coord)

        elif dim.dtype == "str" or dim.dtype == "bytes":
            sr.set_dim_points_string_or_bytes(dim.name, coord)

        elif (
            dim.dtype == "datetime64[s]"
            or dim.dtype == "datetime64[ms]"
            or dim.dtype == "datetime64[us]"
            or dim.dtype == "datetime64[ns]"
        ):
            if not isinstance(coord, (tuple, list, np.ndarray)):
                raise ValueError(
                    f"unhandled coord type {type(coord)} for index column named {dim.name}"
                )
            icoord = [
                int(e.astype("int64")) if isinstance(e, np.datetime64) else e
                for e in coord
            ]
            sr.set_dim_points_int64(dim.name, icoord)

        # TODO: bool

        else:
            raise ValueError(
                f"unhandled type {dim.dtype} for index column named {dim.name}"
            )

        return True

    def _set_reader_coord_by_numeric_slice(
        self, sr: clib.SOMAReader, dim_idx: int, dim: tiledb.Dim, coord: Slice[Any]
    ) -> bool:

        try:
            lo_hi = _util.slice_to_numeric_range(coord, dim.domain)
        except _util.NonNumericDimensionError:
            return False  # We only handle numeric dimensions here.

        if not lo_hi:
            return True

        elif dim.dtype == np.int64:
            sr.set_dim_ranges_int64(dim.name, [lo_hi])
            return True
        elif dim.dtype == np.int32:
            sr.set_dim_ranges_int32(dim.name, [lo_hi])
            return True
        elif dim.dtype == np.int16:
            sr.set_dim_ranges_int16(dim.name, [lo_hi])
            return True
        elif dim.dtype == np.int8:
            sr.set_dim_ranges_int8(dim.name, [lo_hi])
            return True

        elif dim.dtype == np.uint64:
            sr.set_dim_ranges_uint64(dim.name, [lo_hi])
            return True
        elif dim.dtype == np.uint32:
            sr.set_dim_ranges_uint32(dim.name, [lo_hi])
            return True
        elif dim.dtype == np.uint16:
            sr.set_dim_ranges_uint16(dim.name, [lo_hi])
            return True
        elif dim.dtype == np.uint8:
            sr.set_dim_ranges_uint8(dim.name, [lo_hi])
            return True

        elif dim.dtype == np.float64:
            sr.set_dim_ranges_float64(dim.name, [lo_hi])
            return True
        elif dim.dtype == np.float32:
            sr.set_dim_ranges_float32(dim.name, [lo_hi])
            return True

        # TODO:
        # elif dim.dtype == np.bool_:

        return False


def _canonicalize_schema(
    schema: pa.Schema, index_column_names: Sequence[str]
) -> pa.Schema:
    """Turn a Arrow schema into the canonical version and check for errors.

    Returns a schema, which may be modified by the addition of required columns
    (e.g. ``soma_joinid``).
    """
    _util.check_type("schema", schema, (pa.Schema,))
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
            pa.timestamp("s"),
            pa.timestamp("ms"),
            pa.timestamp("us"),
            pa.timestamp("ns"),
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
        dtype = _arrow_types.tiledb_type_from_arrow_type(
            pa_type, is_indexed_column=True
        )
        domain: Tuple[Any, Any]
        if isinstance(dtype, str):
            domain = None, None
        elif np.issubdtype(dtype, NPInteger):
            iinfo = np.iinfo(cast(NPInteger, dtype))
            domain = iinfo.min, iinfo.max - 1
        elif np.issubdtype(dtype, NPFloating):
            finfo = np.finfo(cast(NPFloating, dtype))
            domain = finfo.min, finfo.max

        elif dtype == "datetime64[s]":
            iinfo = np.iinfo(cast(NPInteger, np.int64))
            domain = np.datetime64(iinfo.min, "s"), np.datetime64(iinfo.max - 1, "s")
        elif dtype == "datetime64[ms]":
            iinfo = np.iinfo(cast(NPInteger, np.int64))
            domain = np.datetime64(iinfo.min, "ms"), np.datetime64(iinfo.max - 1, "ms")
        elif dtype == "datetime64[us]":
            iinfo = np.iinfo(cast(NPInteger, np.int64))
            domain = np.datetime64(iinfo.min, "us"), np.datetime64(iinfo.max - 1, "us")
        elif dtype == "datetime64[ns]":
            iinfo = np.iinfo(cast(NPInteger, np.int64))
            domain = np.datetime64(iinfo.min, "ns"), np.datetime64(iinfo.max - 1, "ns")

        else:
            raise TypeError(f"Unsupported dtype {dtype}")
        if index_column_name == SOMA_JOINID:
            domain = (0, domain[1])

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
            dtype=_arrow_types.tiledb_type_from_arrow_type(
                schema.field(attr_name).type
            ),
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
