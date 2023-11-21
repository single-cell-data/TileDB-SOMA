# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Implementation of a SOMA DataFrame
"""
from typing import Any, Dict, Optional, Sequence, Tuple, Type, Union, cast

import numpy as np
import pandas as pd
import pyarrow as pa
import somacore
import tiledb
from somacore import options
from typing_extensions import Self

from . import _arrow_types, _util
from . import pytiledbsoma as clib
from ._constants import SOMA_JOINID
from ._query_condition import QueryCondition
from ._read_iters import TableReadIter
from ._tiledb_array import TileDBArray
from ._types import NPFloating, NPInteger, OpenTimestamp, Slice, is_slice_of
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context
from .options._tiledb_create_options import TileDBCreateOptions

_UNBATCHED = options.BatchSize()


class DataFrame(TileDBArray, somacore.DataFrame):
    """:class:`DataFrame` is a multi-column table with a user-defined schema. The
    schema is expressed as an
    `Arrow Schema <https://arrow.apache.org/docs/python/generated/pyarrow.Schema.html>`_,
    and defines the column names and value types.

    Every :class:`DataFrame` must contain a column called ``soma_joinid``, of type
    ``int64``, with negative values explicitly disallowed. The ``soma_joinid``
    column contains a unique value for each row in the dataframe, and in some
    cases (e.g., as part of an :class:`Experiment`), acts as a join key for other
    objects, such as :class:`SparseNDArray`.

    Lifecycle:
        Experimental.

    Examples:
        >>> import pyarrow as pa
        >>> import tiledbsoma
        >>> schema = pa.schema(
        ...     [
        ...         ("soma_joinid", pa.int64()),
        ...         ("A", pa.float32()),
        ...         ("B", pa.large_string()),
        ...     ]
        ... )
        >>> with tiledbsoma.DataFrame.create("./test_dataframe", schema=schema) as df:
        ...     data = pa.Table.from_pydict(
        ...         {
        ...             "soma_joinid": [0, 1, 2],
        ...             "A": [1.0, 2.7182, 3.1214],
        ...             "B": ["one", "e", "pi"],
        ...         }
        ...     )
        ...     df.write(data)
        >>> with tiledbsoma.DataFrame.open("./test_dataframe") as df:
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


        >>> import pyarrow as pa
        >>> import tiledbsoma
        >>> schema = pa.schema(
        ...    [
        ...        ("soma_joinid", pa.int64()),
        ...        ("A", pa.float32()),
        ...        ("B", pa.large_string()),
        ...    ]
        ...)
        >>> with tiledbsoma.DataFrame.create(
        ...     "./test_dataframe_2",
        ...     schema=schema,
        ...     index_column_names=["A", "B"],
        ...     domain=[(0.0, 10.0), None],
        ... ) as df:
        ...     data = pa.Table.from_pydict(
        ...         {
        ...             "soma_joinid": [0, 1, 2],
        ...             "A": [1.0, 2.7182, 3.1214],
        ...             "B": ["one", "e", "pi"],
        ...         }
        ...     )
        ...     df.write(data)
        >>> with tiledbsoma.DataFrame.open("./test_dataframe_2") as df:
        ...     print(df.schema)
        ...     print("---")
        ...     print(df.read().concat().to_pandas())
        soma_joinid: int64
        ---
                A    B  soma_joinid
        0  1.0000  one            0
        1  2.7182    e            1
        2  3.1214   pi            2

        Here the index-column names are specified. The domain is entirely optional: if
        it's omitted, defaults will be applied yielding the largest possible domain for
        each index column's datatype.  If the domain is specified, it must be a
        tuple/list of equal length to ``index_column_names``. It can be ``None`` in
        a given slot, meaning use the largest possible domain. For string/bytes types,
        it must be ``None``.
    """

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID,),
        domain: Optional[Sequence[Optional[Tuple[Any, Any]]]] = None,
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
    ) -> "DataFrame":
        """Creates the data structure on disk/S3/cloud.

        Args:
            schema:
                `Arrow schema <https://arrow.apache.org/docs/python/generated/pyarrow.Schema.html>`_
                defining the per-column schema. This schema must define all columns, including
                columns to be named as index columns.  If the schema includes types unsupported by
                the SOMA implementation, an error will be raised.
            index_column_names:
                A list of column names to use as user-defined
                index columns (e.g., ``['cell_type', 'tissue_type']``).
                All named columns must exist in the schema, and at least one
                index column name is required.
            domain:
                An optional sequence of tuples specifying the domain of each
                index column. Each tuple should be a pair consisting of the minimum and
                maximum values storable in the index column. For example, if there is a
                single int64-valued index column, then ``domain`` might be ``[(100,
                200)]`` to indicate that values between 100 and 200, inclusive, can be
                stored in that column.  If provided, this sequence must have the same
                length as ``index_column_names``, and the index-column domain will be as
                specified.  If omitted entirely, or if ``None`` in a given dimension,
                the corresponding index-column domain will use the minimum and maximum
                possible values for the column's datatype.  This makes a
                :class:`DataFrame` growable.
            platform_config:
                Platform-specific options used to create this array.
                This may be provided as settings in a dictionary, with options
                located in the ``{'tiledb': {'create': ...}}`` key,
                or as a :class:`~tiledbsoma.TileDBCreateOptions` object.
            tiledb_timestamp:
                If specified, overrides the default timestamp
                used to open this object. If unset, uses the timestamp provided by
                the context.
            enumeration:
                If specified, enumerate attributes with the given sequence of values.


        Returns:
            The DataFrame.

        Raises:
            TypeError:
                If the ``schema`` parameter specifies an unsupported type,
                or if ``index_column_names`` specifies a non-indexable column.
            ValueError:
                If the ``index_column_names`` is malformed or specifies
                an undefined column name.
            ValueError:
                If the ``schema`` specifies illegal column names.
            TileDBError:
                If unable to create the underlying object.

        Examples:
            >>> df = pd.DataFrame(data={"soma_joinid": [0, 1], "col1": ["a", "b"]})
            ... with tiledbsoma.DataFrame.create(
            ...    "a_dataframe", schema=pa.Schema.from_pandas(df)
            ... ) as soma_df:
            ...     soma_df.write(pa.Table.from_pandas(df, preserve_index=False))
            ...
            >>> with tiledbsoma.open("a_dataframe") as soma_df:
            ...     a_df = soma_df.read().concat().to_pandas()
            ...
            >>> a_df
               soma_joinid col1
            0            0    a
            1            1    b

        Lifecycle:
            Experimental.
        """
        context = _validate_soma_tiledb_context(context)
        schema = _canonicalize_schema(schema, index_column_names)
        tdb_schema = _build_tiledb_schema(
            schema,
            index_column_names,
            domain,
            TileDBCreateOptions.from_platform_config(platform_config),
            context,
        )
        handle = cls._create_internal(uri, tdb_schema, context, tiledb_timestamp)
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    def keys(self) -> Tuple[str, ...]:
        """Returns the names of the columns when read back as a dataframe.

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
    def domain(self) -> Tuple[Tuple[Any, Any], ...]:
        """Returns a tuple of minimum and maximum values, inclusive, storable
        on each index column of the dataframe.

        Lifecycle:
            Experimental.
        """
        return self._tiledb_domain()

    @property
    def count(self) -> int:
        """Returns the number of rows in the dataframe. Same as ``len(df)``.

        Lifecycle:
            Experimental.
        """
        self._check_open_read()
        return cast(int, self._soma_reader().nnz())

    def enumeration(self, name: str) -> Tuple[Any, ...]:
        """Doc place holder.

        Returns:
            Tuple[Any, ...]: _description_
        """
        return tuple(self._soma_reader().get_enum(name))

    def column_to_enumeration(self, name: str) -> str:
        return str(self._soma_reader().get_enum_label_on_attr(name))

    def __len__(self) -> int:
        """Returns the number of rows in the dataframe. Same as ``df.count``."""
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
        """Reads a user-defined subset of data, addressed by the dataframe indexing columns,
        optionally filtered, and return results as one or more `Arrow tables <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_.

        Args:
            coords:
                For each index dimension, which rows to read.
                Defaults to ``None``, meaning no constraint -- all IDs.
            column_names:
                The named columns to read and return.
                Defaults to ``None``, meaning no constraint -- all column names.
            result_order:
                Order of read results.
                This can be one of 'row-major', 'col-major', or 'auto'.
            value_filter:
                An optional [value filter] to apply to the results.
                Defaults to no filter.
            partitions:
                An optional :class:`ReadPartitions` hint to indicate
                how results should be organized.

        Returns:
            A :class:`TableReadIter` that can be used to iterate through the result set.

        Raises:
            SOMAError:
                If ``value_filter`` can not be parsed.
            ValueError:
                If ``coords`` are malformed or do not index this DataFrame.
            SOMAError:
                If the object is not open for reading.

        Notes:
            The ``coords`` parameter will support, per dimension:
            a list of values of the type of the indexed column.

            Acceptable ways to index:

            * A sequence of coordinates is accepted, one per dimension.
            * Sequence length must be <= number of dimensions.
            * If the sequence contains missing coordinates (length less than number of dimensions),
              then ``slice(None)`` -- i.e. no constraint -- is assumed for the missing dimensions.
            * Per-dimension, explicitly specified coordinates can be one of: None, a value, a
              list/``numpy.ndarray``/``pyarrow.Array``/etc of values, a slice, etc.
            * Slices are doubly inclusive: ``slice(2,4)`` means [2,3,4] not [2,3].
              Slice steps are not supported.
              Slices can be ``slice(None)``, meaning select all in that dimension, and may be half-specified,
              e.g.  ``slice(2,None)`` or ``slice(None,4)``.
            * Negative indexing is unsupported.

        Lifecycle:
            Experimental.
        """
        del batch_size, platform_config  # Currently unused.
        _util.check_unpartitioned(partitions)
        self._check_open_read()

        schema = self._handle.schema
        query_condition = None
        if value_filter is not None:
            query_condition = QueryCondition(value_filter)

        sr = self._soma_reader(
            schema=schema,  # query_condition needs this
            column_names=column_names,
            query_condition=query_condition,
            result_order=result_order,
        )

        self._set_reader_coords(sr, coords)

        # TODO: platform_config
        # TODO: batch_size

        return TableReadIter(sr)

    def write(
        self, values: pa.Table, platform_config: Optional[options.PlatformConfig] = None
    ) -> Self:
        """Writes an `Arrow table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
        to the persistent object. As duplicate index values are not allowed, index values already
        present in the object are overwritten
        and new index values are added.

        Args:
            values:
                An `Arrow table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
                containing all columns, including the index columns. The schema for the values must
                match the schema for the :class:`DataFrame`.

                If a column is of categorical type in the schema and a flattened/non-categorical
                column is presented for data on write, a ``ValueError`` is raised.  If a column is
                of non-categorical type in the schema and a categorical column is presented for data
                on write, the data are written as an array of category values, and the category-type
                information is not saved.

        Raises:
            TypeError:
                If the ``values`` parameter is an unsupported type.
            ValueError:
                If the ``values`` parameter is an empty table.
            SOMAError:
                If the object is not open for writing.

        Lifecycle:
            Experimental.
        """
        _util.check_type("values", values, (pa.Table,))

        dim_cols_map: Dict[str, pd.DataFrame] = {}
        attr_cols_map: Dict[str, pd.DataFrame] = {}
        dim_names_set = self.index_column_names
        n = None

        for col_info in values.schema:
            name = col_info.name
            col = values.column(name)
            n = len(col)

            if self._handle.schema.has_attr(name):
                attr = self._handle.schema.attr(name)

                # Add the enumeration values to the TileDB Array from ArrowArray
                if attr.enum_label is not None and col.num_chunks != 0:
                    if not pa.types.is_dictionary(col_info.type):
                        raise ValueError(
                            "Expected dictionary type for enumerated attribute "
                            f"{name} but saw {col_info.type}"
                        )

                    enmr = self._handle.enum(attr.name)

                    # get new enumeration values, maintain original ordering
                    update_vals = []
                    for new_val in col.chunk(0).dictionary.tolist():
                        if new_val not in enmr.values():
                            update_vals.append(new_val)

                    # only extend if there are new values
                    if update_vals:
                        se = tiledb.ArraySchemaEvolution(self.context.tiledb_ctx)
                        if np.issubdtype(enmr.dtype.type, np.str_):
                            extend_vals = np.array(update_vals, "U")
                        elif np.issubdtype(enmr.dtype.type, np.bytes_):
                            extend_vals = np.array(update_vals, "S")
                        else:
                            extend_vals = np.array(update_vals, enmr.dtype)
                        new_enmr = enmr.extend(extend_vals)
                        se.extend_enumeration(new_enmr)
                        se.array_evolve(uri=self.uri)

            cols_map = dim_cols_map if name in dim_names_set else attr_cols_map
            if pa.types.is_dictionary(col.type):
                if col.num_chunks != 0:
                    if name in dim_names_set:
                        # Dims are never categorical. Decategoricalize for them.
                        cols_map[name] = pa.chunked_array(
                            [chunk.dictionary_decode() for chunk in col.chunks]
                        )
                    else:
                        attr = self._handle.schema.attr(name)
                        if attr.enum_label is not None:
                            # Normal case: writing categorical data to categorical schema.
                            cols_map[name] = col.chunk(0).indices.to_pandas()
                        else:
                            # Schema is non-categorical but the user is writing categorical.
                            # Simply decategoricalize for them.
                            cols_map[name] = pa.chunked_array(
                                [chunk.dictionary_decode() for chunk in col.chunks]
                            )
                else:
                    cols_map[name] = col.to_pandas()

            else:
                if name not in dim_names_set:
                    attr = self._handle.schema.attr(name)
                    if attr.enum_label is not None:
                        raise ValueError(
                            f"Categorical column {name} must be presented with categorical data"
                        )

                cols_map[name] = col.to_pandas()

        if n is None:
            raise ValueError(f"did not find any column names in {values.schema.names}")

        # We need to produce the dim cols in the same order as they're present in the TileDB schema
        # (tracked by self.index_column_names). This is important in the multi-index case.  Suppose
        # the Arrow schema has two index columns in the order "burger" and "meister", and suppose
        # the user set index_column_names = ["meister", "burger"] when creating the TileDB schema.
        # Then the above for-loop over the Arrow schema will find the former ordering, but for the
        # ``writer[dims] = attrs`` below we must have dims with the latter ordering.
        dim_cols_list = [dim_cols_map[name] for name in self.index_column_names]
        dim_cols_tuple = tuple(dim_cols_list)
        self._handle.writer[dim_cols_tuple] = attr_cols_map
        tiledb_create_options = TileDBCreateOptions.from_platform_config(
            platform_config
        )
        if tiledb_create_options.consolidate_and_vacuum:
            self._consolidate_and_vacuum()

        return self

    def _set_reader_coord(
        self, sr: clib.SOMAArray, dim_idx: int, dim: tiledb.Dim, coord: object
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
            # A ``None`` or empty start is always equivalent to empty str/bytes.
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
        sr: clib.SOMAArray,
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
        self, sr: clib.SOMAArray, dim_idx: int, dim: tiledb.Dim, coord: Slice[Any]
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
    """Turns an Arrow schema into the canonical version and checks for errors.

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
                f"DataFrame schema may not contain fields with name prefix ``soma_``: got ``{field_name}``"
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


def _build_tiledb_schema(
    schema: pa.Schema,
    index_column_names: Sequence[str],
    domain: Optional[Sequence[Optional[Tuple[Any, Any]]]],
    tiledb_create_options: TileDBCreateOptions,
    context: SOMATileDBContext,
) -> tiledb.ArraySchema:
    """Converts an Arrow schema into a TileDB ArraySchema for creation."""

    if domain is None:
        domain = tuple(None for _ in index_column_names)
    else:
        ndom = len(domain)
        nidx = len(index_column_names)
        if ndom != nidx:
            raise ValueError(
                f"if domain is specified, it must have the same length as index_column_names; got {ndom} != {nidx}"
            )

    dims = []
    for index_column_name, slot_domain in zip(index_column_names, domain):
        pa_type = schema.field(index_column_name).type
        dtype = _arrow_types.tiledb_type_from_arrow_type(
            pa_type, is_indexed_column=True
        )

        slot_domain = _fill_out_slot_domain(
            slot_domain, index_column_name, pa_type, dtype
        )

        extent = _find_extent_for_domain(
            index_column_name, tiledb_create_options, dtype, slot_domain
        )

        dim = tiledb.Dim(
            name=index_column_name,
            domain=slot_domain,
            tile=extent,
            dtype=dtype,
            filters=tiledb_create_options.dim_filters_tiledb(
                index_column_name,
                [
                    dict(
                        _type="ZstdFilter",
                        level=tiledb_create_options.dataframe_dim_zstd_level,
                    )
                ],
            ),
        )
        dims.append(dim)

    dom = tiledb.Domain(dims, ctx=context.tiledb_ctx)

    attrs = []
    enums = []
    metadata = schema.metadata or {}
    for pa_attr in schema:
        attr_name = pa_attr.name

        if attr_name in index_column_names:
            continue

        has_enum = pa.types.is_dictionary(pa_attr.type)

        if has_enum:
            enmr_dtype: np.dtype[Any]
            vtype = pa_attr.type.value_type
            if pa.types.is_large_string(vtype) or pa.types.is_string(vtype):
                enmr_dtype = np.dtype("U")
            elif pa.types.is_large_binary(vtype) or pa.types.is_binary(vtype):
                enmr_dtype = np.dtype("S")
            else:
                enmr_dtype = np.dtype(vtype.to_pandas_dtype())
            enums.append(
                tiledb.Enumeration(
                    name=attr_name,
                    ordered=pa_attr.type.ordered,
                    dtype=enmr_dtype,
                )
            )

        attr = tiledb.Attr(
            name=attr_name,
            dtype=_arrow_types.tiledb_type_from_arrow_type(
                schema.field(attr_name).type
            ),
            nullable=metadata.get(attr_name.encode("utf-8")) == b"nullable",
            filters=tiledb_create_options.attr_filters_tiledb(
                attr_name, ["ZstdFilter"]
            ),
            enum_label=attr_name if has_enum else None,
            ctx=context.tiledb_ctx,
        )
        attrs.append(attr)

    cell_order, tile_order = tiledb_create_options.cell_tile_orders()

    return tiledb.ArraySchema(
        domain=dom,
        attrs=attrs,
        enums=enums,
        sparse=True,
        allows_duplicates=tiledb_create_options.allows_duplicates,
        offsets_filters=tiledb_create_options.offsets_filters_tiledb(),
        validity_filters=tiledb_create_options.validity_filters_tiledb(),
        capacity=tiledb_create_options.capacity,
        cell_order=cell_order,
        # As of TileDB core 2.8.2, we cannot consolidate string-indexed sparse arrays with
        # col-major tile order: so we write ``X`` with row-major tile order.
        tile_order=tile_order,
        ctx=context.tiledb_ctx,
    )


def _fill_out_slot_domain(
    slot_domain: Optional[Tuple[Any, Any]],
    index_column_name: str,
    pa_type: pa.DataType,
    dtype: Any,
) -> Tuple[Any, Any]:
    """Helper function for _build_tiledb_schema. Given a user-specified domain for a
    dimension slot -- which may be ``None``, or a two-tuple of which either element
    may be ``None`` -- return either what the user specified (if adequate) or
    sensible type-inferred values appropriate to the datatype.
    """
    if slot_domain is not None:
        # User-specified; go with it when possible
        if (
            pa_type == pa.string()
            or pa_type == pa.large_string()
            or pa_type == pa.binary()
            or pa_type == pa.large_binary()
        ):
            # TileDB Embedded won't raise an error if the user asks for, say
            # domain=[("a", "z")].  But it will simply _ignore_ the request and
            # use [("", "")]. The decision here is to explicitly reject an
            # unsupported operation.
            raise ValueError(
                "TileDB str and bytes index-column types do not support domain specfication"
            )
        if index_column_name == SOMA_JOINID:
            lo = slot_domain[0]
            hi = slot_domain[1]
            if lo is not None and lo < 0:
                raise ValueError(
                    f"soma_joinid indices cannot be negative; got lower bound {lo}"
                )
            if hi is not None and hi < 0:
                raise ValueError(
                    f"soma_joinid indices cannot be negative; got upper bound {hi}"
                )

    elif isinstance(dtype, str):
        slot_domain = None, None
    elif np.issubdtype(dtype, NPInteger):
        iinfo = np.iinfo(cast(NPInteger, dtype))
        slot_domain = iinfo.min, iinfo.max - 1
        # Here the slot_domain isn't specified by the user; we're setting it.
        # The SOMA spec disallows negative soma_joinid.
        if index_column_name == SOMA_JOINID:
            slot_domain = (0, 2**31 - 2)  # R-friendly, which 2**63-1 is not
    elif np.issubdtype(dtype, NPFloating):
        finfo = np.finfo(cast(NPFloating, dtype))
        slot_domain = finfo.min, finfo.max

    # The `iinfo.min+1` is necessary as of tiledb core 2.15 / tiledb-py 0.21.1 since
    # `iinfo.min` maps to `NaT` (not a time), resulting in
    #   TypeError: invalid domain extent, domain cannot be safely cast to dtype dtype('<M8[s]')
    #
    # The `iinfo.max-delta` is necessary since with iinfo.min being bumped by 1, without subtracting
    # we would get
    #   tiledb.cc.TileDBError: [TileDB::Dimension] Error: Tile extent check failed; domain max
    #   expanded to multiple of tile extent exceeds max value representable by domain type. Reduce
    #   domain max by 1 tile extent to allow for expansion.
    elif dtype == "datetime64[s]":
        iinfo = np.iinfo(cast(NPInteger, np.int64))
        slot_domain = np.datetime64(iinfo.min + 1, "s"), np.datetime64(
            iinfo.max - 1000000, "s"
        )
    elif dtype == "datetime64[ms]":
        iinfo = np.iinfo(cast(NPInteger, np.int64))
        slot_domain = np.datetime64(iinfo.min + 1, "ms"), np.datetime64(
            iinfo.max - 1000000, "ms"
        )
    elif dtype == "datetime64[us]":
        iinfo = np.iinfo(cast(NPInteger, np.int64))
        slot_domain = np.datetime64(iinfo.min + 1, "us"), np.datetime64(
            iinfo.max - 1000000, "us"
        )
    elif dtype == "datetime64[ns]":
        iinfo = np.iinfo(cast(NPInteger, np.int64))
        slot_domain = np.datetime64(iinfo.min + 1, "ns"), np.datetime64(
            iinfo.max - 1000000, "ns"
        )

    else:
        raise TypeError(f"Unsupported dtype {dtype}")

    return slot_domain


def _find_extent_for_domain(
    index_column_name: str,
    tiledb_create_options: TileDBCreateOptions,
    dtype: Any,
    slot_domain: Tuple[Any, Any],
) -> Any:
    """Helper function for _build_tiledb_schema. Returns a tile extent that is
    small enough for the index-column type, and that also fits within the
    user-specified slot domain (if any).
    """

    # Default 2048 mods to 0 for 8-bit types and 0 is an invalid extent
    extent = tiledb_create_options.dim_tile(index_column_name)
    if isinstance(dtype, np.dtype) and dtype.itemsize == 1:
        extent = 64

    if isinstance(dtype, str):
        return extent

    lo, hi = slot_domain
    if lo is None or hi is None:
        return extent

    if np.issubdtype(dtype, NPInteger) or np.issubdtype(dtype, NPFloating):
        return min(extent, hi - lo + 1)

    if dtype == "datetime64[s]":
        ilo = int(lo.astype("int64"))
        ihi = int(hi.astype("int64"))
        iextent = min(extent, ihi - ilo + 1)
        return np.datetime64(iextent, "s")

    if dtype == "datetime64[ms]":
        ilo = int(lo.astype("int64"))
        ihi = int(hi.astype("int64"))
        iextent = min(extent, ihi - ilo + 1)
        return np.datetime64(iextent, "ms")

    if dtype == "datetime64[us]":
        ilo = int(lo.astype("int64"))
        ihi = int(hi.astype("int64"))
        iextent = min(extent, ihi - ilo + 1)
        return np.datetime64(iextent, "us")

    if dtype == "datetime64[ns]":
        ilo = int(lo.astype("int64"))
        ihi = int(hi.astype("int64"))
        iextent = min(extent, ihi - ilo + 1)
        return np.datetime64(iextent, "ns")

    return extent
