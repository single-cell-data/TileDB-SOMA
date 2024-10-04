# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""
Implementation of a SOMA DataFrame
"""
from typing import Any, List, Optional, Sequence, Tuple, Type, Union, cast

import numpy as np
import pyarrow as pa
import somacore
from somacore import options
from typing_extensions import Self

from tiledbsoma._flags import NEW_SHAPE_FEATURE_FLAG_ENABLED

from . import _arrow_types, _util
from . import pytiledbsoma as clib
from ._constants import SOMA_JOINID
from ._exception import SOMAError, map_exception_for_create
from ._query_condition import QueryCondition
from ._read_iters import TableReadIter
from ._soma_array import SOMAArray
from ._tdb_handles import DataFrameWrapper
from ._types import NPFloating, NPInteger, OpenTimestamp, Slice, is_slice_of
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context
from .options._tiledb_create_write_options import (
    TileDBCreateOptions,
    TileDBWriteOptions,
)

_UNBATCHED = options.BatchSize()
AxisDomain = Union[None, Tuple[Any, Any], List[Any]]
Domain = Sequence[AxisDomain]


class DataFrame(SOMAArray, somacore.DataFrame):
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
        Maturing.

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

    _wrapper_type = DataFrameWrapper

    @classmethod
    def create(
        cls,
        uri: str,
        *,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID,),
        domain: Optional[Domain] = None,
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
            tiledbsoma.AlreadyExistsError:
                If the underlying object already exists at the given URI.
            tiledbsoma.NotCreateableError:
                If the URI is malformed for a particular storage backend.
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
            Maturing.
        """
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
            clib.SOMADataFrame.create(
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
            Maturing.
        """
        return self._tiledb_array_keys()

    @property
    def index_column_names(self) -> Tuple[str, ...]:
        """Returns index (dimension) column names.

        Lifecycle:
            Maturing.
        """
        return self._tiledb_dim_names()

    @property
    def domain(self) -> Tuple[Tuple[Any, Any], ...]:
        """Returns tuples of minimum and maximum values, one tuple per index column, currently storable
        on each index column of the dataframe. These can be resized up to ``maxdomain``.

        Lifecycle:
            Maturing.
        """
        return self._domain()

    @property
    def maxdomain(self) -> Tuple[Tuple[Any, Any], ...]:
        """Returns tuples of minimum and maximum values, one tuple per index column, to which the dataframe
        can have its domain resized.

        Lifecycle:
            Maturing.
        """
        return self._maxdomain()

    @property
    def count(self) -> int:
        """Returns the number of rows in the dataframe. Same as ``len(df)``.

        Lifecycle:
            Maturing.
        """
        self._check_open_read()
        # if is it in read open mode, then it is a DataFrameWrapper
        return cast(DataFrameWrapper, self._handle).count

    @property
    def _maybe_soma_joinid_shape(self) -> Optional[int]:
        """An internal helper method that returns the shape
        value along the ``soma_joinid`` index column, if the ``DataFrame
        has one, else ``None``.


        Lifecycle:
            Experimental.
        """
        return self._handle.maybe_soma_joinid_shape

    @property
    def _maybe_soma_joinid_maxshape(self) -> Optional[int]:
        """An internal helper method that returns the maxshape
        value along the ``soma_joinid`` index column, if the ``DataFrame
        has one, else ``None``.

        Lifecycle:
            Experimental.
        """
        return self._handle.maybe_soma_joinid_maxshape

    @property
    def tiledbsoma_has_upgraded_domain(self) -> bool:
        """Returns true if the array has the upgraded resizeable domain feature
        from TileDB-SOMA 1.15: the array was created with this support, or it has
        had ``.tiledbsoma_upgrade_domain`` applied to it.

        Lifecycle:
            Maturing.
        """
        return self._handle.tiledbsoma_has_upgraded_domain

    def resize_soma_joinid_shape(self, newshape: int) -> None:
        """Increases the shape of the dataframe on the ``soma_joinid`` index
        column, if it indeed is an index column, leaving all other index columns
        as-is. If the ``soma_joinid`` is not an index column, no change is made.
        This is a special case of ``upgrade_domain`` (WIP for 1.15), but simpler
        to keystroke, and handles the most common case for dataframe domain
        expansion.  Raises an error if the dataframe doesn't already have a
        domain: in that case please call ``tiledbsoma_upgrade_domain`` (WIP for
        1.15).
        """
        self._handle._handle.resize_soma_joinid_shape(newshape)

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
            Maturing.
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

        sr = clib.SOMADataFrame.open(
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
        """Writes an `Arrow table <https://arrow.apache.org/docs/python/generated/pyarrow.Table.html>`_
        to the persistent object. As duplicate index values are not allowed, index values already
        present in the object are overwritten and new index values are added.

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
            platform_config:
                Pass in parameters for tuning writes. Example:
                platform_config = tiledbsoma.TileDBWriteOptions(
                    **{"sort_coords": False, "consolidate_and_vacuum": True}
                )

        Raises:
            TypeError:
                If the ``values`` parameter is an unsupported type.
            ValueError:
                If the ``values`` parameter is an empty table.
            SOMAError:
                If the object is not open for writing.

        Lifecycle:
            Maturing.
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
        joinid_type = schema.field(SOMA_JOINID).type
        if joinid_type != pa.int64():
            raise ValueError(
                f"{SOMA_JOINID} field must be of type Arrow int64 but is {joinid_type}"
            )
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


def _fill_out_slot_soma_domain(
    slot_domain: AxisDomain,
    index_column_name: str,
    pa_type: pa.DataType,
    dtype: Any,
) -> Tuple[Tuple[Any, Any], bool]:
    """Helper function for _build_tiledb_schema. Given a user-specified domain for a
    dimension slot -- which may be ``None``, or a two-tuple of which either element
    may be ``None`` -- return either what the user specified (if adequate) or
    sensible type-inferred values appropriate to the datatype.

    Returns a boolean for whether the underlying datatype's max range was used.
    """
    saturated_range = False
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
        if len(slot_domain) != 2:
            raise ValueError(
                f"domain must be a two-tuple; got {len(slot_domain)} elements"
            )
        slot_domain = slot_domain[0], slot_domain[1]
    elif isinstance(dtype, str):
        # Core string dims have no extent and no (core) domain.  We return "" here
        # simply so we can pass libtiledbsoma "" for domain and extent, while it
        # will (and must) ignore these when creating the TileDB schema.
        slot_domain = "", ""
    elif np.issubdtype(dtype, NPInteger):
        iinfo = np.iinfo(cast(NPInteger, dtype))
        slot_domain = iinfo.min, iinfo.max - 1
        # Here the slot_domain isn't specified by the user; we're setting it.
        # The SOMA spec disallows negative soma_joinid.
        if index_column_name == SOMA_JOINID:
            slot_domain = (0, 2**63 - 2)
        saturated_range = True
    elif np.issubdtype(dtype, NPFloating):
        finfo = np.finfo(cast(NPFloating, dtype))
        slot_domain = finfo.min, finfo.max
        saturated_range = True

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

    return (slot_domain, saturated_range)


def _find_extent_for_domain(
    index_column_name: str,
    tiledb_create_write_options: TileDBCreateOptions,
    dtype: Any,
    slot_domain: Tuple[Any, Any],
) -> Any:
    """Helper function for _build_tiledb_schema. Returns a tile extent that is
    small enough for the index-column type, and that also fits within the
    user-specified slot domain (if any).
    """

    # Default 2048 mods to 0 for 8-bit types and 0 is an invalid extent
    extent = tiledb_create_write_options.dim_tile(index_column_name)
    if isinstance(dtype, np.dtype) and dtype.itemsize == 1:
        extent = 1

    # Core string dims have no extent and no (core) domain.  We return "" here
    # simply so we can pass libtiledbsoma "" for domain and extent, while it
    # will (and must) ignore these when creating the TileDB schema.
    if isinstance(dtype, str):
        return ""

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


# We need to do this to avoid this error at array-creation time:
#
# Error: Tile extent check failed; domain max expanded to multiple of tile
# extent exceeds max value representable by domain type. Reduce domain max
# by 1 tile extent to allow for expansion.
def _revise_domain_for_extent(
    domain: Tuple[Any, Any], extent: Any, saturated_range: bool
) -> Tuple[Any, Any]:
    if saturated_range:
        return (domain[0], domain[1] - extent)
    else:
        return domain
