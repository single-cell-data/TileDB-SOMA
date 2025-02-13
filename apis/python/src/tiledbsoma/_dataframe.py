# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""
Implementation of a SOMA DataFrame
"""

from __future__ import annotations

import inspect
from typing import (
    Any,
    Dict,
    List,
    Sequence,
    Tuple,
    Union,
    cast,
)

import numpy as np
import pyarrow as pa
import somacore
from somacore import options
from typing_extensions import Self

from . import _arrow_types, _util
from . import pytiledbsoma as clib
from ._constants import SOMA_GEOMETRY, SOMA_JOINID
from ._exception import SOMAError, map_exception_for_create
from ._read_iters import ManagedQuery, TableReadIter
from ._soma_array import SOMAArray
from ._tdb_handles import DataFrameWrapper
from ._types import (
    NPFInfo,
    NPFloating,
    NPIInfo,
    NPInteger,
    OpenTimestamp,
    StatusAndReason,
)
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
        domain: Domain | None = None,
        platform_config: options.PlatformConfig | None = None,
        context: SOMATileDBContext | None = None,
        tiledb_timestamp: OpenTimestamp | None = None,
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
                index column. Each tuple must be a pair consisting of the
                minimum and maximum values storable in the index column. For
                example, if there is a single int64-valued index column, then
                ``domain`` might be ``[(100, 200)]`` to indicate that values
                between 100 and 200, inclusive, can be stored in that column.
                If provided, this sequence must have the same length as
                ``index_column_names``, and the index-column domain will be as
                specified.  If omitted entirely, or if ``None`` in a given
                dimension, the corresponding index-column domain will use an
                empty range, and data writes after that will fail with "A range
                was set outside of the current domain". Unless you have a
                particular reason not to, you should always provide the desired
                `domain` at create time: this is an optional but strongly
                recommended parameter. See also ``change_domain`` which allows
                you to expand the domain after create.
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

            if index_column_name == "soma_joinid":
                lower = slot_core_current_domain[0]
                upper = slot_core_current_domain[1]
                if lower < 0 or upper < 0 or upper < lower:
                    raise ValueError(
                        f"domain for soma_joinid must be non-negative with lower <= upper; got ({lower}, {upper})"
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
    def _maybe_soma_joinid_shape(self) -> int | None:
        """An internal helper method that returns the shape
        value along the ``soma_joinid`` index column, if the ``DataFrame
        has one, else ``None``.


        Lifecycle:
            Experimental.
        """
        return self._handle.maybe_soma_joinid_shape

    @property
    def _maybe_soma_joinid_maxshape(self) -> int | None:
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
        had ``tiledbsoma_upgrade_domain`` applied to it.

        Lifecycle:
            Maturing.
        """
        return self._handle.tiledbsoma_has_upgraded_domain

    def tiledbsoma_resize_soma_joinid_shape(
        self, newshape: int, check_only: bool = False
    ) -> StatusAndReason:
        """Increases the shape of the dataframe on the ``soma_joinid`` index
        column, if it indeed is an index column, leaving all other index columns
        as-is.

        If the ``soma_joinid`` is not an index column, no change is made.  This
        is a special case of ``upgrade_domain``, but simpler to
        keystroke, and handles the most common case for dataframe domain
        expansion.

        Raises an error if the dataframe doesn't already have a domain: in that
        case please call ``tiledbsoma_upgrade_domain``.

        If ``check_only`` is ``True``, returns whether the operation would
        succeed if attempted, and a reason why it would not.
        """
        frame = inspect.currentframe()
        function_name_for_messages = frame.f_code.co_name if frame else "tiledbsoma"

        if check_only:
            return cast(
                StatusAndReason,
                self._handle._handle.can_resize_soma_joinid_shape(
                    newshape,
                    function_name_for_messages=function_name_for_messages,
                ),
            )
        else:
            self._handle._handle.resize_soma_joinid_shape(
                newshape,
                function_name_for_messages=function_name_for_messages,
            )
            return (True, "")

    def tiledbsoma_upgrade_soma_joinid_shape(
        self, newshape: int, check_only: bool = False
    ) -> StatusAndReason:
        """This is like ``upgrade_domain``, but it only applies the specified
        domain update to the ``soma_joinid`` index column. (It's a
        keystroke-saver.) Any other index columns have their domain set to match
        the maxdomain. If the ``soma_joinid`` column is not an index column at
        all, then no action is taken.  If ``check_only`` is ``True``, returns
        whether the operation would succeed if attempted, and a reason why it
        would not.
        """
        frame = inspect.currentframe()
        function_name_for_messages = frame.f_code.co_name if frame else "tiledbsoma"

        if check_only:
            return cast(
                StatusAndReason,
                self._handle._handle.can_upgrade_soma_joinid_shape(
                    newshape,
                    function_name_for_messages=function_name_for_messages,
                ),
            )
        else:
            self._handle._handle.upgrade_soma_joinid_shape(
                newshape,
                function_name_for_messages=function_name_for_messages,
            )
            return (True, "")

    def _upgrade_or_change_domain_helper(
        self, newdomain: Domain, function_name_for_messages: str
    ) -> Any:
        """Converts the user-level tuple of low/high pairs into a pyarrow table suitable for calling libtiledbsoma."""

        # Check user-provided domain against dataframe domain.
        dim_names = self._tiledb_dim_names()
        if len(dim_names) != len(newdomain):
            raise ValueError(
                f"{function_name_for_messages}: requested domain has length {len(dim_names)} but the dataframe's schema has index-column count {len(newdomain)}"
            )

        if any([slot is not None and len(slot) != 2 for slot in newdomain]):
            raise ValueError(
                f"{function_name_for_messages}: requested domain must have low,high pairs, or `None`, in each slot"
            )

        # From the dataframe's schema, extract the subschema for only index columns (TileDB dimensions).
        full_schema = self.schema
        dim_schema_list = []
        for dim_name in dim_names:
            dim_schema_list.append(full_schema.field(dim_name))
        dim_schema = pa.schema(dim_schema_list)

        # Convert the user's tuple of low/high pairs into a dict keyed by index-column name.
        new_domain_dict: Dict[str, Domain] = {}
        for dim_name, new_dom in zip(dim_names, newdomain):
            # Domain can't be specified for strings (core constraint) so let them keystroke that easily.
            if (
                dim_schema.field(dim_name).type
                in [
                    pa.string(),
                    pa.large_string(),
                    pa.binary(),
                    pa.large_binary(),
                ]
                and new_dom is None
            ):
                new_domain_dict[dim_name] = ("", "")  # type: ignore
            else:
                new_domain_dict[dim_name] = tuple(new_dom)  # type: ignore

        # Return this as a pyarrow table. This has n columns where n is the number of
        # index columns, and two rows: one row for the low values and one for the high values.
        return pa.RecordBatch.from_pydict(new_domain_dict, schema=dim_schema)

    def tiledbsoma_upgrade_domain(
        self, newdomain: Domain, check_only: bool = False
    ) -> StatusAndReason:
        """Allows you to set the domain of a SOMA :class:`DataFrame`, when the
        ``DataFrame`` does not have a domain set yet.

        The argument must be a tuple of pairs of low/high values for the desired
        domain, one pair per index column. For string index columns, you must
        offer the low/high pair as `("", "")`, or as ``None``.  If ``check_only``
        is ``True``, returns whether the operation would succeed if attempted,
        and a reason why it would not.

        The discussion at ``change_domain`` applies here in its entirety,
        with the following exception: The ``tiledbsoma_upgrade_domain``
        method is used to apply a ``domain`` to a dataframe created before
        TileDB-SOMA 1.15. The ``change_domain`` method is used only for
        a dataframe that already has a domain set, whether it's an older
        dataframe that has had ``tiledbsoma_upgrade_domain`` applied to it,
        or it's a newer dataframe created by TileDB-SOMA 1.15 or later.
        """
        frame = inspect.currentframe()
        function_name_for_messages = frame.f_code.co_name if frame else "tiledbsoma"

        pyarrow_domain_table = self._upgrade_or_change_domain_helper(
            newdomain,
            function_name_for_messages,
        )

        if check_only:
            return cast(
                StatusAndReason,
                self._handle._handle.can_upgrade_domain(
                    pyarrow_domain_table,
                    function_name_for_messages,
                ),
            )
        else:
            self._handle._handle.upgrade_domain(
                pyarrow_domain_table,
                function_name_for_messages,
            )
            return (True, "")

    def change_domain(
        self, newdomain: Domain, check_only: bool = False
    ) -> StatusAndReason:
        """Allows you to enlarge the domain of a SOMA :class:`DataFrame`, when
        the ``DataFrame`` already has a domain.

        The argument must be a tuple of pairs of low/high values for the desired
        domain, one pair per index column. For string index columns, you must
        offer the low/high pair as `("", "")`, or as ``None``.  If ``check_only``
        is ``True``, returns whether the operation would succeed if attempted,
        and a reason why it would not.

        For example, suppose the dataframe's sole index-column name is
        ``"soma_joinid"`` (which is the default at create).  If the dataframe's
        ``.maxdomain`` is ``((0, 999999),)`` and its ``.domain`` is ``((0,
        2899),)``, this means that ``soma_joinid`` values between 0 and 2899 can
        be read or written; any attempt to read or write ``soma_joinid`` values
        outside this range will result in an error. If you then apply
        ``.change_domain([(0, 5700)])``, then ``.domain`` will
        report ``((0, 5699),)``, and now ``soma_joinid`` values in the range 0
        to 5699 can now be written to the dataframe.

        If you use non-default ``index_column_names`` in the dataframe's
        ``create`` then you need to specify the (low, high) pairs for each
        index column. For example, if the dataframe's ``index_column_names``
        is ``["soma_joinid", "cell_type"]``, then you can upgrade domain using
        ``[(0, 5699), ("", "")]``.

        Lastly, it is an error to try to set the ``domain`` to be smaller than
        ``maxdomain`` along any index column.  The ``maxdomain`` of a dataframe is
        set at creation time, and cannot be extended afterward.

        Lifecycle:
            Maturing.
        """
        frame = inspect.currentframe()
        function_name_for_messages = frame.f_code.co_name if frame else "tiledbsoma"

        pyarrow_domain_table = self._upgrade_or_change_domain_helper(
            newdomain,
            function_name_for_messages,
        )
        if check_only:
            return cast(
                StatusAndReason,
                self._handle._handle.can_change_domain(
                    pyarrow_domain_table,
                    function_name_for_messages,
                ),
            )
        else:
            self._handle._handle.change_domain(
                pyarrow_domain_table, function_name_for_messages
            )
            return (True, "")

    def __len__(self) -> int:
        """Returns the number of rows in the dataframe. Same as ``df.count``."""
        return self.count

    def read(
        self,
        coords: options.SparseDFCoords = (),
        column_names: Sequence[str] | None = None,
        *,
        result_order: options.ResultOrderStr = options.ResultOrder.AUTO,
        value_filter: str | None = None,
        batch_size: options.BatchSize = _UNBATCHED,
        partitions: options.ReadPartitions | None = None,
        platform_config: options.PlatformConfig | None = None,
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

        # TODO: batch_size
        return TableReadIter(
            array=self,
            coords=coords,
            column_names=column_names,
            result_order=_util.to_clib_result_order(result_order),
            value_filter=value_filter,
            platform_config=platform_config,
        )

    def write(
        self, values: pa.Table, platform_config: options.PlatformConfig | None = None
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
        if isinstance(platform_config, TileDBCreateOptions):
            raise ValueError(
                "As of TileDB-SOMA 1.13, the write method takes "
                "TileDBWriteOptions instead of TileDBCreateOptions"
            )
        write_options = TileDBWriteOptions.from_platform_config(platform_config)
        sort_coords = write_options.sort_coords

        for batch in values.to_batches():
            mq = ManagedQuery(self)
            mq._handle.set_array_data(batch)
            mq._handle.submit_write(sort_coords or False)

        if write_options.consolidate_and_vacuum:
            self._handle._handle.consolidate_and_vacuum()

        return self


def _canonicalize_schema(
    schema: pa.Schema,
    index_column_names: Sequence[str],
    required_columns: Sequence[str] = [SOMA_JOINID],
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
    elif SOMA_JOINID in required_columns:
        # add SOMA_JOINID
        schema = schema.append(pa.field(SOMA_JOINID, pa.int64()))

    if SOMA_GEOMETRY in schema.names:
        geometry_type = schema.field(SOMA_GEOMETRY).type
        if geometry_type != pa.binary() and geometry_type != pa.large_binary():
            raise ValueError(
                f"{SOMA_GEOMETRY} field must be of type Arrow binary or large_binary but is {geometry_type}"
            )
        schema.set(
            schema.get_field_index(SOMA_GEOMETRY),
            schema.field(SOMA_GEOMETRY).with_metadata({"dtype": "WKB"}),
        )
    elif SOMA_GEOMETRY in required_columns:
        # add SOMA_GEOMETRY
        schema = schema.append(
            pa.field(SOMA_GEOMETRY, pa.large_binary(), metadata={"dtype": "WKB"})
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


def _fill_out_slot_soma_domain(
    slot_domain: AxisDomain,
    is_max_domain: bool,
    index_column_name: str,
    pa_type: pa.DataType,
    dtype: Any,
) -> Tuple[Tuple[Any, Any], Union[bool, Tuple[bool, ...]]]:
    """Helper function for _build_tiledb_schema. Given a user-specified domain for a
    dimension slot -- which may be ``None``, or a two-tuple of which either element
    may be ``None`` -- return either what the user specified (if adequate) or
    sensible type-inferred values appropriate to the datatype.

    Returns a boolean for whether the underlying datatype's max range was used.
    """
    saturated_range = False
    if index_column_name == SOMA_GEOMETRY:
        # SOMA_GEOMETRY domain should be either a list of None or a list of tuple[float, float]
        axes_lo = []
        axes_hi = []
        if isinstance(slot_domain, list):
            f64info: NPFInfo = np.finfo(np.float64)
            saturated_multi_range = []
            for axis_domain in slot_domain:
                if axis_domain is None:
                    axes_lo.append(f64info.min)
                    axes_hi.append(f64info.max)
                    saturated_multi_range.append(True)
                elif not isinstance(axis_domain, tuple) or len(axis_domain) != 2:
                    raise ValueError("Axis domain should be a tuple[float, float]")
                else:
                    if np.issubdtype(type(axis_domain[0]), NPFloating) or np.issubdtype(
                        type(axis_domain[1]), NPFloating
                    ):
                        raise ValueError("Axis domain should be a tuple[float, float]")

                    axes_lo.append(axis_domain[0])
                    axes_hi.append(axis_domain[1])
                    saturated_multi_range.append(False)
            slot_domain = tuple(axes_lo), tuple(axes_hi)
        else:
            raise ValueError(
                f"{SOMA_GEOMETRY} domain should be either a list of None or a list of tuple[float, float]"
            )

        return (slot_domain, tuple(saturated_multi_range))

    if slot_domain is not None:
        # User-specified; go with it when possible
        if (
            pa_type == pa.string()
            or pa_type == pa.large_string()
            or pa_type == pa.binary()
            or pa_type == pa.large_binary()
        ) and tuple(slot_domain) != ("", ""):
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
        if is_max_domain:
            # Core max domain is immutable. If unspecified, it should be as big
            # as possible since it can never be resized.
            iinfo: NPIInfo = np.iinfo(cast(NPInteger, dtype))
            slot_domain = iinfo.min, iinfo.max - 1
            # Here the slot_domain isn't specified by the user; we're setting it.
            # The SOMA spec disallows negative soma_joinid.
            if index_column_name == SOMA_JOINID:
                slot_domain = (0, 2**63 - 2)
            saturated_range = True
        else:
            # Core current domain is mutable but not shrinkable. If
            # unspecified, it should be as small as possible since it can only
            # be grown, not shrunk.
            #
            # Core current-domain semantics are (lo, hi) with both inclusive,
            # with lo <= hi. This means smallest is (0, 0) which is shape 1,
            # not 0.
            slot_domain = 0, 0
    elif np.issubdtype(dtype, NPFloating):
        if is_max_domain:
            finfo: NPFInfo = np.finfo(cast(NPFloating, dtype))
            slot_domain = finfo.min, finfo.max
            saturated_range = True
        else:
            slot_domain = 0.0, 0.0

    # The `iinfo.min+1` is necessary as of tiledb core 2.15 / tiledb-py 0.21.1
    # since `iinfo.min` maps to `NaT` (not a time), resulting in
    #
    #   TypeError: invalid domain extent, domain cannot be safely cast to
    #   dtype dtype('<M8[s]')
    #
    # The `iinfo.max-delta` is necessary since with iinfo.min being bumped by
    # 1, without subtracting we would get
    #
    #   tiledb.cc.TileDBError: [TileDB::Dimension] Error: Tile extent check
    #   failed; domain max expanded to multiple of tile extent exceeds max
    #   value representable by domain type. Reduce domain max by 1 tile extent
    #   to allow for expansion.
    elif dtype == "datetime64[s]":
        if is_max_domain:
            iinfo = np.iinfo(cast(NPInteger, np.int64))
            slot_domain = np.datetime64(iinfo.min + 1, "s"), np.datetime64(
                iinfo.max - 1000000, "s"
            )
        else:
            slot_domain = np.datetime64(0, "s"), np.datetime64(0, "s")
    elif dtype == "datetime64[ms]":
        if is_max_domain:
            iinfo = np.iinfo(cast(NPInteger, np.int64))
            slot_domain = np.datetime64(iinfo.min + 1, "ms"), np.datetime64(
                iinfo.max - 1000000, "ms"
            )
        else:
            slot_domain = np.datetime64(0, "ms"), np.datetime64(0, "ms")
    elif dtype == "datetime64[us]":
        if is_max_domain:
            iinfo = np.iinfo(cast(NPInteger, np.int64))
            slot_domain = np.datetime64(iinfo.min + 1, "us"), np.datetime64(
                iinfo.max - 1000000, "us"
            )
        else:
            slot_domain = np.datetime64(0, "us"), np.datetime64(0, "us")
    elif dtype == "datetime64[ns]":
        if is_max_domain:
            iinfo = np.iinfo(cast(NPInteger, np.int64))
            slot_domain = np.datetime64(iinfo.min + 1, "ns"), np.datetime64(
                iinfo.max - 1000000, "ns"
            )
        else:
            slot_domain = np.datetime64(0, "ns"), np.datetime64(0, "ns")

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

    if index_column_name == SOMA_GEOMETRY:
        return extent

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
    domain: Tuple[Any, Any], extent: Any, saturated_range: Union[bool, Tuple[bool, ...]]
) -> Tuple[Any, Any]:
    if isinstance(saturated_range, tuple):
        # Handle SOMA_GEOMETRY domain with is tuple[list[float], list[float]]
        if isinstance(domain[1], tuple):
            if len(saturated_range) != len(domain[1]):
                raise ValueError(
                    "Internal error: Saturatin flag length does not match domain size"
                )

            return (
                domain[0],
                [
                    (dim_max - extent) if saturated_range[idx] else dim_max
                    for idx, dim_max in enumerate(domain[1])
                ],
            )

        raise ValueError("Expected a complex domain")
    elif saturated_range:
        return (domain[0], domain[1] - extent)
    else:
        return domain
