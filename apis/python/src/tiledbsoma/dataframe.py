from typing import Any, Iterator, Literal, Optional, Sequence, TypeVar

import numpy as np
import pandas as pd
import pyarrow as pa
import tiledb

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from . import util, util_arrow
from .collection import CollectionBase
from .constants import SOMA_JOINID, SOMA_ROWID
from .query_condition import QueryCondition  # type: ignore
from .tiledb_array import TileDBArray
from .types import Ids, ResultOrder

Slice = TypeVar("Slice", bound=Sequence[int])


class DataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.

    A ``DataFrame`` contains a "pseudo-column" called ``soma_rowid``, of type int64 and domain [0,num_rows).  The ``soma_rowid`` pseudo-column contains a unique value for each row in the ``DataFrame``, and is intended to act as a join key for other objects, such as a ``SparseNdArray``.
    """

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[CollectionBase[Any]] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        See also the ``TileDBOject`` constructor.
        """
        super().__init__(uri=uri, parent=parent, ctx=ctx)

    @property
    def soma_type(self) -> Literal["SOMADataFrame"]:
        return "SOMADataFrame"

    def create(
        self,
        schema: pa.Schema,
    ) -> "DataFrame":
        """
        :param schema: Arrow Schema defining the per-column schema. This schema must define all columns. The column name ``soma_rowid`` is reserved for the pseudo-column of the same name.  If the schema includes types unsupported by the SOMA implementation, an error will be raised.
        """
        schema = _validate_schema(schema)
        self._create_empty(schema)
        self._is_indexed = False

        self._common_create()  # object-type metadata etc
        return self

    def _create_empty(
        self,
        schema: pa.Schema,
    ) -> None:
        """
        Create a TileDB 1D dense array with int64 ``soma_rowid`` dimension and multiple attributes.
        """

        level = self._tiledb_platform_config.string_dim_zstd_level

        dom = tiledb.Domain(
            tiledb.Dim(
                name=SOMA_ROWID,
                domain=(0, np.iinfo(np.int64).max - 1),
                tile=2048,  # TODO: PARAMETERIZE
                dtype=np.int64,
                filters=[tiledb.ZstdFilter(level=level)],
            ),
            ctx=self._ctx,
        )

        attrs = [
            tiledb.Attr(
                name=attr_name,
                dtype=util_arrow.tiledb_type_from_arrow_type(
                    schema.field(attr_name).type
                ),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
            for attr_name in schema.names
            if attr_name != SOMA_ROWID
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=False,
            allows_duplicates=False,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=100000,
            cell_order="row-major",
            tile_order="row-major",
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

    def keys(self) -> Sequence[str]:
        """
        Returns the names of the columns when read back as a dataframe.
        """
        return self._tiledb_array_keys()

    @property
    def is_indexed(self) -> Literal[False]:
        return False

    def get_index_column_names(self) -> Sequence[str]:
        return ()

    def read(
        self,
        *,
        # TODO: find the right syntax to get the typechecker to accept args like ``ids=slice(0,10)``
        # ids: Optional[Union[Sequence[int], Slice]] = None,
        ids: Optional[Any] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[ResultOrder] = None,
        # TODO: batch_size
        # TODO: partition,
        # TODO: platform_config,
    ) -> Iterator[pa.Table]:
        """
        Read a user-defined subset of data, addressed by the dataframe indexing column, optionally filtered, and return results as one or more ``Arrow.Table``.

        :param ids: Which rows to read. Defaults to ``None``, meaning no constraint -- all rows.

        :param column_names: the named columns to read and return. Defaults to ``None``, meaning no constraint -- all column names.

        :param partitions: an optional ``ReadPartitions`` hint to indicate how results should be organized.

        :param result_order: order of read results.  This can be one of 'row-major', 'col-major', or 'unordered'.

        :param value_filter: an optional [value filter] to apply to the results. Defaults to no filter.

        **Indexing**: the ``ids`` parameter will support, per dimension: a row offset (uint), a row-offset range (slice), or a list of both.
        """
        with self._tiledb_open("r") as A:
            query_condition = None
            if value_filter is not None:
                query_condition = QueryCondition(value_filter)

            sr = clib.SOMAReader(
                self._uri,
                name=self.__class__.__name__,
                schema=A.schema,  # query_condition needs this
                column_names=column_names,
                query_condition=query_condition,
            )

            if ids is not None:
                if isinstance(ids, list):
                    sr.set_dim_points(SOMA_ROWID, ids)
                elif isinstance(ids, pa.ChunkedArray):
                    sr.set_dim_points(SOMA_ROWID, ids)
                elif isinstance(ids, pa.Array):
                    sr.set_dim_points(SOMA_ROWID, pa.chunked_array(ids))
                elif isinstance(ids, slice):
                    lo_hi = util.slice_to_range(ids)
                    if lo_hi is not None:
                        sr.set_dim_ranges(SOMA_ROWID, [lo_hi])

            # TODO: platform_config
            # TODO: batch_size
            # TODO: result_order
            sr.submit()

            # This requires careful handling in the no-data case.
            #
            # When there is at least one table which has non-zero length, we could use the following:
            #   while arrow_table := sr.read_next():
            #     yield arrow_table
            # This respects the caller's expectation that there will always be at least one table in
            # the iterator. (For example, pd.concat(self.read_as_pandas(coords)) will raise an
            # exception if the iterator has zero elements.)
            #
            # But when there is only zero-length data available, the above does *not* respect
            # caller expectations: because a zero-length pyarrow.Table is falsy (not truthy)
            # the while-loop becomes zero-pass.
            #
            # A tempting alternative is to instead write:
            #   for arrow_table in sr.read_next():
            #       yield arrow_table
            # This is correctly one-pass. However, it yields an iterator of pyarrow.ChunkedArray,
            # not an iterator of pyarrow.Table.
            #
            # For this reason, we use the following i > 0 check to guarantee a minimum of one
            # pass through the yielded iterator even when the resulting table is zero-length.
            i = 0
            while True:
                arrow_table = sr.read_next()
                if not arrow_table and i > 0:
                    break
                i += 1
                yield arrow_table

    def read_all(
        self,
        *,
        # TODO: find the right syntax to get the typechecker to accept args like ``ids=slice(0,10)``
        # ids: Optional[Union[Sequence[int], Slice]] = None,
        ids: Optional[Any] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[ResultOrder] = None,
        # TODO: batch_size
        # TODO: partition,
        # TODO: result_order,
        # TODO: platform_config,
    ) -> pa.Table:
        """
        This is a convenience method around ``read``. It concatenates all partial read results into a single Table. Its nominal use is to simplify unit-test cases.
        """
        return pa.concat_tables(
            self.read(
                ids=ids,
                value_filter=value_filter,
                column_names=column_names,
                result_order=result_order,
            )
        )

    def write(self, values: pa.Table) -> None:
        """
        Write an Arrow.Table to the persistent object.

        :param values: An Arrow.Table containing all columns, including the index columns. The schema for the values must match the schema for the ``DataFrame``.

        The ``values`` Arrow Table must contain a ``soma_rowid`` (int64) column, indicating which rows are being written.
        """
        if (
            SOMA_ROWID not in values.schema.names
            or SOMA_JOINID not in values.schema.names
        ):
            raise ValueError(
                f"Write requires both {SOMA_ROWID} and {SOMA_JOINID} to be specified in table."
            )

        # TODO: contiguity check, and/or split into multiple contiguous writes
        # For now: just assert that these _already are_ contiguous and start with 0.
        rowids = [e.as_py() for e in values.column(SOMA_ROWID)]

        attr_cols_map = {}
        for name in values.schema.names:
            if name != SOMA_ROWID:
                attr_cols_map[name] = np.asarray(
                    values.column(name).to_pandas(
                        types_mapper=util_arrow.tiledb_type_from_arrow_type_for_write,
                    )
                )

        # TODO: This was a quick thing to bootstrap some early ingestion tests but needs more thought.
        # In particular, rowids needn't be either zero-up or contiguous.
        assert len(rowids) > 0
        rowids = sorted(rowids)
        assert rowids[0] == 0

        # dense write
        lo = rowids[0]
        hi = rowids[-1]
        with self._tiledb_open("w") as A:
            A[lo : (hi + 1)] = attr_cols_map

    def read_as_pandas(
        self,
        *,
        ids: Optional[Ids] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[ResultOrder] = None,
    ) -> Iterator[pd.DataFrame]:
        for tbl in self.read(
            ids=ids,
            value_filter=value_filter,
            column_names=column_names,
            result_order=result_order,
        ):
            yield tbl.to_pandas()

    def read_as_pandas_all(
        self,
        *,
        ids: Optional[Ids] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[ResultOrder] = None,
    ) -> pd.DataFrame:
        return pd.concat(
            self.read_as_pandas(
                ids=ids,
                value_filter=value_filter,
                column_names=column_names,
                result_order=result_order,
            ),
            ignore_index=True,
        )

    def write_from_pandas(self, dataframe: pd.DataFrame) -> None:
        """
        Write the Pandas DataFrame. The Pandas DataFrame must contain a soma_rowid and soma_joinid
        column of type int64.
        """
        self.write(pa.Table.from_pandas(dataframe))


def _validate_schema(schema: pa.Schema) -> pa.Schema:
    """
    Handle default column additions (eg, soma_rowid) and error checking on required columns.

    Returns a schema, which may be modified by the addition of required columns.
    """
    if SOMA_ROWID in schema.names:
        if schema.field(SOMA_ROWID).type != pa.int64():
            raise TypeError(f"{SOMA_ROWID} field must be of type Arrow int64")
    else:
        # add SOMA_ROWID
        schema = schema.insert(0, pa.field(SOMA_ROWID, pa.int64()))

    if SOMA_JOINID in schema.names:
        if schema.field(SOMA_JOINID).type != pa.int64():
            raise TypeError(f"{SOMA_JOINID} field must be of type Arrow int64")
    else:
        # add SOMA_JOINID
        schema = schema.insert(1, pa.field(SOMA_JOINID, pa.int64()))

    for field_name in schema.names:
        if field_name.startswith("soma_") and field_name not in [
            SOMA_ROWID,
            SOMA_JOINID,
        ]:
            raise ValueError(
                "DataFrame schema may not contain fields with name prefix `soma_`"
            )

    return schema
