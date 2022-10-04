from typing import Any, Iterator, List, Literal, Optional, Sequence, TypeVar

import numpy as np
import pandas as pd
import pyarrow as pa
import tiledb

from . import util, util_arrow, util_tiledb
from .logging import log_io
from .soma_collection import SOMACollectionBase
from .tiledb_array import TileDBArray
from .types import Ids, NTuple, SOMAResultOrder

ROWID = "soma_rowid"

Slice = TypeVar("Slice", bound=Sequence[int])


class SOMADataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.

    A ``SOMADataFrame`` contains a "pseudo-column" called ``soma_rowid``, of type uint64 and domain [0,num_rows).  The ``soma_rowid`` pseudo-column contains a unique value for each row in the ``SOMADataFrame``, and is intended to act as a join key for other objects, such as a ``SOMASparseNdArray``.
    """

    _shape: Optional[NTuple] = None
    _cached_is_sparse: Optional[bool]
    _index_column_names: List[str]

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[SOMACollectionBase[Any]] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        See also the ``TileDBOject`` constructor.
        """
        super().__init__(uri=uri, parent=parent, ctx=ctx)
        self._cached_is_sparse = None

    @property
    def type(self) -> Literal["SOMADataFrame"]:
        return "SOMADataFrame"

    def create(
        self,
        schema: pa.Schema,
    ) -> "SOMADataFrame":
        """
        :param schema: Arrow Schema defining the per-column schema. This schema must define all columns. The column name ``soma_rowid`` is reserved for the pseudo-column of the same name.  If the schema includes types unsupported by the SOMA implementation, an error will be raised.
        """
        assert ROWID not in schema.names

        self._create_empty(schema)
        self._is_indexed = False
        self._index_column_names = []

        self._common_create()  # object-type metadata etc
        return self

    def _create_empty(
        self,
        schema: pa.Schema,
    ) -> None:
        """
        Create a TileDB 1D dense array with uint64 ``soma_rowid`` dimension and multiple attributes.
        """

        level = self._tiledb_platform_config.string_dim_zstd_level

        dom = tiledb.Domain(
            tiledb.Dim(
                name=ROWID,
                domain=(0, np.iinfo(np.uint64).max - 1),
                tile=2048,  # TODO: PARAMETERIZE
                dtype=np.uint64,
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
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            # TODO: pending tiledb issue involving dense dataframes
            # sparse=False,
            sparse=True,
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
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

        self._cached_is_sparse = sch.sparse
        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

    def keys(self) -> Sequence[str]:
        """
        Returns the names of the columns when read back as a dataframe.  TODO: make it clear whether or not this will read back ``soma_rowid`` / ``soma_joinid``.
        """
        return self._tiledb_attr_names()

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndims``.
        """
        if self._shape is None:
            with self._tiledb_open() as A:
                self._shape = A.shape
        return self._shape

    @property
    def ndims(self) -> int:
        """
        Return number of index columns.
        """
        return len(self.keys())

    @property
    def is_indexed(self) -> Literal[False]:
        return False

    def get_index_column_names(self) -> Sequence[str]:
        return []

    def read(
        self,
        *,
        # TODO: find the right syntax to get the typechecker to accept args like ``ids=slice(0,10)``
        # ids: Optional[Union[Sequence[int], Slice]] = None,
        ids: Optional[Any] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[SOMAResultOrder] = None,
        # TODO: batch_size
        # TODO: partition,
        # TODO: platform_config,
    ) -> Iterator[pa.Table]:
        """
        Read a user-defined subset of data, addressed by the dataframe indexing column, optionally filtered, and return results as one or more ``Arrow.Table``.

        :param ids: Which rows to read. Defaults to ``None``, meaning no constraint -- all rows.

        :param column_names: the named columns to read and return. Defaults to ``None``, meaning no constraint -- all column names.

        :param partitions: an optional ``SOMAReadPartitions`` hint to indicate how results should be organized.

        :param result_order: order of read results.  This can be one of 'row-major', 'col-major', or 'unordered'.

        :param value_filter: an optional [value filter] to apply to the results. Defaults to no filter.

        **Indexing**: the ``ids`` parameter will support, per dimension: a row offset (uint), a row-offset range (slice), or a list of both.
        """
        tiledb_result_order = util_tiledb.tiledb_result_order_from_soma_result_order(
            result_order, accept=["rowid-ordered", "unordered"]
        )

        with self._tiledb_open("r") as A:
            dim_names, attr_names = util_tiledb.split_column_names(
                A.schema, column_names
            )
            if value_filter is None:
                query = A.query(
                    return_arrow=True,
                    return_incomplete=True,
                    order=tiledb_result_order,
                    dims=dim_names,
                    attrs=attr_names,
                )
            else:
                qc = tiledb.QueryCondition(value_filter)
                query = A.query(
                    return_arrow=True,
                    return_incomplete=True,
                    attr_cond=qc,
                    order=tiledb_result_order,
                    dims=dim_names,
                    attrs=attr_names,
                )

            if ids is None:
                iterator = query.df[:]
            else:
                iterator = query.df[ids]

            for table in iterator:
                yield table

    def read_all(
        self,
        *,
        # TODO: find the right syntax to get the typechecker to accept args like ``ids=slice(0,10)``
        # ids: Optional[Union[Sequence[int], Slice]] = None,
        ids: Optional[Any] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[SOMAResultOrder] = None,
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

    def _get_is_sparse(self) -> bool:
        if self._cached_is_sparse is None:

            # Simpler would be:
            # if self.exists():
            #     with self._tiledb_open("r") as A:
            #         self._cached_is_sparse = A.schema.sparse
            # but that has _two_ HTTP round trips in the tiledb-cloud case.
            # This way, there is only one.
            try:
                with self._tiledb_open("r") as A:
                    self._cached_is_sparse = A.schema.sparse
            except tiledb.TileDBError as e:
                raise Exception(f"could not read array schema at {self._uri}") from e

        return self._cached_is_sparse

    def write(self, values: pa.Table) -> None:
        """
        Write an Arrow.Table to the persistent object.

        :param values: An Arrow.Table containing all columns, including the index columns. The schema for the values must match the schema for the ``SOMADataFrame``.

        The ``values`` Arrow Table must contain a ``soma_rowid`` (uint64) column, indicating which rows are being written.
        """
        self._shape = None  # cache-invalidate

        assert ROWID in values.schema.names

        # TODO: contiguity check, and/or split into multiple contiguous writes
        # For now: just assert that these _already are_ contiguous and start with 0.
        rowids = [e.as_py() for e in values.column(ROWID)]

        attr_cols_map = {}
        for name in values.schema.names:
            if name != ROWID:
                attr_cols_map[name] = np.asarray(
                    values.column(name).to_pandas(
                        types_mapper=util_arrow.tiledb_type_from_arrow_type,
                    )
                )

        if self._get_is_sparse():
            # sparse write
            with self._tiledb_open("w") as A:
                A[rowids] = attr_cols_map
        else:
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
        result_order: Optional[SOMAResultOrder] = None,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> Iterator[pd.DataFrame]:
        """
        Reads from SOMA storage into memory.  For ``to_anndata``, as well as for any interactive use where the user wants a Pandas dataframe.  Returns a generator over dataframes for batched read. See also ``read_as_pandas_all`` for a convenience wrapper.

        TODO: params-list
        """
        tiledb_result_order = util_tiledb.tiledb_result_order_from_soma_result_order(
            result_order, accept=["rowid-ordered", "unordered"]
        )

        with self._tiledb_open() as A:
            dim_names, attr_names = util_tiledb.split_column_names(
                A.schema, column_names
            )
            if value_filter is None:
                query = A.query(
                    return_incomplete=True,
                    order=tiledb_result_order,
                    dims=dim_names,
                    attrs=attr_names,
                )
            else:
                qc = tiledb.QueryCondition(value_filter)
                query = A.query(
                    return_incomplete=True,
                    attr_cond=qc,
                    order=tiledb_result_order,
                    dims=dim_names,
                    attrs=attr_names,
                )

            if ids is None:
                iterator = query.df[:]
            else:
                iterator = query.df[ids]

            for df in iterator:

                if id_column_name is not None:
                    df.reset_index(inplace=True)
                    df.set_index(id_column_name, inplace=True)

                # Don't materialize soma_rowid on read
                if (
                    ROWID in df.columns
                    and column_names is not None
                    and ROWID not in column_names
                ):
                    yield df.drop(ROWID, axis=1)
                else:
                    yield df

    def read_as_pandas_all(
        self,
        *,
        ids: Optional[Ids] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[SOMAResultOrder] = None,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> pd.DataFrame:
        """
        This is a convenience method around ``read``. It concatenates all partial read results into a single DataFrame. Its nominal use is to simplify unit-test cases.
        """
        dataframes = []
        generator = self.read_as_pandas(
            ids=ids,
            value_filter=value_filter,
            column_names=column_names,
            id_column_name=id_column_name,
            result_order=result_order,
        )
        for dataframe in generator:
            dataframes.append(dataframe)
        return pd.concat(dataframes)

    def write_from_pandas(
        self,
        dataframe: pd.DataFrame,
        *,
        extent: int = 2048,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> None:
        """
        Writes from memory to SOMA storage. Same as ``write_all_from_pandas``, except this method requires the ``soma_rowid`` column to be present (so it knows where to write data), whereas ``write_all_from_pandas``  will populate ``soma_rowid`` for you as zero-up indices.

        :param dataframe: ``anndata.obs`` for example.
        :param extent: TileDB ``extent`` parameter for the array schema.
        """
        self._shape = None  # cache-invalidate

        offsets_filters = tiledb.FilterList(
            [tiledb.PositiveDeltaFilter(), tiledb.ZstdFilter(level=-1)]
        )
        dim_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])
        attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])

        s = util.get_start_stamp()
        log_io(None, f"{self._indent}START  WRITING {self.uri}")

        assert ROWID in dataframe.keys()
        assert len(dataframe.shape) == 2
        # E.g. (80, 7) is 80 rows x 7 columns

        mode = "ingest"
        if self.exists():
            mode = "append"
            log_io(None, f"{self._indent}Re-using existing array")

        # Make obs_id a data column
        dataframe.reset_index(inplace=True)
        if id_column_name is not None:
            dataframe.rename(columns={"index": id_column_name}, inplace=True)

        dataframe.set_index(ROWID, inplace=True)

        # Force ASCII storage if string, in order to make obs/var columns queryable.
        # TODO: when UTF-8 attributes are fully supported we can remove this.
        column_types = {}
        for column_name in dataframe.keys():
            dfc = dataframe[column_name]
            if len(dfc) > 0 and type(dfc[0]) == str:
                column_types[column_name] = "ascii"

        tiledb.from_pandas(
            uri=self.uri,
            dataframe=dataframe,
            # name=self.name,
            sparse=True,  # TODO
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
            offsets_filters=offsets_filters,
            attr_filters=attr_filters,
            dim_filters=dim_filters,
            capacity=100000,
            tile=extent,
            column_types=column_types,
            ctx=self._ctx,
            mode=mode,
        )

        self._common_create()  # object-type metadata etc

        log_io(
            f"Wrote {self.uri}",
            util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.uri}"),
        )

    def write_all_from_pandas(
        self,
        dataframe: pd.DataFrame,
        *,
        extent: int = 2048,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> None:
        """
        Writes from memory to SOMA storage. Same as ``write_from_pandas``, except ``write_from_pandas`` requires the ``soma_rowid`` column to be present (so it knows where to write data), whereas this method will populate ``soma_rowid`` for you as zero-up indices.

        :param dataframe: ``anndata.obs`` for example.
        :param extent: TileDB ``extent`` parameter for the array schema.
        """

        assert ROWID not in dataframe.keys()
        assert len(dataframe.shape) == 2
        # E.g. (80, 7) is 80 rows x 7 columns
        num_rows = dataframe.shape[0]
        dataframe[ROWID] = np.asarray(range(num_rows), dtype=np.uint64)

        self.write_from_pandas(dataframe, extent=extent, id_column_name=id_column_name)
