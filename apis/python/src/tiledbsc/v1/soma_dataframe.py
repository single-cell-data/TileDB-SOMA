from typing import Any, Iterator, List, Optional, Sequence, TypeVar, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import tiledb

import tiledbsc.v1.util as util

from .logging import log_io
from .soma_collection import SOMACollection
from .tiledb_array import TileDBArray
from .types import NTuple
from .util import tiledb_type_from_arrow_type

ROWID = "soma_rowid"

Slice = TypeVar("Slice", bound=Sequence[int])


class SOMADataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.

    A SOMADataFrame contains a "pseudo-column" called soma_rowid, of type uint64 and domain
    [0,num_rows).  The soma_rowid pseudo-column contains a unique value for each row in the
    SOMADataFrame, and is intended to act as a join key for other objects, such as a SOMANdArray.
    """

    _shape: Optional[NTuple] = None

    # XXX TEMP -- pending discussion on value-filtering with dense arrays
    _is_sparse: Optional[bool]

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[SOMACollection] = None,
    ):
        """
        See also the :class:`TileDBOject` constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)
        self._is_sparse = None

    def create(
        self,
        schema: pa.Schema,
    ) -> None:
        """
        :param schema: Arrow Schema defining the per-column schema. This schema must define all
        columns. The column name ``soma_rowid`` is reserved for the pseudo-column of the same name.
        If the schema includes types unsupported by the SOMA implementation, an error will be
        raised.
        """
        assert ROWID not in schema.names

        self._create_empty(schema)
        self._is_indexed = False
        self._index_column_names = []

        self._common_create()  # object-type metadata etc

    def _create_empty(
        self,
        schema: pa.Schema,
    ) -> None:
        """
        Create a TileDB 1D dense array with uint64 soma_rowid dimension and multiple attributes.
        """

        level = self._tiledb_platform_config.string_dim_zstd_level

        dom = tiledb.Domain(
            tiledb.Dim(
                name=ROWID,
                domain=(0, np.iinfo(np.int64).max),
                tile=2048,  # TODO: PARAMETERIZE
                dtype=np.uint64,
                filters=[tiledb.ZstdFilter(level=level)],
            ),
            ctx=self._ctx,
        )

        attrs = [
            tiledb.Attr(
                name=attr_name,
                dtype=tiledb_type_from_arrow_type(schema.field(attr_name).type),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
            for attr_name in schema.names
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            # XXX TEMP
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

        self._is_sparse = sch.sparse
        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

    def get_shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        if self._shape is None:
            with self._tiledb_open() as A:
                self._shape = A.shape
        return self._shape

    def get_ndims(self) -> int:
        """
        Return number of index columns
        """
        return len(self.get_index_column_names())

    # TODO
    #    def get_schema(self) -> Arrow.Schema:
    #        """
    #        Return data schema, in the form of an Arrow Schema
    #        """

    def get_indexed(self) -> bool:
        return False

    def get_index_column_names(self) -> List[str]:
        return []

    def read(
        self,
        *,
        # TODO: find out how to spell this in a way the type-checker will accept :(
        # ids: Optional[Union[Sequence[int], str, Slice]] = "all",
        ids: Optional[Any] = "all",
        column_names: Optional[Union[Sequence[str], str]] = "all",
        # TODO: partitions,
        # TODO: result_order,
        # TODO: value_filter
        _return_incomplete: Optional[bool] = True,  # XXX TEMP
    ) -> Iterator[pa.RecordBatch]:
        """
        Read a user-defined subset of data, addressed by the dataframe indexing column, optionally
        filtered, and return results as one or more Arrow.RecordBatch.

        :param ids: Which rows to read. Defaults to 'all'.

        :param column_names: the named columns to read and return. Defaults to 'all'.

        :param partitions: an optional ``SOMAReadPartitions`` hint to indicate how results should be
        organized.

        :param result_order: order of read results.  This can be one of 'rowid-ordered' or 'unordered'.

        :param value_filter: an optional [value filter] to apply to the results. Defaults to no
        filter.

        **Indexing**: the `ids` parameter will support, per dimension: a row offset (uint), a
        row-offset range (slice), or a list of both.
        """

        use_all_ids = False
        use_all_column_names = False
        if isinstance(ids, str):
            assert ids == "all"  # Enforce runtime type check
            use_all_ids = True
        if isinstance(column_names, str):
            assert column_names == "all"  # Enforce runtime type check
            use_all_column_names = True

        with self._tiledb_open("r") as A:
            q = A.query(return_arrow=True, return_incomplete=_return_incomplete)

            if use_all_ids:
                if use_all_column_names:
                    iterator = q.df[:]
                else:
                    iterator = q.df[:][column_names]
            else:
                if use_all_column_names:
                    iterator = q.df[ids]
                else:
                    iterator = q.df[ids][column_names]

            if _return_incomplete:
                for df in iterator:
                    batches = df.to_batches()
                    for batch in batches:
                        yield batch
            else:
                yield iterator

    def write(self, values: pa.RecordBatch) -> None:
        """
        Write an Arrow.RecordBatch to the persistent object.

        :param values: An Arrow.RecordBatch containing all columns, including the index columns. The
        schema for the values must match the schema for the :class:`SOMADataFrame`.

        The ``values`` Arrow RecordBatch must contain a ``soma_rowid`` (uint64) column, indicating
        which rows are being written.
        """
        self._shape = None  # cache-invalidate

        assert ROWID in values.schema.names

        # TODO: contiguity check, and/or split into multiple contiguous writes
        # For now: just assert that these _already are_ contiguous and start with 0.
        rowids = [e.as_py() for e in values.column(ROWID)]
        assert len(rowids) > 0

        rowids = sorted(rowids)
        assert rowids[0] == 0

        lo = rowids[0]
        hi = rowids[-1]

        attr_cols_map = {}
        for name in values.schema.names:
            if name != ROWID:
                attr_cols_map[name] = np.asarray(
                    values.column(name).to_pandas(
                        types_mapper=tiledb_type_from_arrow_type,
                    )
                )

        if self._is_sparse:
            # sparse write
            with self._tiledb_open("w") as A:
                A[rowids] = attr_cols_map
        else:
            # dense write
            with self._tiledb_open("w") as A:
                A[lo : (hi + 1)] = attr_cols_map

    # ================================================================
    # ================================================================
    # ================================================================
    def keys(self) -> List[str]:
        """
        TODO
        """
        return self._tiledb_attr_names()

    def __repr__(self) -> str:
        """
        Default display of `SOMADataFrame`.
        """
        return "\n".join(self._repr_aux())

    def _repr_aux(self) -> List[str]:
        lines = [
            self.get_name()
            + " "
            + self.__class__.__name__
            + " "
            + str(self.get_shape())
        ]
        return lines

    # TODO: TEMP
    def to_dataframe(
        self,
        attrs: Optional[Sequence[str]] = None,
        id_column_name: Optional[
            str
        ] = None,  # to rename index to 'obs_id' or 'var_id', if desired, for anndata
    ) -> pd.DataFrame:
        """
        TODO: comment
        """
        with self._tiledb_open() as A:
            df = A.df[:]
            if attrs is not None:
                df = df[attrs]
            df = self._ascii_to_unicode_dataframe_readback(df)
            df.reset_index(inplace=True)
            if id_column_name is not None:
                df.set_index(id_column_name, inplace=True)

            return df

    def _ascii_to_unicode_dataframe_readback(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Implements the 'decode on read' partof our logic
        """
        for k in df:
            dfk = df[k]
            if len(dfk) > 0 and type(dfk[0]) == bytes:
                df[k] = dfk.map(lambda e: e.decode())
        return df

    # TODO: TEMP
    def from_dataframe(
        self,
        dataframe: pd.DataFrame,
        *,
        extent: int = 2048,
        id_column_name: Optional[
            str
        ] = None,  # to rename index to 'obs_id' or 'var_id', if desired, for anndata
    ) -> None:
        """
        Populates the `obs` element of a SOMAExperiment object.

        :param dataframe: `anndata.obs`
        :param extent: TileDB `extent` parameter for the array schema.
        """
        self._shape = None  # cache-invalidate

        offsets_filters = tiledb.FilterList(
            [tiledb.PositiveDeltaFilter(), tiledb.ZstdFilter(level=-1)]
        )
        dim_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])
        attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])

        s = util.get_start_stamp()
        log_io(None, f"{self._indent}START  WRITING {self.get_uri()}")

        # This is the non-indexed bit
        assert ROWID not in dataframe.keys()
        assert len(dataframe.shape) == 2
        # E.g. (80, 7) is 80 rows x 7 columns
        num_rows = dataframe.shape[0]
        dataframe[ROWID] = np.asarray(range(num_rows))

        mode = "ingest"
        if self.exists():
            mode = "append"
            log_io(None, f"{self._indent}Re-using existing array {self.get_uri()}")

        # Make obs_id a data column
        # TODO: rename it from 'index' to 'obs_id' or 'var_id' etc as appropriate
        dataframe.reset_index(inplace=True)
        if id_column_name is not None:
            dataframe.rename(columns={"index": id_column_name}, inplace=True)

        dataframe.set_index(ROWID, inplace=True)

        # ISSUE:
        #
        # TileDB attributes can be stored as Unicode but they are not yet queryable via the TileDB
        # QueryCondition API. While this needs to be addressed -- global collaborators will want to
        # write annotation-dataframe values in Unicode -- until then, to make obs/var data possible
        # to query, we need to store these as ASCII.
        #
        # This is (besides collation) a storage-level issue not a presentation-level issue: At write
        # time, this works — "α,β,γ" stores as "\xce\xb1,\xce\xb2,\xce\xb3"; at read time: since
        # SOMA is an API: utf8-decode those strings when a query is done & give the user back
        # "α,β,γ".
        #
        # CONTEXT:
        # https://github.com/single-cell-data/TileDB-SingleCell/issues/99
        # https://github.com/single-cell-data/TileDB-SingleCell/pull/101
        # https://github.com/single-cell-data/TileDB-SingleCell/issues/106
        # https://github.com/single-cell-data/TileDB-SingleCell/pull/117
        #
        # IMPLEMENTATION:
        # Python types -- float, string, what have you -- appear as dtype('O') which is not useful.
        # Also, `tiledb.from_pandas` has `column_types` but that _forces_ things to string to a
        # particular if they shouldn't be.
        #
        # Instead, we use `dataframe.convert_dtypes` to get a little jump on what `tiledb.from_pandas`
        # is going to be doing anyway, namely, type-inferring to see what is going to be a string.
        #
        # TODO: when UTF-8 attributes are queryable using TileDB-Py's QueryCondition API we can remove this.
        column_types = {}
        for column_name in dataframe.keys():
            dfc = dataframe[column_name]
            if len(dfc) > 0 and type(dfc[0]) == str:
                # Force ASCII storage if string, in order to make obs/var columns queryable.
                column_types[column_name] = np.dtype("S")

        tiledb.from_pandas(
            uri=self.get_uri(),
            dataframe=dataframe,
            name=self.get_name(),
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
            f"Wrote {self._nested_name}",
            util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.get_uri()}"),
        )
