from typing import Any, Generator, Iterator, List, Optional, Sequence, TypeVar, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import tiledb

import tiledbsc.v1.util_arrow as util_arrow

from .soma_collection import SOMACollection
from .tiledb_array import TileDBArray
from .types import NTuple

ROWID = "soma_rowid"

Slice = TypeVar("Slice", bound=Sequence[int])


class SOMAIndexedDataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.

    A SOMAIndexedDataFrame contains a "pseudo-column" called soma_rowid, of type uint64 and domain
    [0,num_rows).  The soma_rowid pseudo-column contains a unique value for each row in the
    SOMAIndexedDataFrame, and is intended to act as a join key for other objects, such as a SOMANdArray.
    """

    _index_column_names: Optional[List[str]]
    _shape: Optional[NTuple] = None
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
        self._index_column_names = None
        self._is_sparse = None

    def create(
        self,
        schema: pa.Schema,
        index_column_names: Optional[List[str]] = None,
    ) -> None:
        """
        :param schema: Arrow Schema defining the per-column schema. This schema must define all
        columns, including columns to be named as index columns. The column name ``soma_rowid`` is
        reserved for the pseudo-column of the same name. If the schema includes types unsupported by
        the SOMA implementation, an error will be raised.

        :param index_column_names: A list of column names to use as user-defined index columns
        (e.g., ``['cell_type', 'tissue_type']``). All named columns must exist in the schema, and at
        least one index column name is required.
        """
        assert index_column_names is not None
        assert len(index_column_names) >= 1

        # assert all index_column_names are present in the schema
        schema_names_set = set(schema.names)
        for index_column_name in index_column_names:
            assert index_column_name in schema_names_set

        assert ROWID not in index_column_names
        assert ROWID not in schema_names_set

        self._create_empty(schema, index_column_names)
        self._is_indexed = True
        self._index_column_names = index_column_names

        self._common_create()  # object-type metadata etc

    def _create_empty(
        self,
        schema: pa.Schema,
        index_column_names: List[str],
    ) -> None:
        """
        Create a TileDB 1D sparse array with string dimension and multiple attributes.
        """

        level = self._tiledb_platform_config.string_dim_zstd_level

        dims = []
        for index_column_name in index_column_names:
            dtype = util_arrow.tiledb_type_from_arrow_type(
                schema.field(index_column_name).type
            )
            # We need domain=(None,None) for string dims
            lo = None
            hi = None
            if dtype != str:
                lo = np.iinfo(dtype).min
                hi = np.iinfo(dtype).max - 1
            dim = tiledb.Dim(
                name=index_column_name,
                domain=(lo, hi),
                tile=2048,  # TODO: PARAMETERIZE
                dtype=dtype,
                filters=[tiledb.ZstdFilter(level=level)],
            )
            dims.append(dim)

        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = []

        attr = tiledb.Attr(
            name=ROWID,
            dtype=np.uint64,
            filters=[tiledb.ZstdFilter()],
            ctx=self._ctx,
        )
        attrs.append(attr)

        for attr_name in schema.names:
            if attr_name in index_column_names:
                continue
            attr = tiledb.Attr(
                name=attr_name,
                dtype=util_arrow.tiledb_type_from_arrow_type(
                    schema.field(attr_name).type
                ),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
            attrs.append(attr)

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=100000,
            cell_order="row-major",
            # As of TileDB core 2.8.2, we cannot consolidate string-indexed sparse arrays with
            # col-major tile order: so we write `X` with row-major tile order.
            tile_order="row-major",
            ctx=self._ctx,
        )
        self._is_sparse = sch.sparse

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

    def __repr__(self) -> str:
        """
        Default display of `SOMAIndexedDataFrame`.
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

    def keys(self) -> List[str]:
        """
        Returns the names of the columns when read back as a dataframe.
        TODO: make it clear whether or not this will read back soma_rowid / soma_joinid.
        """
        return self._tiledb_attr_names()

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
        """
        Return index (dimension) column names.
        """
        # If we've cached the answer, skip the storage read. Especially if the storage is on the
        # cloud, where we'll avoid an HTTP request.
        if self._index_column_names is None:
            if self.get_is_indexed():
                names = []
                with self._tiledb_open() as A:
                    dom = A.domain
                    for i in range(dom.ndim):
                        names.append(dom.dim(i).name)
                self._index_column_names = names
            else:
                self._index_column_names = []
        return self._index_column_names

    def read(
        self,
        *,
        # TODO: find out how to spell this in a way the type-checker will accept :(
        # ids: Optional[Union[Sequence[int], str, Slice]] = None,
        ids: Optional[Any] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Union[Sequence[str], str]] = None,
        # TODO: more arguments
    ) -> Iterator[pa.RecordBatch]:
        """
        Read a user-defined subset of data, addressed by the dataframe indexing columns, optionally
        filtered, and return results as one or more Arrow.RecordBatch.

        :param ids: for each index dimension, which rows to read. Defaults to 'all'.

        :param column_names: the named columns to read and return. Defaults to 'all'.

        :param partitions: an optional ``SOMAReadPartitions`` hint to indicate how results should be
        organized.

        :param result_order: order of read results. This can be one of 'row-major', 'col-major', or
        'unordered'.

        :param value_filter: an optional [value filter] to apply to the results. Defaults to no
        filter.

        **Indexing**: the `ids` parameter will support, per dimension: a list of values of the type
        of the indexed column.
        """
        # TODO: more about index_column_names
        with self._tiledb_open("r") as A:
            if value_filter is None:
                query = A.query(return_arrow=True, return_incomplete=True)
            else:
                qc = tiledb.QueryCondition(value_filter)
                query = A.query(return_arrow=True, return_incomplete=True, attr_cond=qc)

            if ids is None:
                iterator = query.df[:]
                if column_names is not None:
                    iterator = query.df[:][column_names]
            else:
                iterator = query.df[ids]
                if column_names is not None:
                    iterator = query.df[ids][column_names]

            for df in iterator:
                batches = df.to_batches()
                for batch in batches:
                    # XXX COMMENT MORE
                    # This is the 'decode on read' part of our logic; in dim_select we have the
                    # 'encode on write' part.
                    # Context: # https://github.com/single-cell-data/TileDB-SingleCell/issues/99.
                    yield util_arrow.ascii_to_unicode_pyarrow_readback(batch)

    def read_all(
        self,
        *,
        # TODO: find the right syntax to get the typechecker to accept args like `ids=slice(0,10)`
        # ids: Optional[Union[Sequence[int], Slice]] = None,
        ids: Optional[Any] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        # TODO: batch_size
        # TODO: partition,
        # TODO: result_order,
        # TODO: platform_config,
    ) -> pa.RecordBatch:
        """
        This is a convenience method around `read`. It iterates the return value from `read`
        and returns a concatenation of all the record batches found. Its nominal use is to
        simply unit-test cases.
        """
        return util_arrow.concat_batches(
            self.read(ids=ids, value_filter=value_filter, column_names=column_names)
        )

    def write(self, values: pa.RecordBatch) -> None:
        """
        Write an Arrow.RecordBatch to the persistent object. As duplicate index values are not allowed,
        index values already present in the object are overwritten and new index values are added.

        :param values: An Arrow.RecordBatch containing all columns, including the index columns. The
        schema for the values must match the schema for the :class:`SOMAIndexedDataFrame`.
        """
        self._shape = None  # cache-invalidate

        dim_cols_list = []
        attr_cols_map = {}
        dim_names_set = self.get_index_column_names()
        n = None

        for name in values.schema.names:
            n = len(values.column(name))
            if name in dim_names_set:
                dim_cols_list.append(values.column(name).to_pandas())
            else:
                attr_cols_map[name] = values.column(name).to_pandas()
        assert n is not None

        # TODO: temporary
        attr_cols_map[ROWID] = np.asarray(range(n))

        dim_cols_list = [list(dim_col) for dim_col in dim_cols_list]
        with self._tiledb_open("w") as A:
            # TODO: find the right syntax for vardims ... it's not the `*` operator ...
            # A[*dim_cols_list] = attr_cols_map
            if len(dim_cols_list) == 1:
                A[dim_cols_list[0]] = attr_cols_map
            elif len(dim_cols_list) == 2:
                A[dim_cols_list[0], dim_cols_list[1]] = attr_cols_map
            else:
                raise Exception("ndims >= 2 not currently supported")

    def to_pandas(
        self,
        attrs: Optional[Sequence[str]] = None,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> Generator:
        """
        For `to_anndata`, as well as for any interactive use where the user wants a Pandas dataframe.
        """
        raise Exception("indexed to_pandas not implemented yet")

    def from_pandas(
        self,
        dataframe: pd.DataFrame,
        index_column_names: List[str],
        *,
        extent: int = 2048,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> None:
        """
        Populates the `obs` element of a SOMAExperiment object.

        :param dataframe: `anndata.obs`
        :param extent: TileDB `extent` parameter for the array schema.
        """
        raise Exception("indexed from_pandas not implemented yet")
