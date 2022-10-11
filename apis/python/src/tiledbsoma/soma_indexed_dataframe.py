from typing import Any, Iterator, List, Literal, Optional, Sequence, TypeVar

import numpy as np
import pandas as pd
import pyarrow as pa
import tiledb

from . import util_arrow, util_tiledb
from .soma_collection import SOMACollectionBase
from .tiledb_array import TileDBArray
from .types import Ids, NTuple, SOMAResultOrder

ROWID = "soma_rowid"

Slice = TypeVar("Slice", bound=Sequence[int])


class SOMAIndexedDataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.

    All ``SOMAIndexedDataFrame`` must contain a column called ``soma_joinid``, of type ``int64``. The ``soma_joinid`` column contains a unique value for each row in the ``SOMAIndexedDataFrame``, and intended to act as a joint key for other objects, such as ``SOMASparseNdArray``.
    """

    _index_column_names: Optional[List[str]]
    _shape: Optional[NTuple] = None
    _is_sparse: Optional[bool]

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
        self._index_column_names = None
        self._is_sparse = None

    @property
    def type(self) -> Literal["SOMAIndexedDataFrame"]:
        return "SOMAIndexedDataFrame"

    def create(
        self,
        schema: pa.Schema,
        index_column_names: Optional[List[str]] = None,
    ) -> "SOMAIndexedDataFrame":
        """
        :param schema: Arrow Schema defining the per-column schema. This schema must define all columns, including columns to be named as index columns. If the schema includes types unsupported by the SOMA implementation, an error will be raised.

        :param index_column_names: A list of column names to use as user-defined index columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must exist in the schema, and at least one index column name is required.
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
        return self

    def _create_empty(
        self,
        schema: pa.Schema,
        index_column_names: List[str],
    ) -> None:
        """
        Create a TileDB 1D sparse array with dimensions and attributes
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
                if np.issubdtype(dtype, np.integer):
                    lo = np.iinfo(dtype).min
                    hi = np.iinfo(dtype).max - 1
                elif np.issubdtype(dtype, np.floating):
                    lo = np.finfo(dtype).min
                    hi = np.finfo(dtype).max
                else:
                    raise TypeError(f"Unsupported dtype {dtype}")
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
            # col-major tile order: so we write ``X`` with row-major tile order.
            tile_order="row-major",
            ctx=self._ctx,
        )
        self._is_sparse = sch.sparse

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
        return len(self.get_index_column_names())

    @property
    def is_indexed(self) -> Literal[True]:
        return True

    def get_index_column_names(self) -> Sequence[str]:
        """
        Return index (dimension) column names.
        """
        # If we've cached the answer, skip the storage read. Especially if the storage is on the
        # cloud, where we'll avoid an HTTP request.
        if self._index_column_names is None:
            assert self.is_indexed
            names = []
            with self._tiledb_open() as A:
                dom = A.domain
                for i in range(dom.ndim):
                    names.append(dom.dim(i).name)
            self._index_column_names = names

        return self._index_column_names

    def read(
        self,
        *,
        # TODO: find out how to spell this in a way the type-checker will accept :(
        # ids: Optional[Union[Sequence[int], str, Slice]] = None,
        ids: Optional[Any] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[SOMAResultOrder] = None,
        # TODO: more arguments
    ) -> Iterator[pa.Table]:
        """
        Read a user-defined subset of data, addressed by the dataframe indexing columns, optionally filtered, and return results as one or more Arrow.Table.

        :param ids: for each index dimension, which rows to read. Defaults to ``None``, meaning no constraint -- all IDs.

        :param column_names: the named columns to read and return. Defaults to ``None``, meaning no constraint -- all column names.

        :param partitions: an optional ``SOMAReadPartitions`` hint to indicate how results should be organized.

        :param result_order: order of read results. This can be one of 'row-major', 'col-major', or 'unordered'.

        :param value_filter: an optional [value filter] to apply to the results. Defaults to no filter.

        **Indexing**: the ``ids`` parameter will support, per dimension: a list of values of the type of the indexed column.
        """
        tiledb_result_order = util_tiledb.tiledb_result_order_from_soma_result_order(
            result_order, accept=["row-major", "column-major", "unordered"]
        )

        # TODO: more about index_column_names
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
        # TODO: platform_config,
    ) -> pa.Table:
        """
        This is a convenience method around ``read``. It iterates the return value from ``read`` and returns a concatenation of all the table-pieces found. Its nominal use is to simply unit-test cases.
        """
        return pa.concat_tables(
            self.read(ids=ids, value_filter=value_filter, column_names=column_names)
        )

    def write(self, values: pa.Table) -> None:
        """
        Write an Arrow.Table to the persistent object. As duplicate index values are not allowed, index values already present in the object are overwritten and new index values are added.

        :param values: An Arrow.Table containing all columns, including the index columns. The schema for the values must match the schema for the ``SOMAIndexedDataFrame``.
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

        dim_cols_list = [list(dim_col) for dim_col in dim_cols_list]
        with self._tiledb_open("w") as A:
            # TODO: find the right syntax for vardims ... it's not the ``*`` operator ...
            # A[*dim_cols_list] = attr_cols_map
            if len(dim_cols_list) == 1:
                A[dim_cols_list[0]] = attr_cols_map
            elif len(dim_cols_list) == 2:
                A[dim_cols_list[0], dim_cols_list[1]] = attr_cols_map
            else:
                raise Exception("ndims >= 2 not currently supported")

    def read_as_pandas(
        self,
        attrs: Optional[Sequence[str]] = None,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> Iterator[pd.DataFrame]:
        """
        For ``to_anndata``, as well as for any interactive use where the user wants a Pandas dataframe.
        """
        raise NotImplementedError("indexed read_as_pandas not implemented yet")

    def read_as_pandas_all(
        self,
        *,
        ids: Optional[Ids] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[str] = None,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> pd.DataFrame:
        raise NotImplementedError("read_as_pandas_all not implemented.")

    def write_from_pandas(
        self,
        dataframe: pd.DataFrame,
        index_column_names: List[str],
        *,
        extent: int = 2048,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> None:
        """
        Populates the ``obs`` element of a SOMAExperiment object.

        :param dataframe: ``anndata.obs``
        :param extent: TileDB ``extent`` parameter for the array schema.
        """
        raise NotImplementedError("indexed write_from_pandas not implemented yet")

    def write_all_from_pandas(
        self,
        dataframe: pd.DataFrame,
        *,
        extent: int = 2048,
        # to rename index to 'obs_id' or 'var_id', if desired, for anndata
        id_column_name: Optional[str] = None,
    ) -> None:
        raise NotImplementedError("write_all_from_pandas not implemented")
