from typing import Any, Iterator, Literal, Optional, Sequence, Tuple, TypeVar, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import tiledb

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from . import util, util_arrow
from .constants import SOMA_JOINID
from .query_condition import QueryCondition  # type: ignore
from .soma_collection import SOMACollectionBase
from .tiledb_array import TileDBArray
from .types import Ids, SOMAResultOrder

Slice = TypeVar("Slice", bound=Sequence[int])


class SOMAIndexedDataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.

    All ``SOMAIndexedDataFrame`` must contain a column called ``soma_joinid``, of type ``int64``. The ``soma_joinid`` column contains a unique value for each row in the ``SOMAIndexedDataFrame``, and intended to act as a joint key for other objects, such as ``SOMASparseNdArray``.
    """

    _index_column_names: Union[Tuple[()], Tuple[str, ...]]
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
        self._index_column_names = ()
        self._is_sparse = None

    @property
    def soma_type(self) -> Literal["SOMAIndexedDataFrame"]:
        return "SOMAIndexedDataFrame"

    def create(
        self,
        schema: pa.Schema,
        index_column_names: Sequence[str],
    ) -> "SOMAIndexedDataFrame":
        """
        :param schema: Arrow Schema defining the per-column schema. This schema must define all columns, including columns to be named as index columns. If the schema includes types unsupported by the SOMA implementation, an error will be raised.

        :param index_column_names: A list of column names to use as user-defined index columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must exist in the schema, and at least one index column name is required.
        """
        schema = _validate_schema(schema, index_column_names)
        self._create_empty(schema, index_column_names)
        self._is_indexed = True
        self._index_column_names = tuple(index_column_names)

        self._common_create()  # object-type metadata etc
        return self

    def _create_empty(
        self,
        schema: pa.Schema,
        index_column_names: Sequence[str],
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
            lo: Any = None
            hi: Any = None
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
            allows_duplicates=False,
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
        Returns the names of the columns when read back as a dataframe.
        """
        return self._tiledb_array_keys()

    @property
    def is_indexed(self) -> Literal[True]:
        return True

    def get_index_column_names(self) -> Sequence[str]:
        """
        Return index (dimension) column names.
        """
        # If we've cached the answer, skip the storage read. Especially if the storage is on the
        # cloud, where we'll avoid an HTTP request.
        if self._index_column_names == ():
            assert self.is_indexed
            self._index_column_names = self._tiledb_dim_names()

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
        with self._tiledb_open("r") as A:
            query_condition = None
            if value_filter is not None:
                query_condition = QueryCondition(value_filter)

            # As an arg to this method, `column_names` is optional-None. For the pybind11
            # code it's optional-[].
            lib_column_names = [] if column_names is None else column_names

            sr = clib.SOMAReader(
                self._uri,
                name=self.__class__.__name__,
                schema=A.schema,  # query_condition needs this
                column_names=lib_column_names,
                query_condition=query_condition,
            )

            if ids is not None:
                sr.set_dim_points(A.schema.domain.dim(0).name, util.ids_to_list(ids))

            # TODO: platform_config
            # TODO: batch_size
            # TODO: result_order
            sr.submit()

            while arrow_table := sr.read_next():
                yield arrow_table  # XXX what other post-processing

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
            self.read(
                ids=ids,
                value_filter=value_filter,
                column_names=column_names,
                result_order=result_order,
            )
        )

    def write(self, values: pa.Table) -> None:
        """
        Write an Arrow.Table to the persistent object. As duplicate index values are not allowed, index values already present in the object are overwritten and new index values are added.

        :param values: An Arrow.Table containing all columns, including the index columns. The schema for the values must match the schema for the ``SOMAIndexedDataFrame``.
        """
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
                raise Exception("ndim >= 2 not currently supported")

    def read_as_pandas(
        self,
        ids: Optional[Any] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[SOMAResultOrder] = None,
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
        ids: Optional[Ids] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[SOMAResultOrder] = None,
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

    def write_from_pandas(
        self,
        dataframe: pd.DataFrame,
    ) -> None:
        self.write(pa.Table.from_pandas(dataframe))


def _validate_schema(schema: pa.Schema, index_column_names: Sequence[str]) -> pa.Schema:
    """
    Handle default column additions (eg, soma_joinid) and error checking on required columns.

    Returns a schema, which may be modified by the addition of required columns.
    """
    if not index_column_names:
        raise ValueError("SOMAIndexedDataFrame requires one or more index columns")

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
                "SOMAIndexedDataFrame schema may not contain fields with name prefix `soma_`"
            )

    # verify that all index_column_names are present in the schema
    schema_names_set = set(schema.names)
    for index_column_name in index_column_names:
        if index_column_name.startswith("soma_") and index_column_name != SOMA_JOINID:
            raise ValueError(
                "SOMAIndexedDataFrame schema may not contain fields with name prefix `soma_`"
            )
        if index_column_name not in schema_names_set:
            raise ValueError("All index names must be dataframe schema")

    return schema
