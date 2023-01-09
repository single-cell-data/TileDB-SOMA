import collections.abc
from typing import Any, Optional, Sequence, Tuple, TypeVar, Union, cast

import numpy as np
import pyarrow as pa
import tiledb
from typing_extensions import Final, get_args

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from . import tiledb_platform_config as tdbpc
from . import util, util_arrow
from .collection import CollectionBase
from .constants import SOMA_JOINID
from .query_condition import QueryCondition  # type: ignore
from .tiledb_array import TileDBArray
from .types import PlatformConfig, ResultOrder, SparseDataFrameCoordinates
from .util_iter import TableReadIter

Slice = TypeVar("Slice", bound=Sequence[int])


class DataFrame(TileDBArray):
    """
    Represents ``obs``, ``var``, and others.

    All ``DataFrame`` must contain a column called ``soma_joinid``, of type ``int64``. The ``soma_joinid`` column contains a unique value for each row in the ``DataFrame``, and intended to act as a joint key for other objects, such as ``SparseNDArray``.
    """

    _index_column_names: Union[Tuple[()], Tuple[str, ...]]
    _is_sparse: Optional[bool]

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
        self._index_column_names = ()
        self._is_sparse = None

    soma_type: Final = "SOMADataFrame"

    def create(
        self,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID,),
        platform_config: Optional[PlatformConfig] = None,
    ) -> "DataFrame":
        """
        :param schema: Arrow Schema defining the per-column schema. This schema must define all columns, including columns to be named as index columns. If the schema includes types unsupported by the SOMA implementation, an error will be raised.

        :param index_column_names: A list of column names to use as user-defined index columns (e.g., ``['cell_type', 'tissue_type']``). All named columns must exist in the schema, and at least one index column name is required.
        """
        schema = _validate_schema(schema, index_column_names)
        config_wrapper = tdbpc.from_param(platform_config)
        self._create_empty(schema, index_column_names, config_wrapper.create_options())
        self._index_column_names = tuple(index_column_names)

        self._common_create()  # object-type metadata etc
        return self

    def _create_empty(
        self,
        schema: pa.Schema,
        index_column_names: Sequence[str],
        create_options: tdbpc.CreateWrapper,
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

            if dtype == "ascii":
                # Special "dtype" for TileDB-Py
                pass
            elif dtype != str:
                if np.issubdtype(dtype, np.integer):
                    lo = np.iinfo(dtype).min
                    hi = np.iinfo(dtype).max - 1
                elif np.issubdtype(dtype, np.floating):
                    lo = np.finfo(dtype).min
                    hi = np.finfo(dtype).max
                else:
                    raise TypeError(f"Unsupported dtype {dtype}")

            # Default 2048 mods to 0 for 8-bit types and 0 is an invalid extent
            extent = create_options.dim_tile(index_column_name)
            if isinstance(dtype, np.dtype) and dtype.itemsize == 1:
                extent = 64

            dim = tiledb.Dim(
                name=index_column_name,
                domain=(lo, hi),
                tile=extent,
                dtype=dtype,
                filters=create_options.dim_filters(
                    index_column_name, [dict(_type="ZstdFilter", level=level)]
                ),
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
                filters=create_options.attr_filters(attr_name, ["ZstdFilter"]),
                ctx=self._ctx,
            )
            attrs.append(attr)

        cell_order, tile_order = create_options.cell_tile_orders()

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=False,
            offsets_filters=create_options.offsets_filters(),
            capacity=create_options.get("capacity", 100000),
            cell_order=cell_order,
            # As of TileDB core 2.8.2, we cannot consolidate string-indexed sparse arrays with
            # col-major tile order: so we write ``X`` with row-major tile order.
            tile_order=tile_order,
            ctx=self._ctx,
        )
        self._is_sparse = sch.sparse

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

    def keys(self) -> Sequence[str]:
        """
        Returns the names of the columns when read back as a dataframe.
        """
        return self._tiledb_array_keys()

    def index_column_names(self) -> Sequence[str]:
        """
        Return index (dimension) column names.
        """
        # If we've cached the answer, skip the storage read. Especially if the storage is on the
        # cloud, where we'll avoid an HTTP request.
        if self._index_column_names == ():
            self._index_column_names = self._tiledb_dim_names()

        return self._index_column_names

    @property
    def count(self) -> int:
        """
        Return the number of rows in the dataframe. Same as `len(df)`.
        """

        # A.domain.shape at the tiledb level gives us the 0..2^63 range which is not what we want
        num_rows = cast(
            int,
            clib.SOMAReader(
                self.uri,
                platform_config={} if self._ctx is None else self._ctx.config().dict(),
            ).nnz(),
        )

        return num_rows

    def __len__(self) -> int:
        """
        Return the number of rows in the dataframe. Same as `df.count`.
        """
        return self.count

    def read(
        self,
        *,
        ids: Optional[SparseDataFrameCoordinates] = None,
        value_filter: Optional[str] = None,
        column_names: Optional[Sequence[str]] = None,
        result_order: Optional[ResultOrder] = None,
        # TODO: more arguments
    ) -> TableReadIter:
        """
        Read a user-defined subset of data, addressed by the dataframe indexing columns, optionally filtered, and return results as one or more Arrow.Table.

        :param ids: for each index dimension, which rows to read. Defaults to ``None``, meaning no constraint -- all IDs.

        :param column_names: the named columns to read and return. Defaults to ``None``, meaning no constraint -- all column names.

        :param partitions: an optional ``ReadPartitions`` hint to indicate how results should be organized.

        :param result_order: order of read results. This can be one of 'row-major', 'col-major', or 'auto'.

        :param value_filter: an optional [value filter] to apply to the results. Defaults to no filter.

        **Indexing**: the ``ids`` parameter will support, per dimension: a list of values of the type of the indexed column.

        Acceptable ways to index
        ------------------------

        * None
        * A sequence of coordinates is accepted, one per dimension.
        * Sequence length must be at least one and <= number of dimensions.
        * If the sequence contains missing coordinates (length less than number of dimensions),
          then `slice(None)` -- i.e. no constraint -- is assumed for the missing dimensions.
        * Per-dimension, explicitly specified coordinates can be one of: None, a value, a
          list/ndarray/paarray/etc of values, a slice, etc.
        * Slices are doubly inclusive: slice(2,4) means [2,3,4] not [2,3]. Slice steps can only be +1.
          Slices can be `slice(None)`, meaning select all in that dimension, but may not be half-specified:
          `slice(2,None)` and `slice(None,4)` are both unsupported.
        * Negative indexing is unsupported.
        """

        with self._tiledb_open("r") as A:
            query_condition = None
            if value_filter is not None:
                query_condition = QueryCondition(value_filter)

            if result_order is not None and result_order not in get_args(ResultOrder):
                raise ValueError(
                    "result_order must be one of " + ", ".join(get_args(ResultOrder))
                )

            sr = clib.SOMAReader(
                self._uri,
                name=self.__class__.__name__,
                schema=A.schema,  # query_condition needs this
                column_names=column_names,
                query_condition=query_condition,
                platform_config={} if self._ctx is None else self._ctx.config().dict(),
                result_order=(result_order or "auto"),
            )

            if ids is not None:
                if not isinstance(ids, (list, tuple)):
                    raise TypeError(
                        f"ids type {type(ids)} unsupported; expected list or tuple"
                    )
                if len(ids) < 1 or len(ids) > A.schema.domain.ndim:
                    raise ValueError(
                        f"ids {ids} must have length between 1 and ndim ({A.schema.domain.ndim}); got {len(ids)}"
                    )

                for i, dim_ids in enumerate(ids):
                    # Example: ids = [None, 3, slice(4,5)]
                    # dim_ids takes on values None, 3, and slice(4,5) in this loop body.
                    dim_name = A.schema.domain.dim(i).name
                    if dim_ids is None:
                        pass  # No constraint; select all in this dimension
                    elif isinstance(dim_ids, int) or isinstance(dim_ids, str):
                        # TO DO: Support index types other than int and string when we have support
                        # in libtiledbsoma's SOMAReader. See also
                        # https://github.com/single-cell-data/TileDB-SOMA/issues/419
                        sr.set_dim_points(dim_name, [dim_ids])
                    elif isinstance(dim_ids, np.ndarray):
                        if dim_ids.ndim != 1:
                            raise ValueError(
                                f"only 1D numpy arrays may be used to index; got {dim_ids.ndim}"
                            )
                        sr.set_dim_points(dim_name, dim_ids)
                    elif isinstance(dim_ids, slice):
                        ned = A.nonempty_domain()  # None iff the array has no data
                        lo_hi = util.slice_to_range(dim_ids, ned[i]) if ned else None
                        if lo_hi is not None:
                            lo, hi = lo_hi
                            if lo < 0 or hi < 0:
                                raise ValueError(
                                    f"slice start and stop may not be negative; got ({lo}, {hi})"
                                )
                            assert lo <= hi
                            sr.set_dim_ranges(dim_name, [lo_hi])
                        # Else, no constraint in this slot. This is `slice(None)` which is like
                        # Python indexing syntax `[:]`.
                    elif isinstance(
                        dim_ids, (collections.abc.Sequence, pa.Array, pa.ChunkedArray)
                    ):
                        sr.set_dim_points(dim_name, dim_ids)
                    else:
                        raise TypeError(
                            f"dim_ids type {type(dim_ids)} at slot {i} unsupported"
                        )

            # TODO: platform_config
            # TODO: batch_size

        sr.submit()
        return TableReadIter(sr)

    def write(self, values: pa.Table) -> None:
        """
        Write an Arrow.Table to the persistent object. As duplicate index values are not allowed, index values already present in the object are overwritten and new index values are added.

        :param values: An Arrow.Table containing all columns, including the index columns. The schema for the values must match the schema for the ``DataFrame``.
        """
        dim_cols_list = []
        attr_cols_map = {}
        dim_names_set = self.index_column_names()
        n = None

        for name in values.schema.names:
            n = len(values.column(name))
            if name in dim_names_set:
                dim_cols_list.append(values.column(name).to_pandas())
            else:
                attr_cols_map[name] = values.column(name).to_pandas()
        assert n is not None

        dim_cols_tuple = tuple(dim_cols_list)
        with self._tiledb_open("w") as A:
            A[dim_cols_tuple] = attr_cols_map


def _validate_schema(schema: pa.Schema, index_column_names: Sequence[str]) -> pa.Schema:
    """
    Handle default column additions (e.g., soma_joinid) and error checking on required columns.

    Returns a schema, which may be modified by the addition of required columns.
    """
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
        assert (
            not index_column_name.startswith("soma_")
            or index_column_name == SOMA_JOINID
        )
        if index_column_name not in schema_names_set:
            schema_names_string = "{}".format(list(schema_names_set))
            raise ValueError(
                f"All index names must be dataframe schema: '{index_column_name}' not in {schema_names_string}"
            )
        # TODO: Pending
        # https://github.com/single-cell-data/TileDB-SOMA/issues/418
        # https://github.com/single-cell-data/TileDB-SOMA/issues/419
        if not schema.field(index_column_name).type in [
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
            pa.string(),
            pa.large_string(),
        ]:
            raise TypeError("Unsupported index type - pending fix #418 and #419")

    return schema
