from concurrent.futures import ThreadPoolExecutor
from typing import Optional, Sequence, Set, Tuple, Union

import pandas as pd
import pyarrow as pa
import tiledb

import tiledbsoma.util as util

from .logging import log_io
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .types import Ids


class AnnotationDataFrame(TileDBArray):
    """
    Nominally for ``obs`` and ``var`` data within a soma. These have one string dimension, and multiple attributes.
    """

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,
        *,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        assert name in ["obs", "var"]
        super().__init__(uri=uri, name=name, parent=parent)
        self.dim_name = name + "_id"

    # ----------------------------------------------------------------
    def shape(self) -> Tuple[int, int]:
        """
        Returns a tuple with the number of rows and number of columns of the ``AnnotationDataFrame``.
        The row-count is the number of obs_ids (for ``obs``) or the number of var_ids (for ``var``).
        The column-count is the number of columns/attributes in the dataframe.
        """
        with self._open("r") as A:
            self.dim_name = A.domain.dim(0).name
            # These TileDB arrays are string-dimensioned sparse arrays so there is no '.shape'.
            # Instead we compute it ourselves.  See also:
            # * https://github.com/single-cell-data/TileDB-SingleCell/issues/10
            # * https://github.com/TileDB-Inc/TileDB-Py/pull/1055
            #
            # Also note that this row-count for obs/var is used by the .shape() methods
            # for X, raw.X, obsp, and varp -- see the AssayMatrix.shape method.
            if self.uri.startswith("tiledb://"):
                num_rows = len(
                    A.query(attrs=[], dims=[self.dim_name])[:][self.dim_name].tolist()
                )
            else:
                # This is quicker than the query -- we can use it safely off TileDB Cloud,
                # and if there's just one fragment written.
                fragment_info = tiledb.array_fragments(self.uri, ctx=self._ctx)
                if len(fragment_info) == 1:
                    num_rows = sum(fragment_info.cell_num)
                else:
                    num_rows = len(
                        A.query(attrs=[], dims=[self.dim_name])[:][
                            self.dim_name
                        ].tolist()
                    )
            num_cols = A.schema.nattr
            return (num_rows, num_cols)

    # ----------------------------------------------------------------
    def ids(self) -> Sequence[str]:
        """
        Returns the ``obs_ids`` in the matrix (for ``obs``) or the ``var_ids`` (for ``var``).
        """
        with self._open("r") as A:
            self.dim_name = A.domain.dim(0).name

            # TileDB string dims are ASCII not UTF-8. Decode them so they readback not like
            # `b"AKR1C3"` but rather like `"AKR1C3"`. Update as of
            # https://github.com/TileDB-Inc/TileDB-Py/pull/1304 these dims will read back OK.
            retval = A.query(attrs=[], dims=[self.dim_name])[:][self.dim_name].tolist()

            retval = [e.decode() for e in retval]

            if len(retval) > 0 and isinstance(retval[0], bytes):
                return [e.decode() for e in retval]
            else:
                # list(...) is there to appease the linter which thinks we're returning `Any`
                return list(retval)

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display of soma.obs and soma.var.
        """
        return ", ".join(f"'{key}'" for key in self.keys())

    # ----------------------------------------------------------------
    def __len__(self) -> int:
        """
        Implements ``len(soma.obs)`` and ``len(soma.var)``.
        """
        return len(self.ids())

    # ----------------------------------------------------------------
    def keys(self) -> Sequence[str]:
        """
        Returns the column names for the ``obs`` or ``var`` dataframe.  For obs and varp, ``.keys()`` is a
        keystroke-saver for the more general array-schema accessor ``attr_names``.
        """
        return self.attr_names()

    # ----------------------------------------------------------------
    def keyset(self) -> Set[str]:
        """
        Same as ``.keys`` but returns as set.
        """
        return set(self.keys())

    # ----------------------------------------------------------------
    def dim_select(
        self,
        ids: Optional[Ids],
        attrs: Optional[Sequence[str]] = None,
        *,
        return_arrow: bool = False,
    ) -> Union[pd.DataFrame, pa.Table]:
        """
        Selects a slice out of the dataframe with specified ``obs_ids`` (for ``obs``) or ``var_ids`` (for
        ``var``).  If ``ids`` is ``None``, the entire dataframe is returned.  Similarly, if ``attrs`` are
        provided, they're used for the query; else, all attributes are returned.
        """
        with self._open("r") as A:
            self.dim_name = A.domain.dim(0).name
            query = A.query(return_arrow=return_arrow, attrs=attrs)
            if ids is None:
                df = query.df[:]
            else:
                df = query.df[ids]

        # We do not need this:
        #   df.set_index(self.dim_name, inplace=True)
        # as long as these arrays (for this class) are written using tiledb.from_pandas which
        # sets this metadata:
        #   >>> A.meta.items()
        #   (('__pandas_index_dims', '{"obs_id": "<U0"}'),)
        # so the set_index is already done for us.
        #
        # However if the data was written somehow else (e.g. by tiledbsoma-r) then we do.
        if not return_arrow:
            if isinstance(df.index, pd.RangeIndex) and self.dim_name in df.columns:
                df.set_index(self.dim_name, inplace=True)

        # TODO: when UTF-8 attributes are queryable using TileDB-Py's QueryCondition API we can remove this.
        # This is the 'decode on read' part of our logic; in from_dataframe we have the 'encode on write' part.
        # Context: https://github.com/single-cell-data/TileDB-SingleCell/issues/99.
        if return_arrow:
            return self._ascii_to_unicode_arrow_readback(df)
        else:
            return self._ascii_to_unicode_pandas_readback(df)

    # ----------------------------------------------------------------
    def df(
        self,
        ids: Optional[Ids] = None,
        attrs: Optional[Sequence[str]] = None,
        *,
        return_arrow: bool = False,
    ) -> Union[pd.DataFrame, pa.Table]:
        """
        Keystroke-saving alias for ``.dim_select()``. If ``ids`` are provided, they're used
        to subselect; if not, the entire dataframe is returned. If ``attrs`` are provided,
        they're used for the query; else, all attributes are returned.
        """
        return self.dim_select(ids, attrs, return_arrow=return_arrow)

    # ----------------------------------------------------------------
    def query(
        self,
        query_string: Optional[str],
        ids: Optional[Ids] = None,
        attrs: Optional[Sequence[str]] = None,
        *,
        return_arrow: bool = False,
    ) -> Union[pd.DataFrame, pa.Table]:
        """
        Selects from obs/var using a TileDB-Py ``QueryCondition`` string such as ``cell_type ==
        "blood"``.  If ``attrs`` is ``None``, returns all column names in the dataframe; use ``[]`` for
        ``attrs`` to select none of them.  Any column names specified in the ``query_string`` must be
        included in ``attrs`` if ``attrs`` is not ``None``.  Returns ``None`` if the slice is empty.
        """
        if query_string is None:
            return self.dim_select(ids, attrs=attrs, return_arrow=return_arrow)

        return self._query_aux(
            query_string=query_string, ids=ids, attrs=attrs, return_arrow=return_arrow
        )

    def _query_aux(
        self,
        query_string: Optional[str],
        ids: Optional[Ids] = None,
        attrs: Optional[Sequence[str]] = None,
        *,
        return_arrow: bool = False,
    ) -> Union[pd.DataFrame, pa.Table]:
        """
        Helper method for `query`: as this has multiple `return` statements, it's easiest to track
        elapsed-time stats in a call to this helper.
        """

        with self._open() as A:
            self.dim_name = A.domain.dim(0).name
            if attrs is None:
                slice_query = A.query(cond=query_string, return_arrow=return_arrow)
                if ids is None:
                    df = slice_query.df[:]
                else:
                    df = slice_query.df[ids]
            else:
                slice_query = A.query(
                    cond=query_string, attrs=attrs, return_arrow=return_arrow
                )
                if ids is None:
                    df = slice_query.df[:]
                else:
                    df = slice_query.df[ids]

            # We do not need this:
            #   df.set_index(self.dim_name, inplace=True)
            # as long as these arrays (for this class) are written using tiledb.from_pandas which
            # sets this metadata:
            #   >>> A.meta.items()
            #   (('__pandas_index_dims', '{"obs_id": "<U0"}'),)
            # so the set_index is already done for us.
            #
            # However if the data was written somehow else (e.g. by tiledbsoma-r) then we do.
            if not return_arrow:
                if isinstance(df.index, pd.RangeIndex) and self.dim_name in df.columns:
                    df.set_index(self.dim_name, inplace=True)
                # This is the 'decode on read' part of our logic; in dim_select we have the 'encode on write' part.
                # Context: https://github.com/single-cell-data/TileDB-SingleCell/issues/99.

            if return_arrow:
                return self._ascii_to_unicode_arrow_readback(df)
            else:
                return self._ascii_to_unicode_pandas_readback(df)

    # ----------------------------------------------------------------
    def _ascii_to_unicode_pandas_series_readback(
        self, field_name: str, series: pd.Series
    ) -> Tuple[str, bool, Optional[pd.Series]]:
        """
        Helper method for ``_ascii_to_unicode_pandas_readback``
        """
        if len(series) > 0 and type(series[0]) == bytes:
            return (field_name, True, series.map(lambda e: e.decode()))
        else:
            return (field_name, False, None)

    def _ascii_to_unicode_pandas_readback(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Implements the 'decode on read' part of our logic as noted in ``dim_select()``.
        """
        futures = []
        # Empirically we find this has a bit of a speed-up. Presumably that's because of some NumPy
        # C++ code releasing the GIL.
        with ThreadPoolExecutor() as executor:
            for k in df:
                future = executor.submit(
                    self._ascii_to_unicode_pandas_series_readback, k, df[k]
                )
                futures.append(future)

        for future in futures:
            (k, modified, dfk) = future.result()
            if modified:
                df[k] = dfk

        return df

    # ----------------------------------------------------------------
    def _ascii_to_unicode_arrow_series_readback(
        self, array_number: int, series: Union[pa.Array, pa.ChunkedArray]
    ) -> Tuple[int, bool, Optional[Union[pa.Array, pa.ChunkedArray]]]:
        """
        Helper method for ``_ascii_to_unicode_arrow_readback``
        """
        # pyarrow's way of handling 'bytes'
        if len(series) > 0 and (
            type(series[0]) == pa.LargeBinaryArray
            or type(series[0]) == pa.LargeStringScalar
        ):
            return (array_number, True, series.cast(pa.string()))
        else:
            return (array_number, False, None)

    def _ascii_to_unicode_arrow_readback(self, df: pa.Table) -> pa.Table:
        """
        Implements the 'decode on read' part of our logic as noted in ``dim_select()``.
        """

        array_names = df.column_names
        futures = []
        # Empirically we find this doesn't have much of a speed-up. Presumably that's because of
        # PyArrow Python code holding the GIL. Nonetheless, experiments show it isn't slower so
        # we'll keep the ThreadPoolExecutor logic, which will only get faster pending (hypothetical)
        # future PyArrow C++ work.

        with ThreadPoolExecutor() as executor:
            for array_number in range(df.num_columns):
                future = executor.submit(
                    self._ascii_to_unicode_arrow_series_readback,
                    array_number,
                    df[array_number],
                )
                futures.append(future)

        new_arrays = [None] * df.num_columns
        for future in futures:
            array_number, modified, new_array = future.result()
            if modified:
                new_arrays[array_number] = new_array
            else:
                new_arrays[array_number] = df[array_number]
        return pa.Table.from_arrays(new_arrays, names=array_names)

    # ----------------------------------------------------------------
    def from_dataframe(
        self,
        dataframe: pd.DataFrame,
        *,
        extent: int = 2048,
        ingest_mode: str = "write",
    ) -> None:
        """
        Populates the ``obs`` or ``var`` subgroup for a SOMA object.

        :param dataframe: ``anndata.obs``, ``anndata.var``, ``anndata.raw.var``.
        :param extent: TileDB ``extent`` parameter for the array schema.
        """
        assert ingest_mode in util.INGEST_MODES

        offsets_filters = tiledb.FilterList(
            [tiledb.PositiveDeltaFilter(), tiledb.ZstdFilter(level=-1)]
        )
        dim_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])
        attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])

        s = util.get_start_stamp()
        log_io(None, f"{self._indent}START  WRITING {self.nested_name}")

        # Make the row-names column (barcodes for obs, gene names for var) explicitly named.
        # Otherwise it'll be called '__tiledb_rows'.
        #
        # Before:
        #
        #   >>> anndata.obs
        #                  orig.ident nCount_RNA nFeature_RNA ...
        #   ATGCCAGAACGACT 0          70.0       47           ...
        #   CATGGCCTGTGCAT 0          85.0       52           ...
        #   ...            ...        ...        ...          ...
        #   GGAACACTTCAGAC 0          150.0      30           ...
        #   CTTGATTGATCTTC 0          233.0      76           ...
        #
        # After:
        #
        #   >>> anndata.obs.rename_axis('obs_id')
        #                  orig.ident nCount_RNA nFeature_RNA ...
        #   obs_id
        #   ATGCCAGAACGACT 0          70.0       47           ...
        #   CATGGCCTGTGCAT 0          85.0       52           ...
        #   ...            ...        ...        ...          ...
        #   GGAACACTTCAGAC 0          150.0      30           ...
        #   CTTGATTGATCTTC 0          233.0      76           ...
        dataframe = dataframe.rename_axis(self.dim_name)

        if ingest_mode == "schema_only":
            from_pandas_mode = "schema_only"
        else:
            from_pandas_mode = "ingest"
            if self.exists():
                from_pandas_mode = "append"
                log_io(
                    None, f"{self._indent}Re-using existing array {self.nested_name}"
                )

        if ingest_mode == "schema_only" and self.exists():
            log_io(
                f"Reused {self.nested_name}",
                util.format_elapsed(s, f"{self._indent}REUSED {self.nested_name}"),
            )
            return

        # ISSUE:
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
                column_types[column_name] = "ascii"

        if ingest_mode == "resume" and self.exists():
            # This lets us check for already-ingested chunks, when in resume-ingest mode.
            ned = self._get_non_empty_domain_as_strings(1)
            sorted_dim_values = sorted(list(dfc.index))
            dim_range = ((sorted_dim_values[0], sorted_dim_values[-1]),)
            if self._chunk_is_contained_in(dim_range, ned):
                self._set_object_type_metadata()
                log_io(
                    f"Skipped {self.nested_name}",
                    util.format_elapsed(
                        s, f"{self._indent}SKIPPED WRITING {self.nested_name}"
                    ),
                )
                return

        tiledb.from_pandas(
            uri=self.uri,
            dataframe=dataframe,
            name=self.name,
            sparse=True,
            allows_duplicates=False,
            offsets_filters=offsets_filters,
            attr_filters=attr_filters,
            dim_filters=dim_filters,
            capacity=self._soma_options.df_capacity,
            cell_order=self._soma_options.df_cell_order,
            tile_order=self._soma_options.df_tile_order,
            tile=extent,
            column_types=column_types,
            ctx=self._ctx,
            mode=from_pandas_mode,
        )

        self._set_object_type_metadata()

        log_io(
            f"Wrote {self.nested_name}",
            util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.nested_name}"),
        )
