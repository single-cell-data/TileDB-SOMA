from typing import List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import tiledb

import tiledbsc.util as util

from .logging import logger
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup


class AnnotationDataFrame(TileDBArray):
    """
    Nominally for `obs` and `var` data within a soma. These have one string dimension, and multiple attributes.
    """

    dim_name: str

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,
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
        Returns a tuple with the number of rows and number of columns of the `AnnotationDataFrame`.
        The row-count is the number of obs_ids (for `obs`) or the number of var_ids (for `var`).
        The column-count is the number of columns/attributes in the dataframe.
        """
        with self._open("r") as A:
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
    def ids(self) -> List[str]:
        """
        Returns the `obs_ids` in the matrix (for `obs`) or the `var_ids` (for `var`).
        """
        with self._open("r") as A:
            # TileDB string dims are ASCII not UTF-8. Decode them so they readback
            # not like `b"AKR1C3"` but rather like `"AKR1C3"`.
            retval = A.query(attrs=[], dims=[self.dim_name])[:][self.dim_name].tolist()
            return [e.decode() for e in retval]

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display of soma.obs and soma.var.
        """
        return ", ".join(f"'{key}'" for key in self.keys())

    # ----------------------------------------------------------------
    def __len__(self) -> int:
        """
        Implements `len(soma.obs)` and `len(soma.var)`.
        """
        return len(self.ids())

    # ----------------------------------------------------------------
    def keys(self) -> List[str]:
        """
        Returns the column names for the `obs` or `var` dataframe.  For obs and varp, `.keys()` is a
        keystroke-saver for the more general array-schema accessor `attr_names`.
        """
        return self.attr_names()

    # ----------------------------------------------------------------
    def keyset(self) -> Set[str]:
        """
        Same as `.keys` but returns as set.
        """
        return set(self.keys())

    # ----------------------------------------------------------------
    def dim_select(self, ids, attrs=None) -> pd.DataFrame:
        """
        Selects a slice out of the dataframe with specified `obs_ids` (for `obs`) or `var_ids` (for
        `var`).  If `ids` is `None`, the entire dataframe is returned.  Similarly, if `attrs` are
        provided, they're used for the query; else, all attributes are returned.
        """
        with self._open("r") as A:
            if ids is None:
                if attrs is None:
                    df = A.df[:]
                else:
                    df = A.df[:][attrs]
            else:
                if attrs is None:
                    df = A.df[ids]
                else:
                    df = A.df[ids][attrs]

        # We do not need this:
        #   df.set_index(self.dim_name, inplace=True)
        # as long as these arrays (for this class) are written using tiledb.from_pandas which
        # sets this metadata:
        #   >>> A.meta.items()
        #   (('__pandas_index_dims', '{"obs_id": "<U0"}'),)
        # so the set_index is already done for us.

        # TODO: when UTF-8 attributes are queryable using TileDB-Py's QueryCondition API we can remove this.
        # This is the 'decode on read' part of our logic; in dim_select we have the 'encode on write' part.
        # Context: https://github.com/single-cell-data/TileDB-SingleCell/issues/99.
        return self._ascii_to_unicode_dataframe_readback(df)

    # ----------------------------------------------------------------
    def df(self, ids=None, attrs=None) -> pd.DataFrame:
        """
        Keystroke-saving alias for `.dim_select()`. If `ids` are provided, they're used
        to subselect; if not, the entire dataframe is returned. If `attrs` are provided,
        they're used for the query; else, all attributes are returned.
        """
        return self.dim_select(ids, attrs)

    # ----------------------------------------------------------------
    def query(self, query_string, ids=None, attrs=None) -> pd.DataFrame:
        """
        Selects from obs/var using a TileDB-Py `QueryCondition` string such as `cell_type ==
        "blood"`.  If `attrs` is `None`, returns all column names in the dataframe; use `[]` for
        `attrs` to select none of them.  Any column names specified in the `query_string` must be
        included in `attrs` if `attrs` is not `None`.  Returns `None` if the slice is empty.
        """
        if query_string is None:
            return self.dim_select(ids)

        with self._open() as A:
            qc = tiledb.QueryCondition(query_string)
            if attrs is None:
                slice_query = A.query(attr_cond=qc)
                if ids is None:
                    slice_df = slice_query.df[:][:]
                else:
                    slice_df = slice_query.df[ids][:]
            else:
                slice_query = A.query(attr_cond=qc, attrs=attrs)
                if ids is None:
                    slice_df = slice_query.df[:]
                else:
                    slice_df = slice_query.df[ids]
            # This is the 'decode on read' part of our logic; in dim_select we have the 'encode on write' part.
            # Context: https://github.com/single-cell-data/TileDB-SingleCell/issues/99.
            return self._ascii_to_unicode_dataframe_readback(slice_df)

    # ----------------------------------------------------------------
    def _ascii_to_unicode_dataframe_readback(self, df) -> pd.DataFrame:
        """
        Implements the 'decode on read' partof our logic as noted in `dim_select()`.
        """
        for k in df:
            dfk = df[k]
            if len(dfk) > 0 and type(dfk[0]) == bytes:
                df[k] = dfk.map(lambda e: e.decode())
        return df

    # ----------------------------------------------------------------
    def from_dataframe(self, dataframe: pd.DataFrame, extent: int = 2048) -> None:
        """
        Populates the `obs` or `var` subgroup for a SOMA object.

        :param dataframe: `anndata.obs`, `anndata.var`, `anndata.raw.var`.
        :param extent: TileDB `extent` parameter for the array schema.
        """

        offsets_filters = tiledb.FilterList(
            [tiledb.PositiveDeltaFilter(), tiledb.ZstdFilter(level=-1)]
        )
        dim_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])
        attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  WRITING {self.uri}")

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

        mode = "ingest"
        if self.exists():
            mode = "append"
            logger.info(f"{self._indent}Re-using existing array {self.uri}")

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
                column_types[column_name] = np.dtype("S")

        tiledb.from_pandas(
            uri=self.uri,
            dataframe=dataframe,
            name=self.name,
            sparse=True,
            allows_duplicates=self._soma_options.allows_duplicates,
            offsets_filters=offsets_filters,
            attr_filters=attr_filters,
            dim_filters=dim_filters,
            capacity=100000,
            tile=extent,
            column_types=column_types,
            ctx=self._ctx,
            mode=mode,
        )

        self._set_object_type_metadata()

        logger.info(util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.uri}"))
