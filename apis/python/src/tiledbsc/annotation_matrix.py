from typing import Optional, Tuple

import pandas as pd
import tiledb

import tiledbsc.util as util

from .logging import logger
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup


class AnnotationMatrix(TileDBArray):
    """
    Nominally for obsm and varm group elements within a soma.
    """

    dim_name: str  # e.g. 'obs_id' or 'var_id' -- the name of the one string dimension

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,
        dim_name: str,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)
        self.dim_name = dim_name

    # ----------------------------------------------------------------
    def shape(self) -> Tuple[int, int]:
        """
        Returns a tuple with the number of rows and number of columns of the `AnnotationMatrix`.
        The row-count is the number of obs_ids (for `obsm` elements) or the number of var_ids (for
        `varm` elements).  The column-count is the number of columns/attributes in the dataframe.

        Note: currently implemented via data scan -- will be optimized for TileDB core 2.10.
        """
        with self._open() as A:
            # These TileDB arrays are string-dimensioned sparse arrays so there is no '.shape'.
            # Instead we compute it ourselves.  See also:
            # * https://github.com/single-cell-data/TileDB-SingleCell/issues/10
            # * https://github.com/TileDB-Inc/TileDB-Py/pull/1055
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
    def dim_select(self, ids) -> pd.DataFrame:
        """
        Selects a slice out of the array with specified `obs_ids` (for `obsm` elements) or
        `var_ids` (for `varm` elements).  If `ids` is `None`, the entire array is returned.
        """
        if ids is None:
            with self._open() as A:
                df = A.df[:]
        else:
            with self._open() as A:
                df = A.df[ids]
        df.set_index(self.dim_name, inplace=True)
        return df

    # ----------------------------------------------------------------
    def df(self, ids=None) -> pd.DataFrame:
        """
        Keystroke-saving alias for `.dim_select()`. If `ids` are provided, they're used
        to subselect; if not, the entire dataframe is returned.
        """
        return self.dim_select(ids)

    # ----------------------------------------------------------------
    def from_matrix_and_dim_values(self, matrix, dim_values) -> None:
        """
        Populates an array in the obsm/ or varm/ subgroup for a SOMA object.

        :param matrix: anndata.obsm['foo'], anndata.varm['foo'], or anndata.raw.varm['foo'].
        :param dim_values: anndata.obs_names, anndata.var_names, or anndata.raw.var_names.
        """

        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  WRITING {self.uri}")

        if isinstance(matrix, pd.DataFrame):
            self._from_pandas_dataframe(matrix, dim_values)
        else:
            self._numpy_ndarray_or_scipy_sparse_csr_matrix(matrix, dim_values)

        self._set_object_type_metadata()

        logger.info(util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.uri}"))

    # ----------------------------------------------------------------
    def _numpy_ndarray_or_scipy_sparse_csr_matrix(self, matrix, dim_values) -> None:
        # We do not have column names for anndata-provenance annotation matrices.
        # So, if say we're looking at anndata.obsm['X_pca'], we create column names
        # 'X_pca_1', 'X_pca_2', etc.
        (nrow, nattr) = matrix.shape
        attr_names = [self.name + "_" + str(j) for j in range(1, nattr + 1)]

        # Ingest annotation matrices as 1D/multi-attribute sparse arrays
        if self.exists():
            logger.info(f"{self._indent}Re-using existing array {self.uri}")
        else:
            self._create_empty_array([matrix.dtype] * nattr, attr_names)

        self._ingest_data(matrix, dim_values, attr_names)

    # ----------------------------------------------------------------
    def _from_pandas_dataframe(self, df, dim_values) -> None:
        (nrow, nattr) = df.shape
        attr_names = df.columns.values.tolist()

        # Ingest annotation matrices as 1D/multi-attribute sparse arrays
        if self.exists():
            logger.info(f"{self._indent}Re-using existing array {self.uri}")
        else:
            self._create_empty_array(list(df.dtypes), attr_names)

        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:
            A[dim_values] = df.to_dict(orient="list")

    # ----------------------------------------------------------------
    def _create_empty_array(self, matrix_dtypes, attr_names) -> None:
        """
        Create a TileDB 1D sparse array with string dimension and multiple attributes.

        :param matrix_dtypes: For numpy.ndarray, there is a single dtype and this must be
        repeated once per column. For pandas.DataFrame, there is a dtype per column.
        :param attr_names: column names for the dataframe
        """

        # Nominally 'obs_id' or 'var_id'
        level = self._soma_options.string_dim_zstd_level
        dom = tiledb.Domain(
            tiledb.Dim(
                name=self.dim_name,
                domain=(None, None),
                dtype="ascii",
                filters=[tiledb.ZstdFilter(level=level)],
            ),
            ctx=self._ctx,
        )

        attrs = [
            tiledb.Attr(
                attr_name,
                dtype=matrix_dtypes[j],
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
            for j, attr_name in enumerate(attr_names)
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=self._soma_options.allows_duplicates,
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

        tiledb.Array.create(self.uri, sch, ctx=self._ctx)

    # ----------------------------------------------------------------
    def _ingest_data(self, matrix, dim_values, col_names) -> None:
        """
        Convert ndarray/(csr|csc)matrix to a dataframe and ingest into TileDB.

        :param matrix: Matrix-like object coercible to a pandas dataframe.
        :param dim_values: barcode/gene IDs from anndata.obs_names or anndata.var_names
        :param col_names: List of column names.
        """

        assert len(col_names) == matrix.shape[1]

        df = pd.DataFrame(matrix, columns=col_names)

        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:
            A[dim_values] = df.to_dict(orient="list")
