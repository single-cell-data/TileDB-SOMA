import tiledb
from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .soma_options import SOMAOptions
import tiledbsc.util as util

import pandas as pd

from typing import Optional


class AnnotationPairwiseMatrix(TileDBArray):
    """
    Nominally for obsp and varp group elements within a soma.
    """

    row_dim_name: str  # e.g. 'obs_id_i' or 'var_id_i'
    col_dim_name: str  # e.g. 'obs_id_j' or 'var_id_j'

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,
        row_dim_name: str,
        col_dim_name: str,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)
        self.row_dim_name = row_dim_name
        self.col_dim_name = col_dim_name

    # ----------------------------------------------------------------
    def shape(self):
        """
        Returns a tuple with the number of rows and number of columns of the `AnnotationPairwiseMatrix`.
        The row-count and column-counts should be the same: dimensions are `obs_id_i,obs_id_j` for
        `obsm` elements, and `var_id_i,var_id_j` for `varm` elements.
        """
        with self._open() as A:
            # These TileDB arrays are string-dimensioned sparse arrays so there is no '.shape'.
            # Instead we compute it ourselves.  See also:
            # * https://github.com/single-cell-data/TileDB-SingleCell/issues/10
            # * https://github.com/TileDB-Inc/TileDB-Py/pull/1055
            return A.df[:].shape  # nnz x 3 -- id_i, id_j, and value

    # ----------------------------------------------------------------
    def dim_select(self, ids):
        """
        Selects a slice out of the array with specified `obs_ids` (for `obsp` elements) or
        `var_ids` (for `varp` elements).  If `ids` is `None`, the entire array is returned.
        """
        with self._open() as A:
            if ids is None:
                return A.df[:, :]
            else:
                return A.df[ids, ids]

    # ----------------------------------------------------------------
    def df(self, ids=None) -> pd.DataFrame:
        """
        Keystroke-saving alias for `.dim_select()`. If `ids` are provided, they're used
        to subselect; if not, the entire dataframe is returned.
        """
        return self.dim_select(ids)

    # ----------------------------------------------------------------
    def from_anndata(self, matrix, dim_values):
        """
        Populates an array in the obsp/ or varp/ subgroup for a SOMA object.

        :param matrix: anndata.obsp['foo'], anndata.varp['foo'], or anndata.raw.varp['foo'].
        :param dim_values: anndata.obs_names, anndata.var_names, or anndata.raw.var_names.
        """

        if self._verbose:
            s = util.get_start_stamp()
            print(f"{self._indent}START  WRITING {self.uri}")

        # We do not have column names for anndata-provenance annotation matrices.
        # So, if say we're looking at anndata.obsm['X_pca'], we create column names
        # 'X_pca_1', 'X_pca_2', etc.
        (nrow, nattr) = matrix.shape
        attr_names = [self.name + "_" + str(j) for j in range(1, nattr + 1)]

        # Ingest annotation matrices as 1D/multi-attribute sparse arrays
        if self.exists():
            if self._verbose:
                print(f"{self._indent}Re-using existing array {self.uri}")
        else:
            self._create_empty_array(matrix.dtype, attr_names)

        self._ingest_data(matrix, dim_values, attr_names)

        if self._verbose:
            print(util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.uri}"))

    # ----------------------------------------------------------------
    def _create_empty_array(self, matrix_dtype, attr_names):
        """
        Create a TileDB 1D sparse array with string dimension and multiple attributes.

        :param matrix_dtype: e.g. anndata.obsm['X_pca'].dtype
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
                dtype=matrix_dtype,
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
            for attr_name in attr_names
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=True,
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
    def _ingest_data(self, matrix, dim_values, col_names):
        """
        Convert ndarray/(csr|csc)matrix to a dataframe and ingest into TileDB.

        :param matrix: Matrix-like object coercible to a pandas dataframe.
        :param dim_values: barcode/gene IDs from anndata.obs_names or anndata.var_names
        :param col_names: List of column names.
        """

        assert len(col_names) == matrix.shape[1]

        df = pd.DataFrame(matrix, columns=col_names)

        with self._open("w") as A:
            A[dim_values] = df.to_dict(orient="list")
