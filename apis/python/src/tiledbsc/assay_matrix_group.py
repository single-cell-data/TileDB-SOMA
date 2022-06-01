import tiledb
from .assay_matrix import AssayMatrix
from .annotation_dataframe import AnnotationDataFrame
from .tiledb_group import TileDBGroup
from .soma_options import SOMAOptions

from typing import Optional, List

import os


class AssayMatrixGroup(TileDBGroup):
    """
    Nominally for `X`, `raw/X`, `obsp` elements, and `varp` elements.
    """

    row_dim_name: str
    col_dim_name: str
    row_dataframe: AnnotationDataFrame
    col_dataframe: AnnotationDataFrame

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,  # Nominally "X"
        row_dim_name: str,  # obs_id for X, obs_id_i for obsp; var_id_i for varp
        col_dim_name: str,  # var_id for X, obs_id_j for obsp; var_id_j for varp
        row_dataframe: AnnotationDataFrame,  # Nominally a reference to soma.obs
        col_dataframe: AnnotationDataFrame,  # Nominally a reference to soma.var
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the `TileDBObject` constructor.

        See `AssayMatrix` for the rationale behind retaining references to the `row_dataframe` and
        `col_dataframe` objects.
        """
        super().__init__(uri=uri, name=name, parent=parent)

        self.row_dim_name = row_dim_name
        self.col_dim_name = col_dim_name
        self.row_dataframe = row_dataframe
        self.col_dataframe = col_dataframe

    # ----------------------------------------------------------------
    def keys(self):
        """
        For `obsm` and `varm`, `.keys()` is a keystroke-saver for the more general group-member
        accessor `._get_member_names()`.
        """
        return self._get_member_names()

    # ----------------------------------------------------------------
    def __iter__(self) -> List[AssayMatrix]:
        """
        Implements `for matrix in soma.obsm: ...` and `for matrix in soma.varm: ...`
        """
        retval = []
        for name, uri in self._get_member_names_to_uris().items():
            matrix = AssayMatrix(
                uri=uri,
                name=name,
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                row_dataframe=self.row_dataframe,
                col_dataframe=self.col_dataframe,
                parent=self,
            )
            retval.append(matrix)
        return iter(retval)

    # At the tiledb-py API level, *all* groups are name-indexable.  But here at the tiledbsc-py
    # level, we implement name-indexing only for some groups:
    #
    # * Most soma member references are done using Python's dot syntax. For example, rather than
    #   soma['X'], we have simply soma.X, and likewise, soma.raw.X.  Likewise soma.obs and soma.var.
    #
    # * Index references are supported for obsm, varm, obsp, varp, and uns. E.g.
    #   soma.obsm['X_pca'] or soma.uns['neighbors']['params']['method']
    #
    # * Overloading the `[]` operator at the TileDBGroup level isn't necessary -- e.g. we don't need
    #   soma['X'] when we have soma.X -- but also it causes circular-import issues in Python.
    #
    # * Rather than doing a TileDBIndexableGroup which overloads the `[]` operator, we overload
    #   the `[]` operator separately in the various classes which need indexing. This is again to
    #   avoid circular-import issues, and means that [] on `AnnotationMatrixGroup` will return an
    #   `AnnotationMatrix, [] on `UnsGroup` will return `UnsArray` or `UnsGroup`, etc.
    def __getitem__(self, name):
        """
        Returns an `AnnotationMatrix` element at the given name within the group, or None if no such
        member exists.  Overloads the `[...]` operator.
        """

        # TODO: If TileDB-Py were to support `name in G` the line-count could reduce here.
        with self._open("r") as G:
            try:
                obj = G[name]  # This returns a tiledb.object.Object.
            except:
                return None

            if obj.type == tiledb.tiledb.Group:
                raise Exception(
                    "Internal error: found group element where array element was expected."
                )
            elif obj.type == tiledb.libtiledb.Array:
                return AssayMatrix(
                    uri=obj.uri,
                    name=name,
                    row_dim_name=self.row_dim_name,
                    col_dim_name=self.col_dim_name,
                    row_dataframe=self.row_dataframe,
                    col_dataframe=self.col_dataframe,
                    parent=self,
                )

            else:
                raise Exception(
                    f"Internal error: found group element neither subgroup nor array: type is {str(obj.type)}"
                )

    def __contains__(self, name):
        """
        Implements the `in` operator, e.g. `"data" in soma.X`.
        """
        # TODO: this will get easier once TileDB.group.Group supports `name` in `__contains__`.
        # See SC-18057 and https://github.com/single-cell-data/TileDB-SingleCell/issues/113.
        with self._open("r") as G:
            answer = False
            try:
                # This returns a tiledb.object.Object.
                G[name]
                return True
            except:
                return False

    # ----------------------------------------------------------------
    def from_matrix_and_dim_values(self, matrix, row_names, col_names) -> None:
        """
        Populates the `X` or `raw.X` subgroup for a `SOMA` object.  For `X` and `raw.X`, nominally `row_names` will be `anndata.obs_names` and `col_names` will be `anndata.var_names` or `anndata.raw.var_names`.  For `obsp` elements, both will be `anndata.obs_names`; for `varp elements, both will be `anndata.var_names`.
        """

        if matrix is not None:
            # Must be done first, to create the parent directory
            self._create()

            assay_matrix_uri = os.path.join(self.uri, "data")
            X_data = AssayMatrix(
                uri=assay_matrix_uri,
                name="data",
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                row_dataframe=self.row_dataframe,
                col_dataframe=self.col_dataframe,
                parent=self,
            )
            X_data.from_matrix_and_dim_values(matrix, row_names, col_names)

            self._add_object(X_data)
