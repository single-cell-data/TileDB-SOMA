from typing import List, Optional

import tiledb

from .annotation_dataframe import AnnotationDataFrame
from .assay_matrix import AssayMatrix
from .tiledb_group import TileDBGroup


class AssayMatrixGroup(TileDBGroup):
    """
    Nominally for `X` and `raw/X` elements.  You can find element names using soma.X.keys(); you
    access elements using soma.X['data'] etc., or soma.X.data if you prefer.  (The latter syntax is
    possible when the element name doesn't have dashes, dots, etc. in it.)
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
    def keys(self) -> List[str]:
        """
        For `obsm` and `varm`, `.keys()` is a keystroke-saver for the more general group-member
        accessor `._get_member_names()`.
        """
        return self._get_member_names()

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display of soma.X.
        """
        return ", ".join(f"'{key}'" for key in self.keys())

    # ----------------------------------------------------------------
    def __getattr__(self, name) -> AssayMatrix:
        """
        This is called on `soma.X.name` when `name` is not already an attribute.
        This way you can do `soma.X.data` as an alias for `soma.X['data']`.
        """
        with self._open() as G:
            if name not in G:
                raise AttributeError(
                    f"'{self.__class__.__name__}' object has no attribute '{name}'"
                )
        return self[name]

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

    # ----------------------------------------------------------------
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
    def __getitem__(self, name) -> AssayMatrix:
        """
        Returns an `AnnotationMatrix` element at the given name within the group, or None if no such
        member exists.  Overloads the `[...]` operator.
        """

        with self._open("r") as G:
            if name not in G:
                return None

            obj = G[name]  # This returns a tiledb.object.Object.
            if obj.type == tiledb.tiledb.Group:
                raise Exception(
                    "Internal error: found group element where array element was expected."
                )
            if obj.type != tiledb.libtiledb.Array:
                raise Exception(
                    f"Internal error: found group element neither subgroup nor array: type is {str(obj.type)}"
                )
            return AssayMatrix(
                uri=obj.uri,
                name=name,
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                row_dataframe=self.row_dataframe,
                col_dataframe=self.col_dataframe,
                parent=self,
            )

    # ----------------------------------------------------------------
    def __contains__(self, name) -> bool:
        """
        Implements the `in` operator, e.g. `"data" in soma.X`.
        """
        with self._open("r") as G:
            return name in G

    # ----------------------------------------------------------------
    def add_layer_from_matrix_and_dim_values(
        self,
        matrix,
        row_names: str,
        col_names: str,
        layer_name="data",
    ) -> None:
        """
        Populates the `X` or `raw.X` subgroup for a `SOMA` object.  For `X` and `raw.X`, nominally `row_names` will be `anndata.obs_names` and `col_names` will be `anndata.var_names` or `anndata.raw.var_names`.  For `obsp` elements, both will be `anndata.obs_names`; for `varp elements, both will be `anndata.var_names`.
        """

        if matrix is not None:
            # Must be done first, to create the parent directory
            self.create_unless_exists()

            assay_matrix_uri = self._get_child_uri(
                layer_name
            )  # See comments in that function
            assay_matrix = AssayMatrix(
                uri=assay_matrix_uri,
                name=layer_name,
                row_dim_name=self.row_dim_name,
                col_dim_name=self.col_dim_name,
                row_dataframe=self.row_dataframe,
                col_dataframe=self.col_dataframe,
                parent=self,
            )
            assay_matrix.from_matrix_and_dim_values(matrix, row_names, col_names)

            self._add_object(assay_matrix)
