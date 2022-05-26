import tiledb
from .assay_matrix import AssayMatrix
from .annotation_dataframe import AnnotationDataFrame
from .tiledb_group import TileDBGroup
from .soma_options import SOMAOptions

from typing import Optional

import os


class AssayMatrixGroup(TileDBGroup):
    """
    Nominally for soma X and raw/X.
    """

    data: AssayMatrix

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name: str,
        row_dim_name: str,  # obs_id for X, obs_id_i for obsp; var_id_i for varp
        col_dim_name: str,  # var_id for X, obs_id_j for obsp; var_id_j for varp
        row_dataframe: AnnotationDataFrame,  # Nominally a reference to soma.obs
        col_dataframe: AnnotationDataFrame,  # Nominally a reference to soma.var
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        See `AssayMatrix` for the rationale behind retaining references to the `row_dataframe` and
        `col_dataframe` objects.
        """
        super().__init__(uri=uri, name=name, parent=parent)

        assay_matrix_uri = os.path.join(self.uri, "data")
        self.data = AssayMatrix(
            uri=assay_matrix_uri,
            name="data",
            row_dim_name=row_dim_name,
            col_dim_name=col_dim_name,
            row_dataframe=row_dataframe,
            col_dataframe=col_dataframe,
            parent=self,
        )

    # ----------------------------------------------------------------
    def from_matrix_and_dim_values(self, matrix, row_names, col_names) -> None:
        """
        Populates the X/ or raw/X subgroup for a SOMA object.  Nominally row_names will be
        anndata.obs_names and col_names will be anndata.var_names or anndata.raw.var_names.
        """

        if matrix is not None:
            # Must be done first, to create the parent directory
            self._create()
            self.data.from_matrix_and_dim_values(matrix, row_names, col_names)
            self._add_object(self.data)
