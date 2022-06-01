import tiledb
from .soma_options import SOMAOptions
from .tiledb_group import TileDBGroup
from .assay_matrix_group import AssayMatrixGroup
from .annotation_dataframe import AnnotationDataFrame
from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
import tiledbsc.util as util

import anndata as ad

from typing import Optional
import os


class RawGroup(TileDBGroup):
    """
    Nominally for soma raw.
    """

    X: AssayMatrixGroup
    var: AnnotationDataFrame
    varm: AnnotationMatrixGroup
    varp: AnnotationPairwiseMatrixGroup
    parent_obs: AnnotationDataFrame

    def __init__(
        self,
        uri: str,
        name: str,
        obs: AnnotationDataFrame,  # Nominally a reference to soma.obs
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the `TileDBObject` constructor.
        See `AssayMatrix` for the rationale behind retaining a reference to the `parent_obs` object.
        """
        super().__init__(uri=uri, name=name, parent=parent)
        self.parent_obs = obs

        X_uri = os.path.join(self.uri, "X")
        var_uri = os.path.join(self.uri, "var")
        varm_uri = os.path.join(self.uri, "varm")
        varp_uri = os.path.join(self.uri, "varp")

        self.var = AnnotationDataFrame(uri=var_uri, name="var", parent=self)
        self.X = AssayMatrixGroup(
            uri=X_uri,
            name="X",
            row_dim_name="obs_id",
            col_dim_name="var_id",
            row_dataframe=self.parent_obs,
            col_dataframe=self.var,
            parent=self.var,
        )
        self.varm = AnnotationMatrixGroup(uri=varm_uri, name="varm", parent=self)
        self.varp = AnnotationPairwiseMatrixGroup(
            uri=varp_uri,
            name="varp",
            row_dataframe=self.var,
            col_dataframe=self.var,
            parent=self,
        )

    # ----------------------------------------------------------------
    def from_anndata(self, anndata: ad.AnnData):
        """
        Writes `anndata.raw` to a TileDB group structure.
        """
        if self._verbose:
            s = util.get_start_stamp()
            print(f"{self._indent}START  WRITING {self.uri}")

        # Must be done first, to create the parent directory
        self._create()

        self.var.from_dataframe(dataframe=anndata.raw.var, extent=2048)
        self._add_object(self.var)

        self.X.from_matrix_and_dim_values(
            anndata.raw.X, anndata.obs.index, anndata.raw.var.index
        )
        self._add_object(self.X)

        self.varm.from_matrices_and_dim_values(anndata.raw.varm, anndata.raw.var_names)
        self._add_object(self.varm)

        if self._verbose:
            print(util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.uri}"))

    # ----------------------------------------------------------------
    def to_anndata_raw(self, obs_labels):
        """
        Reads TileDB storage and returns the material for an `anndata.Raw` object.
        The `obs_labels` must be from the parent object.
        """

        var_df = self.var.to_dataframe()
        X_mat = self.X["data"].to_csr_matrix(obs_labels, var_df.index)
        varm = self.varm.to_dict_of_csr()

        return (X_mat, var_df, varm)
