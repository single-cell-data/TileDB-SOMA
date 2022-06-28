from typing import Optional

import anndata as ad

import tiledbsc.util as util

from .annotation_dataframe import AnnotationDataFrame
from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
from .assay_matrix_group import AssayMatrixGroup
from .logging import logger
from .tiledb_group import TileDBGroup


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

        X_uri = self._get_child_uri("X")  # See comments in that function
        var_uri = self._get_child_uri("var")
        varm_uri = self._get_child_uri("varm")
        varp_uri = self._get_child_uri("varp")

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
    def from_anndata(self, anndata: ad.AnnData) -> None:
        """
        Writes `anndata.raw` to a TileDB group structure.
        """
        s = util.get_start_stamp()
        logger.info(f"{self._indent}START  WRITING {self.uri}")

        # Must be done first, to create the parent directory
        self.create_unless_exists()

        self.var.from_dataframe(dataframe=anndata.raw.var, extent=2048)
        self._add_object(self.var)

        self.X.add_layer_from_matrix_and_dim_values(
            matrix=anndata.raw.X,
            row_names=anndata.obs.index,
            col_names=anndata.raw.var.index,
            layer_name="data",
        )
        self._add_object(self.X)

        self.varm.create_unless_exists()
        for key in anndata.raw.varm.keys():
            self.varm.add_matrix_from_matrix_and_dim_values(
                anndata.raw.varm[key], anndata.raw.var_names, key
            )
        self._add_object(self.varm)

        logger.info(util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.uri}"))

    # ----------------------------------------------------------------
    def to_anndata_raw(self, obs_labels):
        """
        Reads TileDB storage and returns the material for an `anndata.Raw` object.
        The `obs_labels` must be from the parent object.
        """

        var_df = self.var.df()
        X_mat = self.X["data"].to_csr_matrix(obs_labels, var_df.index)
        varm = self.varm.to_dict_of_csr()

        return (X_mat, var_df, varm)
