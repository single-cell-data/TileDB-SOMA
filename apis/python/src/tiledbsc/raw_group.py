import tiledb
from .soma_options                     import SOMAOptions
from .tiledb_group                     import TileDBGroup
from .assay_matrix_group               import AssayMatrixGroup
from .annotation_dataframe             import AnnotationDataFrame
from .annotation_matrix_group          import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
import tiledbsc.util as util

import anndata as ad

from typing import Optional
import os

class RawGroup(TileDBGroup):
    """
    Nominally for soma raw.
    """

    X:     AssayMatrixGroup
    var:   AnnotationDataFrame
    varm:  AnnotationMatrixGroup
    varp:  AnnotationPairwiseMatrixGroup

    def __init__(
        self,
        uri: str,
        name: str,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent)

        X_uri    = os.path.join(self.uri, "X")
        var_uri  = os.path.join(self.uri, "var")
        varm_uri = os.path.join(self.uri, "varm")
        varp_uri = os.path.join(self.uri, "varp")

        self.X    = AssayMatrixGroup(uri=X_uri, name="X", row_dim_name='obs_id', col_dim_name='var_id', parent=self)
        self.var  = AnnotationDataFrame(uri=var_uri, name="var", parent=self)
        self.varm = AnnotationMatrixGroup(uri=varm_uri, name="varm", parent=self)
        self.varp = AnnotationPairwiseMatrixGroup(uri=varp_uri, name="varp", parent=self)


    # ----------------------------------------------------------------
    def from_anndata(self, anndata: ad.AnnData):
        """
        Writes anndata.raw to a TileDB group structure.
        """
        if self.verbose:
            s = util.get_start_stamp()
            print(f"{self.indent}START  WRITING {self.uri}")

        # Must be done first, to create the parent directory
        self.open("w")

        self.X.from_matrix(anndata.raw.X, anndata.obs.index, anndata.raw.var.index)
        self.add(self.X)

        self.var.from_dataframe(dataframe=anndata.raw.var, extent=2048)
        self.add(self.var)

        self.varm.from_anndata(anndata.raw.varm, anndata.raw.var_names)
        self.add(self.varm)

        if self.verbose:
            print(util.format_elapsed(s, f"{self.indent}FINISH WRITING {self.uri}"))

        self.close()

    # ----------------------------------------------------------------
    def to_anndata_raw(self, obs_labels):
        """
        Reads TileDB storage and returns the material for an anndata.Raw object.
        obs_labels must be from the parent object.
        """

        var_df = self.var.to_dataframe()
        X_mat  = self.X.data.to_csr_matrix(obs_labels, var_df.index)
        varm   = self.varm.to_dict_of_csr()

        return (X_mat, var_df, varm)
