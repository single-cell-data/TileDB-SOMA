from typing import Dict, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import util, util_ann
from .annotation_dataframe import AnnotationDataFrame
from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
from .assay_matrix_group import AssayMatrixGroup
from .logging import log_io
from .tiledb_group import TileDBGroup
from .types import Labels


class RawGroup(TileDBGroup):
    """
    Nominally for soma raw.
    """

    def __init__(
        self,
        uri: str,
        name: str,
        obs: AnnotationDataFrame,  # Nominally a reference to soma.obs
        *,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        See the ``TileDBObject`` constructor.
        See ``AssayMatrix`` for the rationale behind retaining a reference to the ``parent_obs`` object.
        """
        super().__init__(uri=uri, name=name, parent=parent)
        self.parent_obs = obs

        # See comments in _get_child_uris
        child_uris = self._get_child_uris(["X", "var", "varm", "varp"])

        self.var = AnnotationDataFrame(uri=child_uris["var"], name="var", parent=self)
        self.X = AssayMatrixGroup(
            uri=child_uris["X"],
            name="X",
            row_dim_name="obs_id",
            col_dim_name="var_id",
            row_dataframe=self.parent_obs,
            col_dataframe=self.var,
            parent=self,
        )
        self.varm = AnnotationMatrixGroup(
            uri=child_uris["varm"], name="varm", parent=self
        )
        self.varp = AnnotationPairwiseMatrixGroup(
            uri=child_uris["varp"],
            name="varp",
            row_dataframe=self.var,
            col_dataframe=self.var,
            parent=self,
        )

    # ----------------------------------------------------------------
    def from_anndata(
        self,
        anndata: ad.AnnData,
        X_layer_name: str = "data",
        *,
        ingest_mode: str,
    ) -> None:
        """
        Writes ``anndata.raw`` to a TileDB group structure.
        """
        assert ingest_mode in util.INGEST_MODES

        s = util.get_start_stamp()
        log_io(None, f"{self._indent}START  WRITING {self.nested_name}")

        # Must be done first, to create the parent directory
        self.create_unless_exists()

        self.var.from_dataframe(
            dataframe=util_ann._decategoricalize_obs_or_var(anndata.raw.var),
            extent=2048,
            ingest_mode=ingest_mode,
        )
        self._add_object(self.var)

        self.X.add_layer_from_matrix_and_dim_values(
            matrix=anndata.raw.X[:],  # See comments in io.py
            row_names=anndata.obs.index,
            col_names=anndata.raw.var.index,
            layer_name=X_layer_name,
            ingest_mode=ingest_mode,
        )
        self._add_object(self.X)

        self.varm.create_unless_exists()
        for key in anndata.raw.varm.keys():
            self.varm.add_matrix_from_matrix_and_dim_values(
                util._to_tiledb_supported_array_type(anndata.raw.varm[key]),
                anndata.raw.var_names,
                key,
                ingest_mode=ingest_mode,
            )
        self._add_object(self.varm)

        log_io(
            None,
            util.format_elapsed(s, f"{self._indent}FINISH WRITING {self.nested_name}"),
        )

    # ----------------------------------------------------------------
    def to_anndata_raw(
        self,
        obs_labels: Labels,
        X_layer_name: str = "data",
    ) -> Tuple[sp.csr_matrix, pd.DataFrame, Dict[str, np.ndarray]]:
        """
        Reads TileDB storage and returns the material for an ``anndata.Raw`` object.
        The ``obs_labels`` must be from the parent object.
        """
        var_df = self.var.df()
        X_mat = self.X[X_layer_name]
        assert X_mat is not None
        X_mat = X_mat.to_csr_matrix(obs_labels, var_df.index)
        varm = self.varm.to_dict_of_csr()
        return X_mat, var_df, varm
