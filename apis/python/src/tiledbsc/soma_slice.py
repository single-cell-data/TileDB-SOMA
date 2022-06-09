import os
from typing import Optional, Union, Dict

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scanpy
import scipy
import tiledb

from tiledbsc import util
from tiledbsc import util_ann

from .soma_options import SOMAOptions
from .tiledb_group import TileDBGroup
from .assay_matrix_group import AssayMatrixGroup
from .annotation_dataframe import AnnotationDataFrame
from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
from .raw_group import RawGroup
from .uns_group import UnsGroup


class SOMASlice(TileDBGroup):
    """
    In-memory-only object for ephemeral extracting out of a SOMA. Can be used to _construct_ a SOMA
    but is not a SOMA (which would entail out-of-memory storage).  Nothing more than a collection of
    pandas.DataFrame objects. No raw or uns.
    """

    X: pd.DataFrame
    obs: pd.DataFrame
    var: pd.DataFrame
    obsm: Dict[str, pd.DataFrame]
    varm: Dict[str, pd.DataFrame]
    obsp: Dict[str, pd.DataFrame]
    varp: Dict[str, pd.DataFrame]

    # ----------------------------------------------------------------
    def __init__(
        self,
        X: pd.DataFrame,
        obs: pd.DataFrame,
        var: pd.DataFrame,
    ):
        """
        Constructs an in-memory `SOMASlice` object. This is a simple collection of obs, var, and X dataframes.
        """
        assert isinstance(X, pd.DataFrame)
        assert isinstance(obs, pd.DataFrame)
        assert isinstance(var, pd.DataFrame)

        self.X = X
        self.obs = obs
        self.var = var

    # ----------------------------------------------------------------
    def to_anndata(self) -> ad.AnnData:
        """
        Constructs an `AnnData` object from the current `SOMASlice` object.
        """

        # Make obs_id and var_id accessible as columns
        self.X.reset_index(inplace=True)

        sparseX = util.X_and_ids_to_sparse_matrix(
            self.X,
            "obs_id",  # row_dim_name
            "var_id",  # col_dim_name
            "value",  # attr_name
            self.obs.index,
            self.var.index,
        )

        ann = ad.AnnData(
            X=sparseX,
            obs=self.obs,
            var=self.var,
            dtype=sparseX.dtype,
        )

        ann.obs_names_make_unique()
        ann.var_names_make_unique()

        return ann

    # ----------------------------------------------------------------
    @classmethod
    def concat(cls, soma_slices):
        """
        Concatenates multiple `SOMASlice` objects into a single one. Implemented using `AnnData`'s
        `concat`, and in fact maybe `SOMASlice` should itself be simply an `AnnData` object.
        """

        if len(soma_slices) == 0:
            return None

        # Check column names for each dataframe-type are the same
        slice0 = soma_slices[0]
        for i, slicei in enumerate(soma_slices):
            if i == 0:
                continue
            # This works in Python -- not just a reference/pointer compare but a contents-compare :)
            assert list(slicei.X.keys()) == list(slice0.X.keys())
            assert list(slicei.obs.keys()) == list(slice0.obs.keys())
            assert list(slicei.var.keys()) == list(slice0.var.keys())

        result_X_df = pd.concat(
            [soma_slice.X for soma_slice in soma_slices], verify_integrity=True
        )
        result_obs_df = pd.concat(
            [soma_slice.obs for soma_slice in soma_slices], verify_integrity=True
        )
        result_var_df = pd.concat(
            [soma_slice.var for soma_slice in soma_slices], verify_integrity=True
        )

        return cls(
            X=result_X_df,
            obs=result_obs_df,
            var=result_var_df,
        )
