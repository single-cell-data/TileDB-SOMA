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
    In-memory-only object for ephemeral extraction out of a SOMA. Can be used to _construct_ a SOMA
    but is not a SOMA (which would entail out-of-memory storage).  Nothing more than a collection of
    pandas.DataFrame objects. Currently implemented as an `AnnData` object, in order to reuse the
    non-trivial `concat` logic.
    """

    ann: ad.AnnData

    # ----------------------------------------------------------------
    def __init__(
        self,
        X_layer_data: Dict[str, pd.DataFrame],
        obs: pd.DataFrame,
        var: pd.DataFrame,
    ):
        """
        Constructs an in-memory `SOMASlice` object. This is a simple collection of obs, var, and X dataframes.
        """
        assert isinstance(obs, pd.DataFrame)
        assert isinstance(var, pd.DataFrame)
        assert "data" in X_layer_data

        # Find the dtype.
        X_data = X_layer_data["data"]
        if isinstance(X_data, pd.DataFrame):
            X_dtype = X_data.dtypes["value"]
        else:
            X_dtype = X_data.dtype

        self.ann = ad.AnnData(obs=obs, var=var, dtype=X_dtype)

        for name, data in X_layer_data.items():
            # X comes in as a 3-column dataframe with "obs_id", "var_id", and "value".
            # For AnnData we need to make it a sparse matrix.
            if isinstance(data, pd.DataFrame):
                # Make obs_id and var_id accessible as columns.
                data = data.reset_index()
                data = util.X_and_ids_to_sparse_matrix(
                    data,
                    "obs_id",  # row_dim_name
                    "var_id",  # col_dim_name
                    "value",  # attr_name
                    obs.index,
                    var.index,
                )
            # We use AnnData as our in-memory storage. For SOMAs, all X layers are arrays within the
            # soma.X group; for AnnData, the 'data' layer is ann.X and all the others are in
            # ann.layers.
            if name == "data":
                self.ann.X = data
            else:
                self.ann.layers[name] = data

    # ----------------------------------------------------------------
    def __getattr__(self, name):
        """
        Accessors for `obs`, `var`, and `X` attributes of `SOMASlice`.
        """
        if name == "obs":
            return self.ann.obs
        if name == "var":
            return self.ann.var
        if name == "X":
            return self.ann.X
        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    # ----------------------------------------------------------------
    def to_anndata(self) -> ad.AnnData:
        """
        Returns an `AnnData` object from the current `SOMASlice` object.
        """
        return self.ann

    # ----------------------------------------------------------------
    @classmethod
    def concat(cls, soma_slices):
        """
        Concatenates multiple `SOMASlice` objects into a single one. Implemented using `AnnData`'s
        `concat`, and in fact maybe `SOMASlice` should itself be simply an `AnnData` object.
        """

        if len(soma_slices) == 0:
            return None
        anns = [soma_slice.ann for soma_slice in soma_slices]
        annc = ad.concat(anns)
        annc.obs_names_make_unique()
        annc.var_names_make_unique()

        # We use AnnData as our in-memory storage. For SOMAs, all X layers are arrays within the
        # soma.X group; for AnnData, the 'data' layer is ann.X and all the others are in
        # ann.layers.
        X_layer_data = {"data": annc.X}
        for key in annc.layers:
            X_layer_data[key] = annc.layers[key]

        return SOMASlice(X_layer_data=X_layer_data, obs=annc.obs, var=annc.var)
