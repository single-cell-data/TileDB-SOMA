from typing import Dict

import anndata as ad
import pandas as pd

from tiledbsc import util

from .tiledb_group import TileDBGroup


class SOMASlice(TileDBGroup):
    """
    In-memory-only object for ephemeral extracting out of a SOMA. Can be used to _construct_ a SOMA
    but is not a SOMA (which would entail out-of-memory storage).  Nothing more than a collection of
    pandas.DataFrame objects. No raw or uns.
    """

    X: pd.DataFrame
    obs: pd.DataFrame
    var: pd.DataFrame
    # TODO
    # obsm: Dict[str, pd.DataFrame]
    # varm: Dict[str, pd.DataFrame]
    # obsp: Dict[str, pd.DataFrame]
    # varp: Dict[str, pd.DataFrame]

    # ----------------------------------------------------------------
    def __init__(
        self,
        X: Dict[str, pd.DataFrame],
        obs: pd.DataFrame,
        var: pd.DataFrame,
        # TODO
        # obsm: Dict[str, pd.DataFrame],
        # varm: Dict[str, pd.DataFrame],
        # obsp: Dict[str, pd.DataFrame],
        # varp: Dict[str, pd.DataFrame],
        # raw_X: Dict[str, pd.DataFrame],
        # raw_var: pd.DataFrame,
    ):
        """
        Constructs an in-memory `SOMASlice` object. This is a simple collection of obs, var, and X dataframes.
        """
        assert isinstance(obs, pd.DataFrame)
        assert isinstance(var, pd.DataFrame)
        assert "data" in X

        self.obs = obs
        self.var = var
        self.X = X
        # TODO
        # self.obsm = obsm
        # self.varm = varm
        # self.obsp = obsp
        # self.varp = varp
        # self.raw_X = raw_X
        # self.raw_var = raw_var
        assert "data" in X

        # Find the dtype.
        X_data = self.X["data"]
        if isinstance(X_data, pd.DataFrame):
            X_dtype = X_data.dtypes["value"]
        else:
            X_dtype = X_data.dtype

        ann = ad.AnnData(obs=self.obs, var=self.var, dtype=X_dtype)

        for name, data in self.X.items():
            # X comes in from TileDB queries as a 3-column dataframe with "obs_id", "var_id", and
            # "value".  For AnnData we need to make it a sparse matrix.
            if isinstance(data, pd.DataFrame):
                # Make obs_id and var_id accessible as columns.
                data = data.reset_index()
                data = util.X_and_ids_to_sparse_matrix(
                    data,
                    "obs_id",  # row_dim_name
                    "var_id",  # col_dim_name
                    "value",  # attr_name
                    self.obs.index,
                    self.var.index,
                )
            # We use AnnData as our in-memory storage. For SOMAs, all X layers are arrays within the
            # soma.X group; for AnnData, the 'data' layer is ann.X and all the others are in
            # ann.layers.
            if name == "data":
                ann.X = data
            else:
                ann.layers[name] = data

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display of SOMASlice.
        """

        lines = []
        for key in self.X.keys():
            lines.append(f"X/{key}: {self.X[key].shape}")
        lines.append(f"obs: {self.obs.shape}")
        lines.append(f"var: {self.var.shape}")

        return "\n".join(lines)

    # ----------------------------------------------------------------
    def to_anndata(self) -> ad.AnnData:
        """
        Constructs an `AnnData` object from the current `SOMASlice` object.
        """

        # Find the dtype.
        X_data = self.X["data"]
        if isinstance(X_data, pd.DataFrame):
            X_dtype = X_data.dtypes["value"]
        else:
            X_dtype = X_data.dtype

        ann = ad.AnnData(obs=self.obs, var=self.var, dtype=X_dtype)

        # TODO:
        # self.obsm = obsm
        # self.varm = varm
        # self.obsp = obsp
        # self.varp = varp
        # self.raw_X = raw_X
        # self.raw_var = raw_var

        for name, data in self.X.items():
            # X comes in from TileDB queries as a 3-column dataframe with "obs_id", "var_id", and
            # "value".  For AnnData we need to make it a sparse matrix.
            if isinstance(data, pd.DataFrame):
                # Make obs_id and var_id accessible as columns.
                data = data.reset_index()
                data = util.X_and_ids_to_sparse_matrix(
                    data,
                    "obs_id",  # row_dim_name
                    "var_id",  # col_dim_name
                    "value",  # attr_name
                    self.obs.index,
                    self.var.index,
                )
            # We use AnnData as our in-memory storage. For SOMAs, all X layers are arrays within the
            # soma.X group; for AnnData, the 'data' layer is ann.X and all the others are in
            # ann.layers.
            if name == "data":
                ann.X = data
            else:
                ann.layers[name] = data

        return ann

    # ----------------------------------------------------------------
    @classmethod
    def concat(cls, soma_slices):
        """
        Concatenates multiple `SOMASlice` objects into a single one. Implemented using `AnnData`'s
        `concat`.
        """

        if len(soma_slices) == 0:
            return None

        # Check column names for each dataframe-type are the same
        slice0 = soma_slices[0]
        for i, slicei in enumerate(soma_slices):
            if i == 0:
                continue
            # This list-equals works in Python -- not just a reference/pointer compare but a contents-compare :)
            if sorted(list(slicei.X.keys())) != sorted(list(slice0.X.keys())):
                raise Exception(
                    "SOMA slices to be concatenated must have all the same X attributes"
                )
            if sorted(list(slicei.obs.keys())) != sorted(list(slice0.obs.keys())):
                raise Exception(
                    "SOMA slices to be concatenated must have all the same obs attributes"
                )
            if sorted(list(slicei.var.keys())) != sorted(list(slice0.var.keys())):
                raise Exception(
                    "SOMA slices to be concatenated must have all the same var attributes"
                )

        # Use AnnData.concat.
        # TODO: try to remove this dependency.
        anns = [soma_slice.to_anndata() for soma_slice in soma_slices]
        annc = ad.concat(anns, join="outer", merge="first")
        annc.obs_names_make_unique()
        annc.var_names_make_unique()

        # Having leveraged AnnData solely for the concat method, pull out concatenated dataframes
        # into a SOMASlice object.
        X = {}
        # TODO: SHAPE THIS
        X["data"] = annc.X
        for name in annc.layers:
            X[name] = annc.layers[name]

        return cls(
            X=X,
            obs=annc.obs,
            var=annc.var,
        )
