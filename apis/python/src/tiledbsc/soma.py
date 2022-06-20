import os
from typing import Optional, Union, List, Dict

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


class SOMA(TileDBGroup):
    """Single-cell group
    Class for representing a group of TileDB groups/arrays that constitute an SOMA ('stack of matrices, annotated')
    which includes:

    * `X` (`AssayMatrixGroup`): a group of one or more labeled 2D sparse arrays that share the same dimensions.
    * `obs` (`AnnotationDataframe`): 1D labeled array with column labels for `X`
    * `var` (`AnnotationDataframe`): 1D labeled array with row labels for `X`

    Convenience accessors include:

    * `soma.obs_keys()` for `soma.obs_names` for `soma.obs.ids()`
    * `soma.var_keys()` for `soma.var_names` for `soma.var.ids()`
    * `soma.n_obs` for `soma.obs.shape()[0]`
    * `soma.n_var` for `soma.var.shape()[0]`
    """

    X: AssayMatrixGroup
    obs: AnnotationDataFrame
    var: AnnotationDataFrame
    obsm: AnnotationMatrixGroup
    varm: AnnotationMatrixGroup
    obsp: AnnotationPairwiseMatrixGroup
    varp: AnnotationPairwiseMatrixGroup
    raw: RawGroup
    uns: UnsGroup

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name=None,
        soma_options: Optional[SOMAOptions] = None,
        verbose: Optional[bool] = True,
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ):
        """
        Create a new SOMA object. The existing array group is opened at the specified array `uri` if one is present, otherwise a new array group is created.

        :param uri: URI of the TileDB group
        :param verbose: Print status messages
        """

        if ctx is None and config is not None:
            ctx = tiledb.Ctx(config)
        if soma_options is None:
            soma_options = SOMAOptions()  # Use default values from the constructor
        if name is None:
            name = os.path.basename(uri.rstrip("/"))
            if name == "":
                name = "soma"
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            verbose=verbose,
            soma_options=soma_options,
            ctx=ctx,
        )

        obs_uri = os.path.join(self.uri, "obs")
        var_uri = os.path.join(self.uri, "var")
        X_uri = os.path.join(self.uri, "X")
        obsm_uri = os.path.join(self.uri, "obsm")
        varm_uri = os.path.join(self.uri, "varm")
        obsp_uri = os.path.join(self.uri, "obsp")
        varp_uri = os.path.join(self.uri, "varp")
        raw_uri = os.path.join(self.uri, "raw")
        uns_uri = os.path.join(self.uri, "uns")

        self.obs = AnnotationDataFrame(uri=obs_uri, name="obs", parent=self)
        self.var = AnnotationDataFrame(uri=var_uri, name="var", parent=self)
        self.X = AssayMatrixGroup(
            uri=X_uri,
            name="X",
            row_dim_name="obs_id",
            col_dim_name="var_id",
            row_dataframe=self.obs,
            col_dataframe=self.var,
            parent=self,
        )
        self.obsm = AnnotationMatrixGroup(uri=obsm_uri, name="obsm", parent=self)
        self.varm = AnnotationMatrixGroup(uri=varm_uri, name="varm", parent=self)
        self.obsp = AnnotationPairwiseMatrixGroup(
            uri=obsp_uri,
            name="obsp",
            row_dataframe=self.obs,
            col_dataframe=self.obs,
            parent=self,
        )
        self.varp = AnnotationPairwiseMatrixGroup(
            uri=varp_uri,
            name="varp",
            row_dataframe=self.var,
            col_dataframe=self.var,
            parent=self,
        )
        self.raw = RawGroup(uri=raw_uri, name="raw", obs=self.obs, parent=self)
        self.uns = UnsGroup(uri=uns_uri, name="uns", parent=self)

        # If URI is "/something/test1" then:
        # * obs_uri  is "/something/test1/obs"
        # * var_uri  is "/something/test1/var"
        # * data_uri is "/something/test1/X"

        # If URI is "tiledb://namespace/s3://bucketname/something/test1" then:
        # * obs_uri  is "tiledb://namespace/s3://bucketname/something/test1/obs"
        # * var_uri  is "tiledb://namespace/s3://bucketname/something/test1/var"
        # * data_uri is "tiledb://namespace/s3://bucketname/something/test1/X"

    # ----------------------------------------------------------------
    def __str__(self):
        """
        Implements `print(soma)`.
        """
        return f"name={self.name},uri={self.uri}"

    # ----------------------------------------------------------------
    def __getattr__(self, name):
        """
        This is called on `soma.name` when `name` is not already an attribute.
        This is used for `soma.n_obs`, etc.
        """
        if name == "n_obs":
            return self.obs.shape()[0]
        if name == "n_var":
            return self.var.shape()[0]

        if name == "obs_names":
            return self.obs.ids()
        if name == "var_names":
            return self.var.ids()

        raise AttributeError(
            f"'{self.__class__.__name__}' object has no attribute '{name}'"
        )

    # ----------------------------------------------------------------
    def obs_keys(self):
        """
        An alias for `soma.obs.ids()`.
        """
        return self.obs.ids()

    # ----------------------------------------------------------------
    def var_keys(self):
        """
        An alias for `soma.var.ids()`.
        """
        return self.var.ids()

    # ----------------------------------------------------------------
    def cell_count(self) -> int:
        """
        Returns the `obs_id` in `soma.obs`.
        """
        return len(self.obs.ids())

    # ----------------------------------------------------------------
    def get_obs_value_counts(self, obs_label: str) -> pd.DataFrame:
        """
        Given an obs label, e.g. `cell_type`, returns a dataframe count the number of different
        values for that label in the SOMA.
        """
        return self._get_obs_or_var_value_counts(obs_label, True)

    def get_var_value_counts(self, var_label: str) -> pd.DataFrame:
        """
        Given an var label, e.g. `feature_name`, returns a dataframe count the number of different
        values for that label in the SOMA.
        """
        return self._get_obs_or_var_value_counts(var_label, False)

    def _get_obs_or_var_value_counts(
        self, obs_or_var_label: str, use_obs: True
    ) -> pd.DataFrame:
        """
        Supporting method for `get_obs_value_counts` and `get_var_value_counts`.
        """

        attrs = [obs_or_var_label]
        obs_or_var = self.obs.df(attrs=attrs) if use_obs else self.var.df(attrs=attrs)
        if not obs_or_var_label in obs_or_var:
            return

        counts = {}
        obs_label_values = list(obs_or_var[obs_or_var_label])
        for obs_label_value in obs_label_values:
            if obs_label_value in counts:
                counts[obs_label_value] += 1
            else:
                counts[obs_label_value] = 1

        name_column = []
        counts_column = []
        for k, v in dict(
            sorted(counts.items(), reverse=True, key=lambda item: item[1])
        ).items():
            name_column.append(k)
            counts_column.append(v)

        df = pd.DataFrame.from_dict({"name": name_column, "count": counts_column})
        df.set_index("name", inplace=True)
        return df
