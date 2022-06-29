import os
from typing import Dict, List, Optional

import pandas as pd
import tiledb

from .annotation_dataframe import AnnotationDataFrame
from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
from .assay_matrix_group import AssayMatrixGroup
from .raw_group import RawGroup
from .soma_options import SOMAOptions
from .soma_slice import SOMASlice
from .tiledb_group import TileDBGroup
from .uns_group import UnsGroup


class SOMA(TileDBGroup):
    """Single-cell group
    Class for representing a group of TileDB groups/arrays that constitute an SOMA ('stack of matrices, annotated')
    which includes:

    * `X` (group of `AssayMatrixGroup`): a group of one or more labeled 2D sparse arrays that share the same dimensions.
    * `obs` (`AnnotationDataframe`): 1D labeled array with column labels for `X`
    * `var` (`AnnotationDataframe`): 1D labeled array with row labels for `X`
    * `obsm` (group of `AnnotationMatrix`): multi-attribute arrays keyed by IDs of `obs`
    * `varm` (group of `AnnotationMatrix`): multi-attribute arrays keyed by IDs of `var`
    * `obsp` (group of `AnnotationMatrix`): 2D arrays keyed by IDs of `obs`
    * `varp` (group of `AnnotationMatrix`): 2D arrays keyed by IDs of `var`
    * `raw`: contains raw versions of `X` and `varm`
    * `uns`: nested, unstructured data

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
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ):
        """
        Create a new SOMA object. The existing array group is opened at the specified array `uri` if one is present, otherwise a new array group is created.

        :param uri: URI of the TileDB group
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
            soma_options=soma_options,
            ctx=ctx,
        )

        obs_uri = self._get_child_uri("obs")  # See comments in that function
        var_uri = self._get_child_uri("var")
        X_uri = self._get_child_uri("X")
        obsm_uri = self._get_child_uri("obsm")
        varm_uri = self._get_child_uri("varm")
        obsp_uri = self._get_child_uri("obsp")
        varp_uri = self._get_child_uri("varp")
        raw_uri = self._get_child_uri("raw")
        uns_uri = self._get_child_uri("uns")

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
    def __repr__(self) -> str:
        """
        Default display of SOMA.
        """

        lines = [
            "Name:    " + self.name,
            "URI:     " + self.uri,
        ]
        if self.exists():
            lines.append(f"(n_obs, n_var): ({len(self.obs)}, {len(self.var)})")
            lines.append("X:       " + repr(self.X))
            lines.append("obs:     " + repr(self.obs))
            lines.append("var:     " + repr(self.var))
            lines.append("obsm:    " + repr(self.obsm))
            lines.append("varm:    " + repr(self.varm))
            lines.append("obsp:    " + repr(self.obsp))
            lines.append("varp:    " + repr(self.varp))
            if self.raw.exists():
                lines.append("raw/X:   " + repr(self.raw.X))
                lines.append("raw/var: " + repr(self.raw.var))
            # repr(self.uns) is very chatty (too chatty) for some datasets:
            lines.append("uns:     " + ", ".join(self.uns.keys()))
        else:
            lines.append("Unpopulated")

        return "\n".join(lines)

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
    def obs_keys(self) -> List[str]:
        """
        An alias for `soma.obs.ids()`.
        """
        return self.obs.ids()

    # ----------------------------------------------------------------
    def var_keys(self) -> List[str]:
        """
        An alias for `soma.var.ids()`.
        """
        return self.var.ids()

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
        if obs_or_var_label not in obs_or_var:
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

    # ----------------------------------------------------------------
    def dim_slice(self, obs_ids, var_ids) -> Dict:
        """
        Subselects the SOMA's obs, var, and X/data using the specified obs_ids and var_ids.
        Using a value of `None` for obs_ids means use all obs_ids, and likewise for var_ids.
        Returns `None` for empty slice.
        """
        assert obs_ids is not None or var_ids is not None
        if obs_ids is None:
            # Try the var slice first to see if that produces zero results -- if so we don't need to
            # load the obs.
            slice_var_df = self.var.dim_select(var_ids)
            if slice_var_df.shape[0] == 0:
                return None
            slice_obs_df = self.obs.dim_select(obs_ids)
            if slice_obs_df.shape[0] == 0:
                return None

        elif var_ids is None:
            # Try the obs slice first to see if that produces zero results -- if so we don't need to
            # load the var.
            slice_obs_df = self.obs.dim_select(obs_ids)
            if slice_obs_df.shape[0] == 0:
                return None
            slice_var_df = self.var.dim_select(var_ids)
            if slice_var_df.shape[0] == 0:
                return None

        else:
            slice_obs_df = self.obs.dim_select(obs_ids)
            if slice_obs_df.shape[0] == 0:
                return None
            slice_var_df = self.var.dim_select(var_ids)
            if slice_var_df.shape[0] == 0:
                return None

        # TODO:
        # do this here:
        # * raw_var
        # do these in _assemble_soma_slice:
        # * raw_X
        # * obsm
        # * varm
        # * obsp
        # * varp

        return self._assemble_soma_slice(obs_ids, var_ids, slice_obs_df, slice_var_df)

    # ----------------------------------------------------------------
    def query(
        self,
        obs_attrs: Optional[List[str]] = None,
        obs_query_string: Optional[str] = None,
        obs_ids: Optional[List[str]] = None,
        var_attrs: Optional[List[str]] = None,
        var_query_string: Optional[str] = None,
        var_ids: Optional[List[str]] = None,
    ) -> SOMASlice:
        """
        Subselects the SOMA's obs, var, and X/data using the specified queries on obs and var.
        Queries use the TileDB-Py `QueryCondition` API.

        If `obs_query_string` is `None`, the `obs` dimension is not filtered and all of `obs` is
        used; similiarly for `var`.

        If `obs_attrs` or `var_attrs` are unspecified, the slice will take all `obs`/`var` attributes
        from the source SOMAs; if they are specified, the slice will take the specified `obs`/`var`
        """

        slice_obs_df = self.obs.query(
            query_string=obs_query_string, ids=obs_ids, attrs=obs_attrs
        )
        # E.g. querying for 'cell_type == "blood"' and this SOMA does have a cell_type column in its
        # obs, but no rows with cell_type == "blood".
        if slice_obs_df is None:
            return None
        obs_ids = list(slice_obs_df.index)
        if len(obs_ids) == 0:
            return None

        slice_var_df = self.var.query(
            query_string=var_query_string, ids=var_ids, attrs=var_attrs
        )
        # E.g. querying for 'feature_name == "MT-CO3"' and this SOMA does have a feature_name column
        # in its var, but no rows with feature_name == "MT-CO3".
        if slice_var_df is None:
            return None
        var_ids = list(slice_var_df.index)
        if len(var_ids) == 0:
            return None

        # TODO:
        # do this here:
        # * raw_var
        # do these in _assemble_soma_slice:
        # * raw_X
        # * obsm
        # * varm
        # * obsp
        # * varp

        return self._assemble_soma_slice(obs_ids, var_ids, slice_obs_df, slice_var_df)

    # ----------------------------------------------------------------
    def _assemble_soma_slice(
        self,
        obs_ids,
        var_ids,
        slice_obs_df,
        slice_var_df,
    ) -> SOMASlice:
        """
        An internal method for constructing a `SOMASlice` object given query results.
        """

        X = {key: self.X[key].dim_select(obs_ids, var_ids) for key in self.X.keys()}

        return SOMASlice(
            X=X,
            obs=slice_obs_df,
            var=slice_var_df,
        )

    # ----------------------------------------------------------------
    @classmethod
    def from_soma_slice(
        cls,
        soma_slice: SOMASlice,
        uri: str,
        name=None,
        soma_options: Optional[SOMAOptions] = None,
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ) -> None:
        """
        Constructs `SOMA` storage from a given in-memory `SOMASlice` object.
        """

        soma = cls(
            uri=uri,
            name=name,
            soma_options=soma_options,
            config=config,
            ctx=ctx,
            parent=parent,
        )

        soma.create_unless_exists()
        soma.obs.from_dataframe(soma_slice.obs)
        soma.var.from_dataframe(soma_slice.var)
        for name in soma_slice.X.keys():
            soma.X.add_layer_from_matrix_and_dim_values(
                soma_slice.X[name],
                soma.obs.ids(),
                soma.var.ids(),
                layer_name=name,
            )

        return soma
