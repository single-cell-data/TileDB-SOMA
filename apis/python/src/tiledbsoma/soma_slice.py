from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Union

import anndata as ad
import pandas as pd
import pyarrow as pa

from tiledbsoma import util

from .tiledb_group import TileDBGroup
from .types import Matrix


class SOMASlice(TileDBGroup):
    """
    In-memory-only object for ephemeral extracting out of a SOMA. Can be used to *construct* a SOMA
    but is not a SOMA (which would entail out-of-memory storage).  This is simply a collection of
    either ``pandas.DataFrame`` or ``pyarrow.Table`` objects.
    """

    # ----------------------------------------------------------------
    def __init__(
        self,
        X: Dict[str, Union[pd.DataFrame, pa.Table, Matrix]],
        obs: Union[pd.DataFrame, pa.Table],
        var: Union[pd.DataFrame, pa.Table],
        # TODO
        # obsm: Dict[str, pd.DataFrame],
        # varm: Dict[str, pd.DataFrame],
        # obsp: Dict[str, pd.DataFrame],
        # varp: Dict[str, pd.DataFrame],
        # raw_X: Dict[str, pd.DataFrame],
        # raw_var: pd.DataFrame,
    ):
        """
        Constructs an in-memory ``SOMASlice`` object. This is a simple collection of obs, var, and X dataframes.
        """
        assert isinstance(obs, pd.DataFrame) or isinstance(obs, pa.Table)
        assert isinstance(var, pd.DataFrame) or isinstance(obs, pa.Table)

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
        Constructs an ``AnnData`` object from the current ``SOMASlice`` object.
        """

        # Find the dtype.
        X_data = self.X["data"]
        if isinstance(X_data, pd.DataFrame):
            X_dtype = X_data.dtypes["value"]
        elif isinstance(X_data, pa.Table):
            X_dtype = X_data["value"].type.to_pandas_dtype()
        else:
            X_dtype = X_data.dtype

        obs = self.obs
        var = self.var
        if isinstance(obs, pa.Table):
            obs = obs.to_pandas()
            obs.set_index("obs_id", inplace=True)
        if isinstance(var, pa.Table):
            var = var.to_pandas()
            var.set_index("var_id", inplace=True)
        ann = ad.AnnData(obs=obs, var=var, dtype=X_dtype)

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
                    obs.index,
                    var.index,
                )

            if isinstance(data, pa.Table):
                data = data.to_pandas()
                data.set_index(["obs_id", "var_id"], inplace=True)

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
                ann.X = data
            else:
                ann.layers[name] = data

        return ann

    # ----------------------------------------------------------------
    @classmethod
    def _keys_for_concat(cls, df: Union[pd.DataFrame, pa.Table, Matrix]) -> List[str]:
        if isinstance(df, pa.Table):
            return sorted(list(df.column_names))
        else:
            return sorted(list(df.keys()))

    @classmethod
    def concat(cls, soma_slices: Sequence[SOMASlice]) -> Optional[SOMASlice]:
        """
        Concatenates multiple ``SOMASlice`` objects into a single one. Implemented using ``AnnData``'s
        ``concat``. Requires that all slices share the same ``obs`` and ``var`` keys. Please
        see the ``SOMA`` class method ``find_common_obs_and_var_keys``.
        """

        if len(soma_slices) == 0:
            return None

        # Check column names for each dataframe-type are the same
        slice0 = soma_slices[0]
        okeys0 = cls._keys_for_concat(slice0.obs)
        vkeys0 = cls._keys_for_concat(slice0.var)

        for i, slicei in enumerate(soma_slices):
            if i == 0:
                continue
            # This list-equals works in Python -- not just a reference/pointer compare but a contents-compare :)
            if sorted(list(slicei.X.keys())) != sorted(list(slice0.X.keys())):
                raise Exception(
                    "SOMA slices to be concatenated must have all the same X attributes"
                )
            for key in slice0.X.keys():
                if cls._keys_for_concat(slicei.X[key]) != cls._keys_for_concat(
                    slice0.X[key]
                ):
                    raise Exception(
                        "SOMA slices to be concatenated must have all the same obs attributes"
                    )

            if cls._keys_for_concat(slicei.obs) != okeys0:
                raise Exception(
                    "SOMA slices to be concatenated must have all the same obs attributes"
                )
            if cls._keys_for_concat(slicei.var) != vkeys0:
                raise Exception(
                    "SOMA slices to be concatenated must have all the same var attributes"
                )

        # Use AnnData.concat.
        # TODO: try to remove this dependency.

        # This is a possible candidate for threading -- one can imagine to_anndata() being at least
        # partially implemented in C++ and releasing the GIL. However it appears to be pure Python;
        # experiments have shown that threading here actually increases execution time by a little
        # bit.

        anns = [soma_slice.to_anndata() for soma_slice in soma_slices]

        # @classmethod
        # def _concat_aux(cls, index, soma_slice):
        #     ann = soma_slice.to_anndata()
        #     return (index, ann)

        # futures = []
        # max_thread_pool_workers = 8
        # with ThreadPoolExecutor(max_workers=max_thread_pool_workers) as executor:
        #     for i, soma_slice in enumerate(soma_slices):
        #         future = executor.submit(
        #             cls._concat_aux,
        #             i,
        #             soma_slice,
        #         )
        #         futures.append(future)
        # anns = []
        # for future in futures:
        #     index, ann = future.result()
        #     if ann is not None:
        #         anns.append(ann)

        # We find that the ad.concat is relatively quick.
        annc = ad.concat(anns, join="outer", merge="first")

        # Problem: https://discourse.scverse.org/t/help-with-concat/676
        #
        # Example of the problem:
        # ---- SOMA1        ---- SOMA2        ---- SLICE (WITH QUERY ALL) AND CONCAT
        # X:                X:                X:
        # [[1. 2. 3.]       [[ 7.  8.]        [[ 1.  2.  3.    ]
        #  [4. 5. 6.]]       [ 9. 10.]]        [ 4.  5.  6.    ]
        # obs:              obs:               [ 7.          8.]
        #        oa ob             oa ob       [ 9.         10.]]
        # cell1  10  a      cell3  60  f      obs:
        # cell2  20  b      cell4  70  g             oa ob
        # var:              var:              cell1  10  a
        #        va vb             va vb      cell2  20  b
        # gene1  30  c      gene1  80  h      cell3  60  f
        # gene2  40  d      gene4  90  i      cell4  70  g
        # gene3  50  e                        var:
        #                                              va   vb
        #                                     gene1  30.0    c
        #                                     gene2  40.0    d
        #                                     gene3  50.0    e
        #                                     gene4   NaN  NaN <---- make sure this doesn't happen

        # Solution thanks to Bruce Martin at CZI:
        merged_obs = pd.concat([ann.obs for ann in anns], join="outer")
        merged_var = pd.concat([ann.var for ann in anns], join="outer")
        annc.obs = merged_obs[~merged_obs.index.duplicated()]
        annc.var = merged_var[~merged_var.index.duplicated()]

        annc.obs_names_make_unique()
        annc.var_names_make_unique()

        # Having leveraged AnnData solely for the concat method, pull out concatenated dataframes
        # into a SOMASlice object.
        X = {}
        X["data"] = annc.X
        for name in annc.layers:
            X[name] = annc.layers[name]

        retval = cls(X=X, obs=annc.obs, var=annc.var)
        return retval
