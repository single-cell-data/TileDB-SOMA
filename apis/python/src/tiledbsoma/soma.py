from __future__ import annotations

import os
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from typing import Any, List, Optional, Sequence, Tuple, Union

import pandas as pd
import pyarrow as pa
import tiledb

from .annotation_dataframe import AnnotationDataFrame
from .annotation_matrix_group import AnnotationMatrixGroup
from .annotation_pairwise_matrix_group import AnnotationPairwiseMatrixGroup
from .assay_matrix import AssayMatrix
from .assay_matrix_group import AssayMatrixGroup
from .raw_group import RawGroup
from .soma_options import SOMAOptions
from .soma_slice import SOMASlice
from .tiledb_group import TileDBGroup
from .types import Ids, Matrix
from .uns_group import UnsGroup
from .util import INGEST_MODES


class SOMA(TileDBGroup):
    """
    Class for representing a group of TileDB groups/arrays that constitute an SOMA ('stack of matrices, annotated')
    which includes:

    * ``X`` (group of ``AssayMatrixGroup``): a group of one or more labeled 2D sparse arrays that share the same dimensions.
    * ``obs`` (``AnnotationDataframe``): 1D labeled array with column labels for ``X``
    * ``var`` (``AnnotationDataframe``): 1D labeled array with row labels for ``X``
    * ``obsm`` (group of ``AnnotationMatrix``): multi-attribute arrays keyed by IDs of ``obs``
    * ``varm`` (group of ``AnnotationMatrix``): multi-attribute arrays keyed by IDs of ``var``
    * ``obsp`` (group of ``AnnotationMatrix``): 2D arrays keyed by IDs of ``obs``
    * ``varp`` (group of ``AnnotationMatrix``): 2D arrays keyed by IDs of ``var``
    * ``raw``: contains raw versions of ``X`` and ``varm``
    * ``uns``: nested, unstructured data

    Convenience accessors include:

    * ``soma.obs_keys()`` for ``soma.obs_names`` for ``soma.obs.ids()``
    * ``soma.var_keys()`` for ``soma.var_names`` for ``soma.var.ids()``
    * ``soma.n_obs`` for ``soma.obs.shape()[0]``
    * ``soma.n_var`` for ``soma.var.shape()[0]``
    """

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        soma_options: Optional[SOMAOptions] = None,
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ):
        """
        Create a new SOMA object. The existing array group is opened at the specified array ``uri`` if one is present, otherwise a new array group is created.

        :param uri: URI of the TileDB group
        """

        # People can (and should) call by name. However, it's easy to forget. For example,
        # if someone does 'tiledbsoma.SOMA("myuri", ctx)' instead of 'tiledbsoma.SOMA("myury", ctx)',
        # behavior will not be what they expect, and we should let them know sooner than later.
        if name is not None:
            assert isinstance(name, str)
        if soma_options is not None:
            assert isinstance(soma_options, SOMAOptions)
        if config is not None:
            assert isinstance(config, tiledb.Config)
        if ctx is not None:
            assert isinstance(ctx, tiledb.Ctx)
        if parent is not None:
            assert isinstance(parent, TileDBGroup)

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

        # See comments in _get_child_uris
        child_uris = self._get_child_uris(
            ["obs", "var", "X", "obsm", "varm", "obsp", "varp", "raw", "uns"]
        )

        obs_uri = child_uris["obs"]
        var_uri = child_uris["var"]
        X_uri = child_uris["X"]
        obsm_uri = child_uris["obsm"]
        varm_uri = child_uris["varm"]
        obsp_uri = child_uris["obsp"]
        varp_uri = child_uris["varp"]
        raw_uri = child_uris["raw"]
        uns_uri = child_uris["uns"]

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

        # Sample SOMA/member URIs:
        #
        # SOMA                             member
        # "/foo/test1"                     "/foo/test1/obs"
        # "tiledb://ns/s3://bkt/foo/test1" "tiledb://ns/s3://bkt/foo/test1/obs"
        # "tiledb://ns/test1"              "tiledb://ns/95cb11b5-2052-461b-9b99-cfa94e51e23f"

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
            if self.uns.exists():
                # repr(self.uns) is very chatty (too chatty) for some datasets:
                lines.append("uns:     " + ", ".join(self.uns.keys()))
        else:
            lines.append("Unpopulated")

        return "\n".join(lines)

    # ----------------------------------------------------------------
    @property
    def n_obs(self) -> int:
        return self.obs.shape()[0]

    @property
    def n_var(self) -> int:
        return self.var.shape()[0]

    @property
    def obs_names(self) -> Sequence[str]:
        return self.obs.ids()

    @property
    def var_names(self) -> Sequence[str]:
        return self.var.ids()

    # ----------------------------------------------------------------
    def obs_keys(self) -> Sequence[str]:
        """
        An alias for ``soma.obs.ids()``.
        """
        return self.obs.ids()

    # ----------------------------------------------------------------
    def var_keys(self) -> Sequence[str]:
        """
        An alias for ``soma.var.ids()``.
        """
        return self.var.ids()

    # ----------------------------------------------------------------
    def get_obs_value_counts(self, obs_label: str) -> pd.DataFrame:
        """
        Given an obs label, e.g. ``cell_type``, returns a dataframe count the number of different
        values for that label in the SOMA.
        """
        return self._get_obs_or_var_value_counts(obs_label, True)

    def get_var_value_counts(self, var_label: str) -> pd.DataFrame:
        """
        Given a var label, e.g. ``feature_name``, returns a dataframe count the number of different
        values for that label in the SOMA.
        """
        return self._get_obs_or_var_value_counts(var_label, False)

    def _get_obs_or_var_value_counts(
        self, obs_or_var_label: str, use_obs: bool
    ) -> pd.DataFrame:
        """
        Supporting method for ``get_obs_value_counts`` and ``get_var_value_counts``.
        """
        attrs = [obs_or_var_label]
        obs_or_var = self.obs.df(attrs=attrs) if use_obs else self.var.df(attrs=attrs)
        if obs_or_var_label not in obs_or_var:
            return

        counts = Counter(obs_or_var[obs_or_var_label])
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
    def dim_slice(
        self,
        obs_ids: Optional[Ids],
        var_ids: Optional[Ids],
        *,
        return_arrow: bool = False,
    ) -> Optional[SOMASlice]:
        """
        Subselects the SOMA's ``obs``, ``var``, and ``X/data`` using the specified ``obs_ids`` and ``var_ids``.
        Using a value of ``None`` for obs_ids means use all ``obs_ids``, and likewise for ``var_ids``.
        Returns ``None`` for empty slice.
        """
        assert obs_ids is not None or var_ids is not None
        if obs_ids is None:
            # Try the var slice first to see if that produces zero results -- if so we don't need to
            # load the obs.
            slice_var_df = self.var.dim_select(var_ids, return_arrow=return_arrow)
            if slice_var_df.shape[0] == 0:
                return None
            slice_obs_df = self.obs.dim_select(obs_ids, return_arrow=return_arrow)
            if slice_obs_df.shape[0] == 0:
                return None

        elif var_ids is None:
            # Try the obs slice first to see if that produces zero results -- if so we don't need to
            # load the var.
            slice_obs_df = self.obs.dim_select(obs_ids, return_arrow=return_arrow)
            if slice_obs_df.shape[0] == 0:
                return None
            slice_var_df = self.var.dim_select(var_ids, return_arrow=return_arrow)
            if slice_var_df.shape[0] == 0:
                return None

        else:
            slice_obs_df = self.obs.dim_select(obs_ids, return_arrow=return_arrow)
            if slice_obs_df.shape[0] == 0:
                return None
            slice_var_df = self.var.dim_select(var_ids, return_arrow=return_arrow)
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

        return self._assemble_soma_slice(
            obs_ids, var_ids, slice_obs_df, slice_var_df, return_arrow=return_arrow
        )

    # ----------------------------------------------------------------
    def query(
        self,
        *,
        obs_attrs: Optional[Sequence[str]] = None,
        obs_query_string: Optional[str] = None,
        obs_ids: Optional[Ids] = None,
        var_attrs: Optional[Sequence[str]] = None,
        var_query_string: Optional[str] = None,
        var_ids: Optional[Ids] = None,
        X_layer_names: Optional[Sequence[str]] = None,
        return_arrow: bool = False,
    ) -> Optional[SOMASlice]:
        """
        Subselects the SOMA's obs, var, and X/data using the specified queries on obs and var.
        Queries use the TileDB-Py ``QueryCondition`` API.

        If ``obs_query_string`` is ``None``, the ``obs`` dimension is not filtered and all of ``obs`` is
        used; similiarly for ``var``.

        If ``obs_attrs`` or ``var_attrs`` are unspecified, the slice will take all ``obs``/``var`` attributes
        from the source SOMAs; if they are specified, the slice will take the specified ``obs``/``var``

        If ``X_layer_names`` is `None`, they are all returned; otherwise you can specify which layer(s)
        you want to be operated on.
        """

        retval = self._query_aux(
            obs_attrs=obs_attrs,
            obs_query_string=obs_query_string,
            obs_ids=obs_ids,
            var_attrs=var_attrs,
            var_query_string=var_query_string,
            var_ids=var_ids,
            X_layer_names=X_layer_names,
            return_arrow=return_arrow,
        )
        return retval

    # ----------------------------------------------------------------
    def _query_aux(
        self,
        *,
        obs_attrs: Optional[Sequence[str]] = None,
        obs_query_string: Optional[str] = None,
        obs_ids: Optional[Ids] = None,
        var_attrs: Optional[Sequence[str]] = None,
        var_query_string: Optional[str] = None,
        var_ids: Optional[Ids] = None,
        X_layer_names: Optional[Sequence[str]] = None,
        return_arrow: bool = False,
    ) -> Optional[SOMASlice]:
        """
        Helper method for `query`: as this has multiple `return` statements, it's easiest to track
        elapsed-time stats in a call to this helper.
        """

        # This is a good candidate for parallelization because the soma.query bits are independent
        # of one another, and are also in TileDB's core C++ engine which releases the GIL.
        max_thread_pool_workers = SOMAOptions().max_thread_pool_workers
        id_df_futures = []
        with ThreadPoolExecutor(max_workers=max_thread_pool_workers) as executor:

            obs_future = executor.submit(
                self._query_aux_obs,
                obs_attrs=obs_attrs,
                obs_query_string=obs_query_string,
                obs_ids=obs_ids,
                return_arrow=return_arrow,
            )
            id_df_futures.append(obs_future)

            var_future = executor.submit(
                self._query_aux_var,
                var_attrs=var_attrs,
                var_query_string=var_query_string,
                var_ids=var_ids,
                return_arrow=return_arrow,
            )
            id_df_futures.append(var_future)

        for id_future in id_df_futures:
            triple = id_future.result()
            if triple is None:
                return None
            name, df, ids = triple
            if name == "obs":
                slice_obs_df = df
                obs_ids = ids
            elif name == "var":
                slice_var_df = df
                var_ids = ids
            else:
                # ToDO
                # raise SOMAError
                raise Exception("internal coding error")

        return self._assemble_soma_slice(
            obs_ids,
            var_ids,
            slice_obs_df,
            slice_var_df,
            X_layer_names=X_layer_names,
            return_arrow=return_arrow,
        )

    # ----------------------------------------------------------------
    def _query_aux_obs(
        self,
        *,
        obs_attrs: Optional[Sequence[str]] = None,
        obs_query_string: Optional[str] = None,
        obs_ids: Optional[Ids] = None,
        return_arrow: bool = False,
    ) -> Any:  # Optional[SOMASlice]:
        """
        TODO
        """

        slice_obs_df = self.obs.query(
            query_string=obs_query_string,
            ids=obs_ids,
            attrs=obs_attrs,
            return_arrow=return_arrow,
        )
        # E.g. querying for 'cell_type == "blood"' and this SOMA does have a cell_type column in its
        # obs, but no rows with cell_type == "blood".
        if slice_obs_df is None:
            return None
        if return_arrow:
            if len(slice_obs_df["obs_id"]) == 0:
                return None
        else:
            if len(slice_obs_df.index) == 0:
                return None
        # At the tiledb multi-index level, if we're say slicing on obs_ids but not var_ids,
        # we'll do ``A.df[obs_ids, :]``. We can't pass a ``:`` down the callstack to get there,
        # but we pass ``None`` instead.
        #
        # It's important to do this. Say for example the X matrix is nobs=1000 by nvar=2000,
        # and we have a query that has 158 obs_ids. At the tiledb multi-index level, doing
        # ``A.df[{158 obs ids}, {all 2000 var ids}]`` is non-performant while
        # ``A.df[{158 obs ids}, :]`` is performant.
        if obs_ids is not None or obs_query_string is not None:
            if return_arrow:
                obs_ids = [obs_id.as_py() for obs_id in slice_obs_df["obs_id"]]
            else:
                obs_ids = list(slice_obs_df.index)

        return ("obs", slice_obs_df, obs_ids)

    # ----------------------------------------------------------------
    def _query_aux_var(
        self,
        *,
        var_attrs: Optional[Sequence[str]] = None,
        var_query_string: Optional[str] = None,
        var_ids: Optional[Ids] = None,
        return_arrow: bool = False,
    ) -> Any:  # Optional[SOMASlice]:
        """
        TODO
        """

        slice_var_df = self.var.query(
            var_query_string, ids=var_ids, attrs=var_attrs, return_arrow=return_arrow
        )
        # E.g. querying for 'feature_name == "MT-CO3"' and this SOMA does have a feature_name column
        # in its var, but no rows with feature_name == "MT-CO3".
        if slice_var_df is None:
            return None
        if return_arrow:
            if len(slice_var_df["var_id"]) == 0:
                return None
        else:
            if len(slice_var_df.index) == 0:
                return None
        # See above comment re keeping obs_ids == None if that's what it came in as.
        if var_ids is not None or var_query_string is not None:
            if return_arrow:
                var_ids = [var_id.as_py() for var_id in slice_var_df["var_id"]]
            else:
                var_ids = list(slice_var_df.index)

        return ("var", slice_var_df, var_ids)

    # ----------------------------------------------------------------
    @classmethod
    def queries(
        cls,
        somas: Sequence[SOMA],
        *,
        obs_attrs: Optional[Sequence[str]] = None,
        obs_query_string: Optional[str] = None,
        obs_ids: Optional[Ids] = None,
        var_attrs: Optional[Sequence[str]] = None,
        var_query_string: Optional[str] = None,
        var_ids: Optional[Ids] = None,
        X_layer_names: Optional[Sequence[str]] = None,
        return_arrow: bool = False,
        max_thread_pool_workers: Optional[int] = None,
    ) -> List[SOMASlice]:
        """
        Subselects the obs, var, and X/data using the specified queries on obs and var,
        concatenating across SOMAs in the list.  Queries use the TileDB-Py ``QueryCondition``
        API.

        If ``obs_query_string`` is ``None``, the ``obs`` dimension is not filtered and all of ``obs`` is
        used; similiarly for ``var``. Return value of ``None`` indicates an empty slice.  If ``obs_ids``
        or ``var_ids`` are not ``None``, they are effectively ANDed into the query.  For example, you
        can pass in a known list of ``obs_ids``, then use ``obs_query_string`` to further restrict the
        query.

        If ``obs_attrs`` or ``var_attrs`` are unspecified, slices will take all ``obs``/``var`` attributes
        from their source SOMAs; if they are specified, slices will take the specified ``obs``/``var``
        attributes.  If all SOMAs in the collection have the same ``obs``/``var`` attributes, then you
        needn't specify these; if they don't, you must.

        If ``X_layer_names`` is `None`, they are all returned; otherwise you can specify which layer(s)
        you want to be operated on.
        """

        # This is a good candidate for parallelization because the soma.query bits are independent
        # of one another, and are also in TileDB's core C++ engine which releases the GIL.
        if max_thread_pool_workers is None:
            max_thread_pool_workers = SOMAOptions().max_thread_pool_workers
        soma_slice_futures = []
        with ThreadPoolExecutor(max_workers=max_thread_pool_workers) as executor:
            for soma in somas:
                # E.g. querying for 'cell_type == "blood"' but this SOMA doesn't have a cell_type
                # column in its obs at all.
                if obs_query_string is not None and not soma.obs.has_attr_names(
                    obs_attrs or []
                ):
                    continue
                # E.g. querying for 'feature_name == "MT-CO3"' but this SOMA doesn't have a feature_name
                # column in its var at all.
                if var_query_string is not None and not soma.var.has_attr_names(
                    var_attrs or []
                ):
                    continue

                soma_slice_future = executor.submit(
                    soma.query,
                    obs_attrs=obs_attrs,
                    var_attrs=var_attrs,
                    obs_query_string=obs_query_string,
                    var_query_string=var_query_string,
                    obs_ids=obs_ids,
                    var_ids=var_ids,
                    X_layer_names=X_layer_names,
                    return_arrow=return_arrow,
                )
                soma_slice_futures.append(soma_slice_future)

        soma_slices = []
        for soma_slice_future in soma_slice_futures:
            soma_slice = soma_slice_future.result()
            if soma_slice is not None:
                soma_slices.append(soma_slice)
        return soma_slices

    # ----------------------------------------------------------------
    def _assemble_soma_slice_aux(
        self,
        layer_name: str,
        X: AssayMatrix,
        obs_ids: Optional[Ids],
        var_ids: Optional[Ids],
        *,
        return_arrow: bool = False,
    ) -> Tuple[str, Union[pd.DataFrame, pa.Table]]:
        return (layer_name, X.dim_select(obs_ids, var_ids, return_arrow=return_arrow))

    def _assemble_soma_slice(
        self,
        obs_ids: Optional[Ids],
        var_ids: Optional[Ids],
        slice_obs_df: Union[pd.DataFrame, pa.Table],
        slice_var_df: Union[pd.DataFrame, pa.Table],
        *,
        return_arrow: bool = False,
        X_layer_names: Optional[Sequence[str]] = None,
    ) -> SOMASlice:
        """
        An internal method for constructing a ``SOMASlice`` object given query results.
        """
        # There aren't always multiple X layers, and if there aren't, this parallelization doesn't
        # help. But neither does it hurt. And if there are, that's good news, since the dim_select is
        # in TileDB's C++ core engine which releases the GIL.
        futures = []
        if X_layer_names is None:
            X_layer_names = self.X.keys()
        with ThreadPoolExecutor(
            max_workers=self._soma_options.max_thread_pool_workers
        ) as executor:
            for X_layer_name in X_layer_names:
                X_layer = self.X[X_layer_name]
                if X_layer is not None:
                    future = executor.submit(
                        self._assemble_soma_slice_aux,
                        X_layer_name,
                        X_layer,
                        obs_ids,
                        var_ids,
                        return_arrow=return_arrow,
                    )
                    futures.append(future)

        X = {}
        for future in futures:
            X_layer_name, df = future.result()
            assert df is not None
            X[X_layer_name] = df

        return SOMASlice(X=X, obs=slice_obs_df, var=slice_var_df)

    # ----------------------------------------------------------------
    @classmethod
    def from_soma_slice(
        cls,
        soma_slice: SOMASlice,
        uri: str,
        name: Optional[str] = None,
        soma_options: Optional[SOMAOptions] = None,
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ) -> SOMA:
        """
        Constructs ``SOMA`` storage from a given in-memory ``SOMASlice`` object.
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
        for layer_name in soma_slice.X.keys():
            soma.X.add_layer_from_matrix_and_dim_values(
                matrix=soma_slice.X[layer_name],
                row_names=soma.obs.ids(),
                col_names=soma.var.ids(),
                layer_name=layer_name,
            )

        return soma

    # ----------------------------------------------------------------
    @classmethod
    def find_common_obs_and_var_keys(
        cls,
        somas: Sequence[SOMA],
    ) -> Tuple[List[str], List[str]]:
        obs_attrs_set = None
        var_attrs_set = None

        if len(somas) == 0:
            return ([], [])

        obs_attrs_set = set(somas[0].obs.keys())
        var_attrs_set = set(somas[0].var.keys())

        for soma in somas[1:]:
            obs_attrs_set = set(soma.obs.keys()).intersection(obs_attrs_set)
            var_attrs_set = set(soma.var.keys()).intersection(var_attrs_set)

        obs_attrs = sorted(list(obs_attrs_set))
        var_attrs = sorted(list(var_attrs_set))

        return (obs_attrs, var_attrs)

    # ----------------------------------------------------------------
    def add_X_layer(
        self,
        matrix: Matrix,
        layer_name: str = "data",
        ingest_mode: str = "write",
    ) -> None:
        """
        Populates the ``X`` or ``raw.X`` subgroup for a ``SOMA`` object.
        """
        assert ingest_mode in INGEST_MODES

        self.X.add_layer_from_matrix_and_dim_values(
            matrix,
            self.obs.ids(),
            self.var.ids(),
            layer_name,
            ingest_mode=ingest_mode,
        )
