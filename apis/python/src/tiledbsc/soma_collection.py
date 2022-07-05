from typing import Iterator, Optional, Sequence, Set, Union

import tiledb

from .soma import SOMA
from .soma_options import SOMAOptions
from .soma_slice import SOMASlice
from .tiledb_group import TileDBGroup
from .types import Ids


class SOMACollection(TileDBGroup):
    """
    Implements a collection of `SOMA` objects.
    """

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        *,
        name: str = "soco",
        soma_options: Optional[SOMAOptions] = None,
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ):
        """
        Create a new `SOMACollection` object. The existing group is opened at the
        specified `uri` if one is present, otherwise a new group will be created upon ingest.

        :param uri: URI of the TileDB group
        """

        # People can (and should) call by name. However, it's easy to forget. For example,
        # if someone does 'tiledbsc.SOMACollection("myuri", ctx)' instead of 'tiledbsc.SOMA("myury", ctx)',
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
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            soma_options=soma_options,
            ctx=ctx,
        )

    # ----------------------------------------------------------------
    def __repr__(self) -> str:
        """
        Default display of SOMACollection.
        """
        lines = [
            "URI:        " + self.uri,
        ]
        if self.exists():
            lines.append(f"SOMA count: {len(self)}")
        else:
            lines.append("Unpopulated")
        return "\n".join(lines)

    # ----------------------------------------------------------------
    def __len__(self) -> int:
        """
        Implements `len(soco)`. Returns the number of elements in the collection.
        """
        return len(self._get_member_names())

    # ----------------------------------------------------------------
    def add(self, soma: SOMA) -> None:
        """
        Adds a `SOMA` to the `SOMACollection`.
        """
        self._add_object(soma)

    # ----------------------------------------------------------------
    def remove(self, soma: Union[SOMA, str]) -> None:
        """
        Removes a `SOMA` from the `SOMACollection`, when invoked as `soco.remove("namegoeshere")`.
        """
        if isinstance(soma, str):
            self._remove_object_by_name(soma)
        else:
            self._remove_object(soma)

    def __delattr__(self, matrix_name: str) -> None:
        """
        Removes a `SOMA` from the `SOMACollection`, when invoked as `del soco.namegoeshere`.
        """
        self.remove(matrix_name)

    def __delitem__(self, matrix_name: str) -> None:
        """
        Removes a `SOMA` from the `SOMACollection`, when invoked as `del soco["namegoeshere"]`.
        """
        self.remove(matrix_name)

    # ----------------------------------------------------------------
    def keys(self) -> Sequence[str]:
        """
        Returns the names of the SOMAs in the collection.
        """
        return self._get_member_names()

    # ----------------------------------------------------------------
    def __iter__(self) -> Iterator[SOMA]:
        """
        Implements `for soma in soco: ...`
        """
        for name, uri in self._get_member_names_to_uris().items():
            yield SOMA(uri=uri, name=name, parent=self, ctx=self._ctx)

    # ----------------------------------------------------------------
    def __contains__(self, name: str) -> bool:
        """
        Implements `name in soco`
        """
        with self._open("r") as G:
            return name in G

    # At the tiledb-py API level, *all* groups are name-indexable.  But here at the tiledbsc-py
    # level, we implement name-indexing only for some groups:
    #
    # * Most soma member references are done using Python's dot syntax. For example, rather than
    #   soma['X'], we have simply soma.X, and likewise, soma.raw.X.  Likewise soma.obs and soma.var.
    #
    # * Index references are supported for obsm, varm, obsp, varp, and uns. E.g.
    #   soma.obsm['X_pca'] or soma.uns['neighbors']['params']['method']
    #
    # * Overloading the `[]` operator at the TileDBGroup level isn't necessary -- e.g. we don't need
    #   soma['X'] when we have soma.X -- but also it causes circular-import issues in Python.
    #
    # * Rather than doing a TileDBIndexableGroup which overloads the `[]` operator, we overload
    #   the `[]` operator separately in the various classes which need indexing. This is again to
    #   avoid circular-import issues, and means that [] on `AnnotationMatrixGroup` will return an
    #   `AnnotationMatrix, [] on `UnsGroup` will return `UnsArray` or `UnsGroup`, etc.
    def __getitem__(self, name: str) -> Optional[SOMA]:
        """
        Returns a `SOMA` element at the given name within the group, or `None` if no such
        member exists.  Overloads the `[...]` operator.
        """
        with self._open("r") as G:
            try:
                obj = G[name]  # This returns a tiledb.object.Object.
            except tiledb.TileDBError:
                return None

            if obj.type != tiledb.group.Group:
                raise Exception(
                    f"Internal error: found element which is not a subgroup: type is {obj.type}"
                )

            return SOMA(uri=obj.uri, name=name, parent=self)

    # ----------------------------------------------------------------
    def query(
        self,
        *,
        obs_attrs: Optional[Sequence[str]] = None,
        obs_query_string: str = None,
        obs_ids: Optional[Ids] = None,
        var_attrs: Optional[Sequence[str]] = None,
        var_query_string: str = None,
        var_ids: Optional[Ids] = None,
    ) -> Optional[SOMASlice]:
        """
        Subselects the obs, var, and X/data using the specified queries on obs and var,
        concatenating across SOMAs in the collection.  Queries use the TileDB-Py `QueryCondition`
        API.

        If `obs_query_string` is `None`, the `obs` dimension is not filtered and all of `obs` is
        used; similiarly for `var`. Return value of `None` indicates an empty slice.  If `obs_ids`
        or `var_ids` are not `None`, they are effectively ANDed into the query.  For example, you
        can pass in a known list of `obs_ids`, then use `obs_query_string` to further restrict the
        query.

        If `obs_attrs` or `var_attrs` are unspecified, slices will take all `obs`/`var` attributes
        from their source SOMAs; if they are specified, slices will take the specified `obs`/`var`
        attributes.  If all SOMAs in the collection have the same `obs`/`var` attributes, then you
        needn't specify these; if they don't, you must.
        """

        soma_slices = []
        for soma in self:
            # E.g. querying for 'cell_type == "blood"' but this SOMA doesn't have a cell_type column in
            # its obs at all.
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

            soma_slice = soma.query(
                obs_attrs=obs_attrs,
                var_attrs=var_attrs,
                obs_query_string=obs_query_string,
                var_query_string=var_query_string,
                obs_ids=obs_ids,
                var_ids=var_ids,
            )
            if soma_slice is not None:
                # print("Slice SOMA from", soma.name, soma.X.data.shape(), "to", soma_slice.ann.X.shape)
                soma_slices.append(soma_slice)

        return SOMASlice.concat(soma_slices)

    # ----------------------------------------------------------------
    def find_unique_obs_values(self, obs_label: str) -> Set:
        """
        Given an `obs` label such as `cell_type` or `tissue`, returns a set of unique
        values for that label among all SOMAs in the collection.
        """
        return self._find_unique_obs_or_var_values(obs_label, True)

    def find_unique_var_values(self, var_label: str) -> Set:
        """
        Given an `var` label such as `feature_name`, returns a set of unique values for
        that label among all SOMAs in the collection.
        """
        return self._find_unique_obs_or_var_values(var_label, False)

    def _find_unique_obs_or_var_values(
        self, obs_or_var_label: str, use_obs: bool
    ) -> Set:
        """
        Helper method for `find_unique_obs_values` and `find_unique_var_values`.
        """
        unique_values_in_soco = set()

        for soma in self:
            annotation_matrix = soma.obs if use_obs else soma.var
            if obs_or_var_label not in annotation_matrix.keys():
                continue
            unique_values_in_soco.update(annotation_matrix.df()[obs_or_var_label])

        return unique_values_in_soco
