from typing import Optional, List

from .soma_options import SOMAOptions
from .soma import SOMA
from .tiledb_group import TileDBGroup

import tiledb
import pandas as pd


class SOMACollection(TileDBGroup):
    """
    Implements a collection of `SOMA` objects.
    """

    # ----------------------------------------------------------------
    def __init__(
        self,
        uri: str,
        name="soco",
        soma_options: Optional[SOMAOptions] = None,
        verbose: Optional[bool] = True,
        config: Optional[tiledb.Config] = None,
        ctx: Optional[tiledb.Ctx] = None,
        parent: Optional[TileDBGroup] = None,  # E.g. a SOMA collection
    ):
        """
        Create a new `SOMACollection` object. The existing group is opened at the specified `uri` if one is present, otherwise a new group will be created upon ingest.

        :param uri: URI of the TileDB group
        :param verbose: Print status messages
        """

        if ctx is None and config is not None:
            ctx = tiledb.Ctx(config)
        if soma_options is None:
            soma_options = SOMAOptions()  # Use default values from the constructor
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            verbose=verbose,
            soma_options=soma_options,
            ctx=ctx,
        )

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
    def remove(self, soma: SOMA) -> None:
        """
        Removes a `SOMA` from the `SOMACollection`.
        """
        self._remove_object(soma)

    # ----------------------------------------------------------------
    def keys(self) -> None:
        """
        Returns the names of the SOMAs in the collection.
        """
        return self._get_member_names()

    # ----------------------------------------------------------------
    def __iter__(self) -> List[SOMA]:
        """
        Implements `for soma in soco: ...`
        """
        retval = []
        for name, uri in self._get_member_names_to_uris().items():
            soma = SOMA(uri=uri, name=name, parent=self, ctx=self._ctx)
            retval.append(soma)
        return iter(retval)

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
    def __getitem__(self, name):
        """
        Returns a `SOMA` element at the given name within the group, or `None` if no such
        member exists.  Overloads the `[...]` operator.
        """

        with self._open("r") as G:
            try:
                obj = G[name]  # This returns a tiledb.object.Object.
            except:
                return None

            if obj.type != tiledb.group.Group:
                raise Exception(
                    f"Internal error: found element which is not a subgroup: type is {str(obj.type)}"
                )

            return SOMA(uri=obj.uri, name=name, parent=self)

    # ----------------------------------------------------------------
    def cell_count(self) -> int:
        """
        Returns sum of `soma.cell_count()` over SOMAs in the collection.
        """
        return sum(soma.cell_count() for soma in self)

    # ----------------------------------------------------------------
    def find_unique_obs_values(self, obs_label: str):
        """
        Given an `obs` label such as `cell_type` or `tissue`, returns a list of unique values for
        that label among all SOMAs in the collection.
        """
        return self._find_unique_obs_or_var_values(obs_label, True)

    def find_unique_var_values(self, var_label: str):
        """
        Given an `var` label such as `feature_name`, returns a list of unique values for
        that label among all SOMAs in the collection.
        """
        return self._find_unique_obs_or_var_values(var_label, False)

    def _find_unique_obs_or_var_values(self, obs_or_var_label: str, use_obs: bool):
        """
        Helper method for `find_unique_obs_values` and `find_unique_var_values`.
        """
        unique_values_in_soco = set()

        for soma in self:
            annotation_matrix = soma.obs if use_obs else soma.var
            if not obs_or_var_label in annotation_matrix.keys():
                continue

            unique_values_in_soma = list(set(annotation_matrix.df()[obs_or_var_label]))

            unique_values_in_soco = unique_values_in_soco.union(unique_values_in_soma)

        return unique_values_in_soco

    # ----------------------------------------------------------------
    def get_obs_value_counts(self, obs_label: str, do_sum: bool):
        """
        For a given obs label, e.g. "cell_type", count the number of occurrences of different values in
        SOMAs in the collection. If `do_sum` is false, count the number of SOMAs having that value. If
        `do_sum` is true, count the total number of instances of that value across the collection.
        """
        return self._get_obs_or_var_value_counts(obs_label, do_sum, True)

    def get_var_value_counts(self, var_label: str, do_sum: bool):
        """
        For a given var label, e.g. "feature_name", count the number of occurrences of different values in
        SOMAs in the collection. If `do_sum` is false, count the number of SOMAs having that value. If
        `do_sum` is true, count the total number of instances of that value across the collection.
        """
        return self._get_obs_or_var_value_counts(var_label, do_sum, False)

    def _get_obs_or_var_value_counts(
        self, obs_or_var_label: str, do_sum: bool, do_obs: bool
    ):
        """
        Supporting method for `get_obs_value_counts` and `get_var_value_counts`.
        """

        counts = {}
        for soma in self:
            attrs = [obs_or_var_label]
            obs_or_var = (
                soma.obs.df(attrs=attrs) if do_obs else soma.var.df(attrs=attrs)
            )

            if not obs_or_var_label in obs_or_var:
                continue

            if do_sum:
                obs_label_values = list(obs_or_var[obs_or_var_label])
            else:
                obs_label_values = sorted(list(set(obs_or_var[obs_or_var_label])))
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
