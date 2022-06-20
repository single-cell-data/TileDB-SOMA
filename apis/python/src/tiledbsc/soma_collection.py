from typing import Optional, List

import tiledb

from .soma_options import SOMAOptions
from .soma import SOMA
from .tiledb_group import TileDBGroup


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
    def __iter__(self) -> List[SOMA]:
        """
        Implements `for soma in soco: ...`
        """
        retval = []
        for name, uri in self._get_member_names_to_uris().items():
            soma = SOMA(uri=uri, name=name, parent=self)
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
