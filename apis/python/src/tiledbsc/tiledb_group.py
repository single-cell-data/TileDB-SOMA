import tiledb
import tiledbsc.util_tiledb
from .soma_options import SOMAOptions
from .tiledb_object import TileDBObject

from contextlib import contextmanager

from typing import Optional, Union, List, Dict
import os


class TileDBGroup(TileDBObject):
    """
    Wraps groups from TileDB-Py by retaining a URI, verbose flag, etc.
    """

    def __init__(
        self,
        uri: str,
        name: str,
        # Non-top-level objects can have a parent to propgate context, depth, etc.
        # What we really want to say is:
        # parent: Optional[TileDBGroup] = None,
        parent=None,
        # Top-level objects should specify these:
        soma_options: Optional[SOMAOptions] = None,
        verbose: Optional[bool] = True,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            soma_options=soma_options,
            verbose=verbose,
            ctx=ctx,
        )

    def _object_type(self):
        """
        This should be implemented by child classes and should return what tiledb.object_type(uri)
        returns for objects of a given type -- nominally 'group' or 'array'.
        """
        return "group"

    def exists(self) -> bool:
        """
        Tells whether or not there is storage for the group. This might be in case a SOMA
        object has not yet been populated, e.g. before calling `from_anndata` -- or, if the
        SOMA has been populated but doesn't have this member (e.g. not all SOMAs have a `varp`).
        """
        return tiledb.object_type(self.uri) == "group"

    def _create(self):
        """
        Creates the TileDB group data structure on disk/S3/cloud.
        """
        if self._verbose:
            print(f"{self._indent}Creating TileDB group {self.uri}")
        tiledb.group_create(uri=self.uri, ctx=self._ctx)
        with self._open("w") as G:
            G.meta[
                tiledbsc.util_tiledb.SOMA_OBJECT_TYPE_METADATA_KEY
            ] = self.__class__.__name__

    def _open_withlessly(self, mode="r"):
        """
        This is just a convenience wrapper around tiledb.open of the tiledb group
        associated with this SOMA element.
        """
        assert mode in ["w", "r"]
        if mode == "w" and not self.exists():
            self._create()
        if mode == "r" and not self.exists():
            raise Exception(f"Does not exist: {self.uri}")
        return tiledb.Group(self.uri, mode=mode, ctx=self._ctx)

    @contextmanager
    def _open(self, mode="r"):
        """
        This is just a convenience wrapper around tiledb.open of the tiledb group
        associated with this SOMA element, supporting Python with-as syntax.
        TODO: One TileDB.Py's Group objects have `__enter__` and `__exit__`
        method, fold this and _open_withlessly together.
        """
        assert mode in ("r", "w")

        # Do this check here, not just in _open_withlessly -- otherwise we get
        #   UnboundLocalError: local variable 'G' referenced before assignment
        # which is super-confusing for the user.
        if mode == "r" and not self.exists():
            raise Exception(f"Does not exist: {self.uri}")

        try:
            G = self._open_withlessly(mode)
            yield G
        finally:
            G.close()

    def _add_object(self, obj: TileDBObject):
        """
        Adds a SOMA group/array to the current SOMA group -- e.g. base SOMA adding
        X, X adding a layer, obsm adding an element, etc.

        Semantics of `relative` from `self._soma_options.member_uris_are_relative`:

        * If `False` then the group will have the absolute path of the member. For populating matrix
        elements within a SOMA in TileDB cloud, this is necessary. For populating SOMA elements within
        a SOMACollection on local disk, this can be useful if you want to be able to move the SOMACollection
        storage around and have it remember the (unmoved) locations of SOMA objects elsewhere.

        * If `True` then the group will have the relative path of the member. For TileDB Cloud, this
        is never the right thing to do. For local-disk storage, this is essential if you want to move
        a SOMA to another directory and have it remember the locations of the members within it.

        * If `None`, then we select `relative=False` if the URI starts with `tiledb://`, else we
        select `relative=True`. This is the default.
        """
        relative = self._soma_options.member_uris_are_relative
        child_uri = obj.uri
        if relative is None:
            relative = not child_uri.startswith("tiledb://")
        if relative:
            child_uri = obj.name
        with self._open("w") as G:
            G.add(uri=child_uri, relative=relative, name=obj.name)

    def _remove_object(self, obj: TileDBObject):
        with self._open("w") as G:
            G.remove(obj.name)

    def _get_member_names(self):
        """
        Returns the names of the group elements. For a SOMACollection, these will SOMA names;
        for a SOMA, these will be matrix/group names; etc.
        """
        return list(self._get_member_names_to_uris().keys())

    def _get_member_uris(self) -> List[str]:
        """
        Returns the URIs of the group elements. For a SOMACollection, these will SOMA URIs;
        for a SOMA, these will be matrix/group URIs; etc.
        """
        return list(self._get_member_names_to_uris().values())

    def _get_member_names_to_uris(self) -> Dict[str, str]:
        """
        Like `_get_member_names()` and `_get_member_uris`, but returns a dict mapping from
        member name to member URI.
        """
        with self._open("r") as G:

            # We aren't able to get a listing of group keys, so as a workaround we parse the group printer.
            names = [
                e.split()[1] for e in G.__repr__().split("\n") if e.startswith("|-- ")
            ]

            return {name: G[name].uri for name in names}
