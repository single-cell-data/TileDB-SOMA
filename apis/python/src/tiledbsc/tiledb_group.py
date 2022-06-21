import tiledb
import tiledbsc.util_tiledb
from .soma_options import SOMAOptions
from .tiledb_object import TileDBObject
from .tiledb_array import TileDBArray

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
        return tiledb.object_type(self.uri, ctx=self._ctx) == "group"

    def _create(self):
        """
        Creates the TileDB group data structure on disk/S3/cloud.
        """
        if self._verbose:
            print(f"{self._indent}Creating TileDB group {self.uri}")
        tiledb.group_create(uri=self.uri, ctx=self._ctx)

        self._set_soma_object_type_metadata()

    def create_unless_exists(self):
        """
        Creates the TileDB group data structure on disk/S3/cloud, unless it already exists.
        """
        if not self.exists():
            self._create()

    def _set_soma_object_type_metadata(self):
        """
        This helps nested-structured traversals (especially those that start at the SOMACollection
        level) confidently navigate with a minimum of introspection on group contents.
        """
        with self._open("w") as G:
            G.meta[
                tiledbsc.util.SOMA_OBJECT_TYPE_METADATA_KEY
            ] = self.__class__.__name__
            G.meta[
                tiledbsc.util.SOMA_ENCODING_VERSION_METADATA_KEY
            ] = tiledbsc.util.SOMA_ENCODING_VERSION

    def _set_soma_object_type_metadata_recursively(self):
        """
        SOMAs/SOCOs written very early on in the development of this project may not have these set.
        Using this method we can after-populate these, without needig to re-ingest entire datasets.
        Any SOMAs/SOCOs ingested from June 2022 onward won't need this -- this metadata will be
        written at ingestion time.
        """
        self._set_soma_object_type_metadata()
        with self._open() as G:
            for O in G:  # This returns a tiledb.object.Object
                # It might appear simpler to have all this code within TileDBObject class,
                # rather than (with a little duplication) in TileDBGroup and TileDBArray.
                # However, getting it to work with a recursive data structure and finding the
                # required methods, it was simpler to split the logic this way.
                object_type = tiledb.object_type(O.uri, ctx=self._ctx)
                if object_type == "group":
                    group = TileDBGroup(uri=O.uri, name=O.name, parent=self)
                    group._set_soma_object_type_metadata_recursively()
                elif object_type == "array":
                    array = TileDBArray(uri=O.uri, name=O.name, parent=self)
                    array._set_soma_object_type_metadata()
                else:
                    raise Exception(
                        f"Unexpected object_type found: {object_type} at {O.uri}"
                    )

    def _open(self, mode="r"):
        """
        This is just a convenience wrapper around tiledb group-open.
        It works asa `with self._open() as G:` as well as `G = self._open(); ...; G.close()`.
        """
        assert mode in ("r", "w")
        if mode == "r" and not self.exists():
            raise Exception(f"Does not exist: {self.uri}")
        # This works in with-open-as contexts because tiledb.Group has __enter__ and __exit__ methods.
        return tiledb.Group(self.uri, mode=mode, ctx=self._ctx)

    def _add_object(self, obj: TileDBObject):
        """
        Adds a SOMA group/array to the current SOMA group -- e.g. base SOMA adding
        X, X adding a layer, obsm adding an element, etc.

        Semantics of `relative` from `self._soma_options.member_uris_are_relative`:

        * If `False` then the group will have the absolute path of the member. For populating matrix
        elements within a SOMA in TileDB cloud, this is necessary. For populating SOMA elements within
        a SOMACollection on local disk, this can be useful if you want to be able to move the SOMACollection
        storage around and have it remember the (unmoved) locations of SOMA objects elsewhere, i.e.
        if the SOMACollectio is in one place while its members are in other places. If the SOMAs
        in the collection are contained within the SOMACollection directory, you probably want `relative=True`.

        * If `True` then the group will have the relative path of the member. For TileDB Cloud, this
        is never the right thing to do. For local-disk storage, this is essential if you want to move
        a SOMA to another directory and have it remember the locations of the members within it.

        * If `None`, then we select `relative=False` if the URI starts with `tiledb://`, else we
        select `relative=True`. This is the default.
        """
        self.create_unless_exists()
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
            return {O.name: O.uri for O in G}

    def show_metadata(self, recursively=True, indent=""):
        """
        Shows metadata for the group, recursively by default.
        """
        print(f"{indent}[{self.name}]")
        for key, value in self.metadata().items():
            print(f"{indent}- {key}: {value}")
        if recursively:
            child_indent = indent + "  "
            with self._open() as G:
                for O in G:  # This returns a tiledb.object.Object
                    # It might appear simpler to have all this code within TileDBObject class,
                    # rather than (with a little duplication) in TileDBGroup and TileDBArray.
                    # However, getting it to work with a recursive data structure and finding the
                    # required methods, it was simpler to split the logic this way.
                    object_type = tiledb.object_type(O.uri, ctx=self._ctx)
                    if object_type == "group":
                        group = TileDBGroup(uri=O.uri, name=O.name, parent=self)
                        group.show_metadata(recursively, indent=child_indent)
                    elif object_type == "array":
                        array = TileDBArray(uri=O.uri, name=O.name, parent=self)
                        array.show_metadata(recursively, indent=child_indent)
                    else:
                        raise Exception(
                            f"Unexpected object_type found: {object_type} at {O.uri}"
                        )
