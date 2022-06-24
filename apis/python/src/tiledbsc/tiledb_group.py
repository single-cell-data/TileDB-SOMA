from typing import Dict, List, Optional

import tiledb

import tiledbsc.util_tiledb

from .logging import logger
from .soma_options import SOMAOptions
from .tiledb_array import TileDBArray
from .tiledb_object import TileDBObject


class TileDBGroup(TileDBObject):
    """
    Wraps groups from TileDB-Py by retaining a URI, options, etc.
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
            ctx=ctx,
        )

    def _object_type(self) -> str:
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

    def _create(self) -> None:
        """
        Creates the TileDB group data structure on disk/S3/cloud.
        """
        logger.info(f"{self._indent}Creating TileDB group {self.uri}")
        tiledb.group_create(uri=self.uri, ctx=self._ctx)

        self._set_object_type_metadata()

    def create_unless_exists(self) -> None:
        """
        Creates the TileDB group data structure on disk/S3/cloud, unless it already exists.
        """
        if not self.exists():
            self._create()

    def _set_object_type_metadata(self) -> None:
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

    def get_object_type(self) -> str:
        """
        Returns the class name associated with the group.
        """
        with self._open("r") as G:
            return G.meta[tiledbsc.util_tiledb.SOMA_OBJECT_TYPE_METADATA_KEY]

    def _set_object_type_metadata_recursively(self) -> None:
        """
        SOMAs/SOCOs written very early on in the development of this project may not have these set.
        Using this method we can after-populate these, without needig to re-ingest entire datasets.
        Any SOMAs/SOCOs ingested from June 2022 onward won't need this -- this metadata will be
        written at ingestion time.
        """
        self._set_object_type_metadata()
        with self._open() as G:
            for obj in G:  # This returns a tiledb.object.Object
                # It might appear simpler to have all this code within TileDBObject class,
                # rather than (with a little duplication) in TileDBGroup and TileDBArray.
                # However, getting it to work with a recursive data structure and finding the
                # required methods, it was simpler to split the logic this way.
                object_type = tiledb.object_type(obj.uri, ctx=self._ctx)
                if object_type == "group":
                    group = TileDBGroup(uri=obj.uri, name=obj.name, parent=self)
                    group._set_object_type_metadata_recursively()
                elif object_type == "array":
                    array = TileDBArray(uri=obj.uri, name=obj.name, parent=self)
                    array._set_object_type_metadata()
                else:
                    raise Exception(
                        f"Unexpected object_type found: {object_type} at {obj.uri}"
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

    def _get_child_uri(self, member_name: str) -> str:
        """
        Computes the URI for a child of the given object. For local disk, S3, and
        tiledb://.../s3://...  pre-creation URIs, this is simply the parent's URI, a slash, and the
        member name.  For post-creation TileDB-Cloud URIs, this is computed from the parent's
        information.  (This is because in TileDB Cloud, members have URIs like
        tiledb://namespace/df584345-28b7-45e5-abeb-043d409b1a97.)
        """
        if not self.exists():
            # TODO: comment
            return self.uri + "/" + member_name
        mapping = self._get_member_names_to_uris()
        if member_name in mapping:
            return mapping[member_name]
        else:
            # Truly a slash, not os.path.join:
            # * If the client is Linux/Un*x/Mac, it's the same of course
            # * On Windows, os.path.sep is a backslash but backslashes are _not_ accepted for S3 or
            #   tiledb-cloud URIs, whereas in Windows versions for years now forward slashes _are_
            #   accepted for local-disk paths.
            # This means forward slash is acceptable in all cases.
            return self.uri + "/" + member_name

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
        # See _get_child_uri. Key point is that, on TileDB Cloud, URIs change from pre-creation to
        # post-creation. Example:
        # * Upload to pre-creation URI tiledb://namespace/s3://bucket/something/something/somaname
        # * Results will be at post-creation URI tiledb://namespace/somaname
        # * Note people can still use the pre-creation URI to read the data if they like.
        # * Member pre-creation URI tiledb://namespace/s3://bucket/something/something/somaname/obs
        # * Member post-creation URI tiledb://somaname/e4de581a-1353-4150-b1f4-6ed12548e497
        obj.uri = self._get_child_uri(obj.name)

    def _remove_object(self, obj: TileDBObject) -> None:
        with self._open("w") as G:
            G.remove(obj.name)

    def _get_member_names(self) -> List[str]:
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
            return {obj.name: obj.uri for obj in G}

    def show_metadata(self, recursively=True, indent="") -> None:
        """
        Shows metadata for the group, recursively by default.
        """
        logger.info(f"{indent}[{self.name}]")
        for key, value in self.metadata().items():
            logger.info(f"{indent}- {key}: {value}")
        if recursively:
            child_indent = indent + "  "
            with self._open() as G:
                for obj in G:  # This returns a tiledb.object.Object
                    # It might appear simpler to have all this code within TileDBObject class,
                    # rather than (with a little duplication) in TileDBGroup and TileDBArray.
                    # However, getting it to work with a recursive data structure and finding the
                    # required methods, it was simpler to split the logic this way.
                    object_type = tiledb.object_type(obj.uri, ctx=self._ctx)
                    if object_type == "group":
                        group = TileDBGroup(uri=obj.uri, name=obj.name, parent=self)
                        group.show_metadata(recursively, indent=child_indent)
                    elif object_type == "array":
                        array = TileDBArray(uri=obj.uri, name=obj.name, parent=self)
                        array.show_metadata(recursively, indent=child_indent)
                    else:
                        raise Exception(
                            f"Unexpected object_type found: {object_type} at {obj.uri}"
                        )
