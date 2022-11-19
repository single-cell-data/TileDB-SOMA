from __future__ import annotations

from typing import Dict, Optional, Sequence

import tiledb

from .logging import logger
from .soma_options import SOMAOptions
from .tiledb_array import TileDBArray
from .tiledb_object import TileDBObject


class TileDBGroup(TileDBObject):
    """
    Wraps groups from TileDB-Py by retaining a URI, options, etc.
    """

    # Cache to avoid repeated calls to the REST server for resolving group-member URIs
    # in the tiledb-cloud case. We invalidate this on add-member or remove-member.
    _cached_member_names_to_uris: Optional[Dict[str, str]]

    def __init__(
        self,
        uri: str,
        name: str,
        *,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional[TileDBGroup] = None,
        # Top-level objects should specify these:
        soma_options: Optional[SOMAOptions] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        See the TileDBObject constructor.
        """
        super().__init__(uri, name, parent=parent, soma_options=soma_options, ctx=ctx)
        self._cached_member_names_to_uris = None

    def exists(self) -> bool:
        """
        Tells whether or not there is storage for the group. This might be in case a SOMA
        object has not yet been populated, e.g. before calling ``from_anndata`` -- or, if the
        SOMA has been populated but doesn't have this member (e.g. not all SOMAs have a ``varp``).
        """
        # For tiledb:// URIs this is a REST-server request which we'd like to cache.
        # However, remove-and-replace use-cases are possible and common in notebooks
        # and it turns out caching the existence-check isn't a robust approach.
        return bool(tiledb.object_type(self.uri, ctx=self._ctx) == "group")

    def create_unless_exists(self) -> None:
        """
        Creates the TileDB group data structure on disk/S3/cloud, unless it already exists.
        """
        if not self.exists():
            logger.debug(f"{self._indent}Creating TileDB group {self.nested_name}")
            try:
                tiledb.group_create(uri=self.uri, ctx=self._ctx)
                self._set_object_type_metadata()
            except tiledb.cc.TileDBError as e:
                stre = str(e)
                # This is fine in case of parallel creates
                if "already exists" not in stre:
                    # bare raise will raise the current exception without rewriting the stack trace
                    raise

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

    def _open(self, mode: str = "r") -> tiledb.Group:
        """
        This is just a convenience wrapper around tiledb group-open.
        It works asa ``with self._open() as G:`` as well as ``G = self._open(); ...; G.close()``.
        """
        assert mode in ("r", "w")
        # This works in with-open-as contexts because tiledb.Group has __enter__ and __exit__ methods.
        # Raises an exception, as desired, if the group does not exist (or, doesn't exist _yet_).
        return tiledb.Group(self.uri, mode=mode, ctx=self._ctx)

    def _get_child_uris(self, member_names: Sequence[str]) -> Dict[str, str]:
        """
        Batched version of ``get_child_uri``. Since there's REST-server latency for getting
        name-to-URI mapping for group-member URIs, in the tiledb://... case, total latency
        is reduced when we ask for all group-element name-to-URI mappings in a single
        request to the REST server.
        """

        try:
            # If the group exists, get URIs for elements. In local disk, S3, etc this is just
            # concatenating a slash and the member name to the self URI; for TileDB Cloud,
            # the member URIs involve UUIDs which are known to the server.
            answer = {}

            mapping = self.get_member_names_to_uris()
            for member_name in member_names:
                if member_name in mapping:
                    answer[member_name] = mapping[member_name]
                else:
                    # Truly a slash, not os.path.join:
                    # * If the client is Linux/Un*x/Mac, it's the same of course
                    # * On Windows, os.path.sep is a backslash but backslashes are _not_ accepted for S3 or
                    #   tiledb-cloud URIs, whereas in Windows versions for years now forward slashes _are_
                    #   accepted for local-disk paths.
                    # This means forward slash is acceptable in all cases.
                    answer[member_name] = self.uri + "/" + member_name

            return answer

        except tiledb.cc.TileDBError as e:
            stre = str(e)
            # Local-disk/S3 does-not-exist exceptions say 'Group does not exist'; TileDB Cloud
            # does-not-exist exceptions are worded less clearly.
            # 'Group' or 'group'
            if "roup does not exist" in stre or "HTTP code 401" in stre:
                # Group not constructed yet, but we must return pre-creation URIs. Here, appending
                # "/" and name is appropriate in all cases: even for tiledb://... URIs,
                # pre-construction URIs are of the form
                # tiledb://namespace/s3://something/something/soma/membername.
                return {
                    member_name: self.uri + "/" + member_name
                    for member_name in member_names
                }
            else:
                raise e

    def _get_child_uri(self, member_name: str) -> str:
        """
        Computes the URI for a child of the given object. For local disk, S3, and
        tiledb://.../s3://...  pre-creation URIs, this is simply the parent's URI, a slash, and the
        member name.  For post-creation TileDB-Cloud URIs, this is computed from the parent's
        information.  (This is because in TileDB Cloud, members have URIs like
        tiledb://namespace/df584345-28b7-45e5-abeb-043d409b1a97.)

        Please use _get_child_uris whenever possible, to reduce the number of REST-server requests
        in the tiledb//... URIs.
        """
        if not self.exists():
            # TODO: comment
            return self.uri + "/" + member_name
        mapping = self.get_member_names_to_uris()
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

    def _add_object(
        self,
        obj: TileDBObject,
        relative: Optional[bool] = None,
        check_is_direct_child: Optional[bool] = False,
    ) -> None:
        """
        Adds a SOMA group/array to the current SOMA group -- e.g. base SOMA adding
        X, X adding a layer, obsm adding an element, etc.

        Semantics of ``relative`` from ``self._soma_options.member_uris_are_relative``:

        * If ``False`` then the group will have the absolute path of the member. For populating matrix
        elements within a SOMA in TileDB cloud, this is necessary. For populating SOMA elements within
        a SOMACollection on local disk, this can be useful if you want to be able to move the SOMACollection
        storage around and have it remember the (unmoved) locations of SOMA objects elsewhere, i.e.
        if the SOMACollectio is in one place while its members are in other places. If the SOMAs
        in the collection are contained within the SOMACollection directory, you probably want ``relative=True``.

        * If ``True`` then the group will have the relative path of the member. For TileDB Cloud, this
        is never the right thing to do. For local-disk storage, this is essential if you want to move
        a SOMA to another directory and have it remember the locations of the members within it.

        * If ``None``, then we select ``relative=False`` if the URI starts with ``tiledb://``, else we
        select ``relative=True``. This is the default.

        If the relative argument is supplied and is not None, it is used; secondly
        ``self._soma_options.member_uris_are_relative`` is consulted; thirdly the URI prefix
        is consulted as described above.
        """
        self.create_unless_exists()
        child_uri = obj.uri
        if relative is None:
            relative = self._soma_options.member_uris_are_relative
        if relative is None:
            relative = not child_uri.startswith("tiledb://")
        if relative:
            child_uri = obj.name

        # Please see https://github.com/single-cell-data/TileDB-SingleCell/issues/258
        if check_is_direct_child and not child_uri.startswith("tiledb://"):
            parent_cleaned_path = self.uri.strip("/")
            child_cleaned_path = obj.uri.strip("/")

            # Windows paths may have ``C:\something\something`` or ``C:/something/something``
            # (including in our CI jobs) and this is user-dependent. Therefore, if we're going
            # to compare for equality, we need to compare both ways.
            expected_paths = [
                parent_cleaned_path + "/" + obj.name,
                parent_cleaned_path + "\\" + obj.name,
            ]

            if child_cleaned_path not in expected_paths:
                raise Exception(
                    f"Member URI {obj.uri} must be a direct child of parent URI {self.uri}"
                )

        self._cached_member_names_to_uris = None  # invalidate on add-member
        with self._open("r") as G:
            exists = obj.name in G
        with self._open("w") as G:
            if not exists:
                G.add(uri=child_uri, relative=relative, name=obj.name)

    def _remove_object(self, obj: TileDBObject) -> None:
        self._remove_object_by_name(obj.name)

    def _remove_object_by_name(self, member_name: str) -> None:
        self._cached_member_names_to_uris = None  # invalidate on remove-member
        if self.uri.startswith("tiledb://"):
            mapping = self.get_member_names_to_uris()
            if member_name not in mapping:
                raise Exception(f"name {member_name} not present in group {self.uri}")
            member_uri = mapping[member_name]
            with self._open("w") as G:
                G.remove(member_uri)
        else:
            with self._open("w") as G:
                G.remove(member_name)

    def get_member_names(self) -> Sequence[str]:
        """
        Returns the names of the group elements. For a SOMACollection, these will SOMA names;
        for a SOMA, these will be matrix/group names; etc.
        """
        return list(self.get_member_names_to_uris().keys())

    def get_member_uris(self) -> Sequence[str]:
        """
        Returns the URIs of the group elements. For a SOMACollection, these will SOMA URIs;
        for a SOMA, these will be matrix/group URIs; etc.
        """
        return list(self.get_member_names_to_uris().values())

    def get_member_names_to_uris(self) -> Dict[str, str]:
        """
        Like ``get_member_names()`` and ``get_member_uris``, but returns a dict mapping from
        member name to member URI.
        """
        if self._cached_member_names_to_uris is None:
            with self._open("r") as G:
                self._cached_member_names_to_uris = {obj.name: obj.uri for obj in G}
        return self._cached_member_names_to_uris

    def show_metadata(self, recursively: bool = True, indent: str = "") -> None:
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

    def _show_uris(self, recursively: bool = True, indent: str = "") -> None:
        """
        Shows URIs for the group, recursively by default.
        """
        print(f"{indent}{self.name} group {self.uri}")
        if recursively:
            child_indent = indent + "  "
            with self._open() as G:
                for obj in G:  # This returns a tiledb.object.Object
                    object_type = tiledb.object_type(obj.uri, ctx=self._ctx)
                    if object_type == "group":
                        group = TileDBGroup(uri=obj.uri, name=obj.name, parent=self)
                        group._show_uris(recursively, indent=child_indent)
                    elif object_type == "array":
                        print(f"{child_indent}{obj.name} array {obj.uri}")
                    else:
                        raise Exception(
                            f"Unexpected object_type found: {object_type} at {obj.uri}"
                        )
