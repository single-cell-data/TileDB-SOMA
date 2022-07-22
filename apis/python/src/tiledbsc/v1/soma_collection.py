from typing import Dict, Iterator, Optional

from .tiledb_group import TileDBGroup
from .tiledb_object import TileDBObject


class SOMACollection(TileDBGroup):
    """
    Contains a key-value mapping where the keys are string names and the values are any SOMA-defined
    foundational or composed type, including SOMACollection, SOMADataFrame, SOMADenseNdArray,
    SOMASparseNdArray or SOMAExperiment.
    """

    # TODO: comment re the performance impact of this cache.
    _members: Dict[str, TileDBObject]

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """

        super().__init__(uri=uri, name=name, parent=parent)
        self._members = {}

    # create is inherited from TileDBGroup

    # delete(uri)
    # Delete the SOMACollection specified with the URI.

    #    # TODO: static/class method?
    #    #    def delete(uri: str) -> None
    #    #        """
    #    #        Delete the SOMADataFrame specified with the URI.
    #    #        """

    # exists(uri) -> bool
    # Return true if object exists and is a SOMACollection.

    #    # TODO: static/class method?
    #    #    def exists(uri: str) -> bool
    #    #        """
    #    #        Return true if object exists and is a SOMADataFrame.
    #    #        """

    # get metadata
    # Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)

    #    #    def get_metadata():
    #    #        """
    #    #        Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)
    #    #        """

    #    # get_type() is inherited from TileDBObject

    def get(self, member_name: str) -> TileDBObject:
        """
        Get the member object associated with the key
        """
        # TODO: KEEP A MEMBER-CACHE ON ADD
        # TODO: UPDATE MEMBER-CACHE ON REMOVE

        if member_name not in self._members:
            # Do this here to avoid a cyclic package dependency:
            from .factory import _construct_member

            member_uri = self._get_child_uri(member_name)
            self._members[member_name] = _construct_member(member_uri, self)

        return self._members[member_name]

    def __len__(self) -> int:
        """
        Returns the number of members in the collection.  Implements Python's `len(collection)`.
        """
        return len(self._get_member_names_to_uris())

    def __contains__(self, member_name: str) -> bool:
        """
        Tests for the existence of key in collection.
        Implements the `in` operator.
        """
        return member_name in self._get_child_uris()

    def set(self, member: TileDBObject, relative: Optional[bool] = None) -> None:
        """
        Adds a member to the collection.
        """
        self._add_object(member, relative)
        self._members[member.get_name()] = member

    # TODO: note for the SOMA v1 spec: it says `del` not `delete` but `del` is a reserved word in Python.
    def delete(self, member_name: str) -> None:
        """
        Removes a member from the collection, when invoked as `collection.delete("namegoeshere")`.
        """
        self._remove_object_by_name(member_name)

    def __delattr__(self, member_name: str) -> None:
        """
        Removes a member from the collection, when invoked as `del collection.namegoeshere`.
        """
        self.delete(member_name)

    def __delitem__(self, member_name: str) -> None:
        """
        Removes a member from the collection, when invoked as `del collection["namegoeshere"]`.
        """
        self.delete(member_name)

    def __iter__(self) -> Iterator[TileDBObject]:
        """
        Iterates over the collection.  Implements Python `for member in collection: ...` syntax.
        """
        for name, uri in self._get_member_names_to_uris().items():
            if name not in self._members:
                # member-constructor cache -- this is an important optimization for
                # cloud storage, as URI-resolvers require HTTP requests
                #
                # TODO: need to read object metadata, including the classname, so we
                # can invoke the right constructor here -- not just TileDBObject --
                # so we can get a polymorphic return value.
                self._members[name] = TileDBObject(
                    uri=uri, name=name, parent=self, ctx=self._ctx
                )
            yield self._members[name]
