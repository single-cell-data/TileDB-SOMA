from typing import Dict, Optional

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

    # TODO: FOR create():
    #    self._common_create() # object-type metadata etc

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

    # TODO: polymorphic return-type annotation
    def get(self, member_name: str):
        """
        Get the member object associated with the key
        """
        # TODO: MEMBER-CACHE ON ADD
        # TODO: UPDATE MEMBER-CACHE ON REMOVE
        from .factory import _construct_member

        member_uri = self._get_child_uri(member_name)
        return _construct_member(member_uri, self)

    def __len__(self) -> int:
        """
        Returns the number of members in the collection.  Implements Python's `len(collection)`.
        """
        return len(self._get_member_names_to_uris())

    # has(string key)
    # Test for the existence of key in collection.

    # ----------------------------------------------------------------
    def set(self, member: TileDBObject, relative: Optional[bool] = None) -> None:
        """
        Adds a member to the collection.
        """
        self._add_object(member, relative)

    # TODO: note for the SOMA v1 spec: it says `del` not `delete` but `del` is a reserved word in Python.
    def delete(self, member_name: str) -> None:
        """
        Removes a member from the collection, when invoked as `collection.remove("namegoeshere")`.
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

    #    # ----------------------------------------------------------------
    #    def __iter__(self) -> Iterator[TileDBObject]:
    #        """
    #        Iterates over the collection.  Implements Python `for member in collection: ...` syntax.
    #        """
    #        for name, uri in self._get_member_names_to_uris().items():
    #            if name not in self._members:
    #                # member-constructor cache -- this is an important optimization for
    #                # cloud storage, as URI-resolvers require HTTP requests
    #                #
    #                # TODO: need to read object metadata, including the classname, so we
    #                # can invoke the right constructor here -- not just TileDBObject --
    #                # so we can get a polymorphic return value.
    #                self._members[name] = TileDBObject(uri=uri, name=name, parent=self, ctx=self._ctx)
    #            yield self._members[name]


#    # TODO: UNION RETURN TYPE
#    # TODO: FIGURE OUT HOW TO DO THIS WITHOUT CIRCULAR IMPORTS
#    def _construct_member(self, member_uri: str):
#        """
#        TODO: COMMENT
#        """
#        # TODO: xref to TileDBObject _set_object_type_metadata/_get_object_type_metadata.
#        # and/or, put some of this there as a class/static method.
#        pass
#        # sketch:
#        # get class name from meta -- with due respect for:
#        # * is-array vs is-group
#
#        # TODO:
#        # * reduce cloud ops necessary to answer that question
#        # need get-object-type ...
#
#        # object_type = tiledb.object_type
#        # if object_type == 'array':
#        # elif object_type == 'array':
#        # else:
#        #     raise Exception(f"object type {object_type} unrecognized")
#
#
#        # if class_name == "SOMAExperiment":
#        #     return SOMAExperiment(uri=member_uri, parent=self)
#        # elif class_name == "SOMAMeasurement":
#        #     return SOMAMeasurement(uri=member_uri, parent=self)
#        # elif class_name == "SOMACollection":
#        #     return SOMACollection(uri=member_uri, parent=self)
#        # elif class_name == "SOMADataFrame":
#        #     return SOMADataFrame(uri=member_uri, parent=self)
#        # elif class_name == "SOMADenseNdArray":
#        #     return SOMADenseNdArray(uri=member_uri, parent=self)
#        # elif class_name == "SOMASparseNdArray":
#        #     return SOMASparseNdArray(uri=member_uri, parent=self)
#        # else:
#        #    raise Exception(f"class name \"{class_name}\" unrecognized")
