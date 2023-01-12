from abc import ABC, abstractmethod, abstractproperty
from typing import Optional, Union

import tiledb

from . import util
from .metadata_mapping import MetadataMapping
from .tiledb_session_context import TileDBSessionContext


class TileDBObject(ABC):
    """
    Base class for ``TileDBArray`` and ``Collection``.

    Accepts a TileDBSessionContext, to enable session state to be shared across SOMA objects.
    """

    _uri: str
    _soma_session_context: TileDBSessionContext
    _metadata: MetadataMapping

    def __init__(
        self,
        # All objects:
        uri: str,
        *,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional["TileDBObject"] = None,
        # Top-level objects should specify this:
        session_context: Optional[TileDBSessionContext] = None,
    ):
        """
        Initialization-handling shared between ``TileDBArray`` and ``Collection``.  Specify ``session_context`` for
        the top-level object; omit it and specify parent for non-top-level objects. Note that the parent reference
        is solely for propagating the session_context
        """

        self._uri = uri

        if parent is not None:
            assert session_context is None, "Only one of `session_context` and `parent` params can be passed as an arg"
            # inherit from parent
            self._soma_session_context = parent._soma_session_context
        else:
            self._soma_session_context = session_context or TileDBSessionContext()

        self._metadata = MetadataMapping(self)

    @property
    def _ctx(self) -> tiledb.Ctx:
        return self._soma_session_context.tiledb_ctx

    @property
    def metadata(self) -> MetadataMapping:
        """Metadata accessor"""
        return self._metadata
        # Note: this seems trivial, like we could just have `metadata` as an attribute.
        # However, we've found that since in `somacore` it's implemented as `@property`,
        # to avoid a static-analysis failure we have to do the same here.

    def delete(self) -> None:
        """
        Delete the storage specified with the URI.

        TODO: should this raise an error if the object does not exist?
        """
        try:
            tiledb.remove(self._uri)
        except tiledb.TileDBError:
            pass
        return

    def __repr__(self) -> str:
        """
        Default repr
        """
        if self.exists():
            return f'{self.soma_type}(uri="{self._uri}")'
        else:
            return f"{self.soma_type}(not created)"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, TileDBObject):
            return False
        return self._uri == other._uri

    @property
    def uri(self) -> str:
        """
        Accessor for the object's storage URI
        """
        return self._uri

    @abstractproperty
    def soma_type(self) -> str:
        """
        Returns the SOMA object type, e.g. "SOMADataFrame".
        """
        ...

    def exists(self) -> bool:
        """
        Returns true if the object exists and has the desired class name.

        This might be in case an object has not yet been populated, or, if a containing object has been populated but doesn't have a particular member (e.g. not all ``Measurement`` objects have a ``varp``).

        For ``tiledb://`` URIs this is a REST-server request which we'd like to cache.  However, remove-and-replace use-cases are possible and common in notebooks and it turns out caching the existence-check isn't a robust approach.
        """

        # Pre-checking if the group exists by calling tiledb.object_type is simple, however, for
        # tiledb-cloud URIs that occurs a penalty of two HTTP requests to the REST server, even
        # before a third, successful HTTP request for group-open.  Instead, we directly attempt the
        # group-open request, checking for an exception.
        try:
            return self._get_soma_type_from_metadata() == self.soma_type
        except tiledb.cc.TileDBError:
            return False

    @abstractmethod
    def _tiledb_open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        """Open the underlying TileDB array or Group"""
        ...

    def _common_create(self) -> None:
        """
        Utility method for various constructors.
        """
        self._set_object_type_metadata()

    def _set_object_type_metadata(self) -> None:
        """
        This helps nested-structure traversals (especially those that start at the Collection level) confidently navigate with a minimum of introspection on group contents.
        """
        # TODO: make a multi-set in MetadataMapping that would above a double-open there.
        with self._tiledb_open("w") as obj:
            obj.meta.update(
                {
                    util.SOMA_OBJECT_TYPE_METADATA_KEY: self.soma_type,
                    util.SOMA_ENCODING_VERSION_METADATA_KEY: util.SOMA_ENCODING_VERSION,
                }
            )

    def _get_soma_type_from_metadata(self) -> str:
        """
        Returns the class name associated with the group/array.
        """
        # mypy says:
        # error: Returning Any from function declared to return "str"  [no-any-return]
        return self._metadata.get(util.SOMA_OBJECT_TYPE_METADATA_KEY)  # type: ignore
