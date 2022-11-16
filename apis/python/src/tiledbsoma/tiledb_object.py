from abc import ABC, abstractmethod, abstractproperty
from typing import Optional, Union

import tiledb

from . import util
from .metadata_mapping import MetadataMapping
from .tiledb_platform_config import TileDBPlatformConfig


class TileDBObject(ABC):
    """
    Base class for ``TileDBArray`` and ``Collection``.

    Manages tiledb_platform_config, context, etc. which are common to both.
    """

    _uri: str
    _tiledb_platform_config: TileDBPlatformConfig
    metadata: MetadataMapping

    def __init__(
        self,
        # All objects:
        uri: str,
        *,
        # Non-top-level objects can have a parent to propgate context, depth, etc.
        parent: Optional["TileDBObject"] = None,
        # Top-level objects should specify these:
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Initialization-handling shared between ``TileDBArray`` and ``Collection``.  Specify ``tiledb_platform_config`` and ``ctx`` for the top-level object; omit them and specify parent for non-top-level objects. Note that the parent reference is solely for propagating options, ctx, display depth, etc.
        """
        self._uri = uri
        if parent is None:
            self._ctx = ctx
            # TODO - this does not belong in a core class.
            self._indent = ""
        else:
            tiledb_platform_config = parent._tiledb_platform_config
            self._ctx = parent._ctx
            # TODO - this does not belong in a core class.
            self._indent = parent._indent + "  "

        self._tiledb_platform_config = tiledb_platform_config or TileDBPlatformConfig()
        # Null ctx is OK if that's what they wanted (e.g. not doing any TileDB-Cloud ops).

        self.metadata = MetadataMapping(self)

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
        return self._uri

    @abstractproperty
    def soma_type(self) -> str:
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
        return self.metadata.get(util.SOMA_OBJECT_TYPE_METADATA_KEY)  # type: ignore
