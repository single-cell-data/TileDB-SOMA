from abc import ABC, abstractmethod
from typing import Optional, Union

import somacore
import tiledb

from . import util
from .metadata_mapping import MetadataMapping
from .options import SOMATileDBContext


class TileDBObject(ABC, somacore.SOMAObject):
    """
    Base class for ``TileDBArray`` and ``Collection``.

    Accepts a SOMATileDBContext, to enable session state to be shared across SOMA objects.
    """

    _uri: str
    _context: SOMATileDBContext
    _metadata: MetadataMapping

    def __init__(
        self,
        # All objects:
        uri: str,
        *,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional["TileDBObject"] = None,
        # Top-level objects should specify this:
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        Initialization-handling shared between ``TileDBArray`` and ``Collection``.  Specify ``context`` for
        the top-level object; omit it and specify parent for non-top-level objects. Note that the parent reference
        is solely for propagating the context
        """

        self._uri = uri

        if parent is not None:
            if context is not None:
                raise TypeError(
                    "Only one of `context` and `parent` params can be passed as an arg"
                )
            # inherit from parent
            self._context = parent._context
        else:
            self._context = context or SOMATileDBContext()

        self._metadata = MetadataMapping(self)

    @property
    def context(self) -> SOMATileDBContext:
        return self._context

    @property
    def _ctx(self) -> tiledb.Ctx:
        return self._context.tiledb_ctx

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
            return (
                self._metadata.get(util.SOMA_OBJECT_TYPE_METADATA_KEY) == self.soma_type
            )
        except tiledb.cc.TileDBError:
            return False

    @abstractmethod
    def _tiledb_open(self, mode: str = "r") -> Union[tiledb.Array, tiledb.Group]:
        """Open the underlying TileDB array or Group"""
        ...

    def _common_create(self, soma_type: str) -> None:
        """
        This helps nested-structure traversals (especially those that start at the
        Collection level) confidently navigate with a minimum of introspection on group
        contents.
        """
        # TODO: make a multi-set in MetadataMapping that would avoid a double-open there.
        with self._tiledb_open("w") as obj:
            obj.meta.update(
                {
                    util.SOMA_OBJECT_TYPE_METADATA_KEY: soma_type,
                    util.SOMA_ENCODING_VERSION_METADATA_KEY: util.SOMA_ENCODING_VERSION,
                }
            )
