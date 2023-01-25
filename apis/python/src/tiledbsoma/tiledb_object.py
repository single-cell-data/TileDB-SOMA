from contextlib import ExitStack
from typing import Generic, Optional, TypeVar

import somacore
import tiledb
from somacore import options

from tiledbsoma import constants

from .metadata_mapping import MetadataMapping
from .options import SOMATileDBContext
from .types import StorageType, TDBHandle
from .util_tiledb import ReadWriteHandle

_HandleType = TypeVar("_HandleType", bound=TDBHandle)


class TileDBObject(somacore.SOMAObject, Generic[_HandleType]):
    """
    Base class for ``TileDBArray`` and ``Collection``.

    Accepts a SOMATileDBContext, to enable session state to be shared across SOMA objects.

    [lifecycle: experimental]
    """

    def __init__(
        self,
        uri: str,
        mode: options.OpenMode,
        handle: ReadWriteHandle[_HandleType],
        context: SOMATileDBContext,
        *,
        _this_is_internal_only: str = "unset",
    ):
        """Common initialization.

        This function is internal; users should open TileDB SOMA object using
        the :meth:`create` and :meth:`open` factory class methods.
        """
        if _this_is_internal_only != "tiledbsoma-internal-code":
            name = type(self).__name__
            raise TypeError(
                f"{name} initializers are intended for internal use only."
                f" To open an existing {name}, use tiledbsoma.open(...)"
                f" or the {name}.open(...) class method."
                f" To create a new {name}, use the {name}.create class method."
            )
        self._uri = uri
        self._mode: options.OpenMode = mode
        self._handle = handle
        self._context = context
        self._metadata = MetadataMapping(self._handle)
        self._close_stack = ExitStack()
        self._close_stack.enter_context(self._handle)
        self._was_deleted = False

    _STORAGE_TYPE: StorageType

    @property
    def context(self) -> SOMATileDBContext:
        return self._context

    @property
    def _ctx(self) -> tiledb.Ctx:
        return self._context.tiledb_ctx

    @property
    def metadata(self) -> MetadataMapping:
        # This needs to be implemented as a @property because Python's ABCs
        # require that abstract properties be implemented on the object itself,
        # rather than being a field (i.e., creating self.whatever in __init__).
        return self._metadata

    # TODO: This needs reconsidering, since it means we're deleting something
    # out from underneath ourselves.
    def delete(self) -> None:
        """
        Delete the storage specified with the URI.

        [lifecycle: experimental]
        """

        # TODO: should this raise an error if the object does not exist?
        self.close()
        try:
            self._was_deleted = True
            tiledb.remove(self._uri)
        except tiledb.TileDBError:
            pass
        return

    def __repr__(self) -> str:
        return f'{self.soma_type}(uri="{self._uri}")'

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

        [lifecycle: experimental]
        """
        # We already opened this, so it should always exist.
        return not self._was_deleted

    def flush(self) -> None:
        """Flushes any pending writes to the TileDB store.

        [lifecycle: experimental]
        """
        self._handle.flush(update_read=True)

    def close(self) -> None:
        """
        Release any resources held while the object is open. Closing an already-closed object is a
        no-op.
        """
        my_cs: Optional[ExitStack] = getattr(self, "_close_stack", None)
        if my_cs:
            my_cs.close()

    @property
    def mode(self) -> options.OpenMode:
        """
        Current open mode: read (r), write (w), or closed (None).
        """
        return self._mode

    @classmethod
    def _set_create_metadata(cls, handle: TDBHandle) -> None:
        """Sets the necessary metadata on a newly-created TileDB object."""
        handle.meta.update(
            {
                constants.SOMA_OBJECT_TYPE_METADATA_KEY: cls.soma_type,
                constants.SOMA_ENCODING_VERSION_METADATA_KEY: constants.SOMA_ENCODING_VERSION,
            }
        )


AnyTileDBObject = TileDBObject[TDBHandle]
