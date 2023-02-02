from contextlib import ExitStack
from typing import ClassVar, Generic, Type, TypeVar

import somacore
import tiledb
from somacore import options

from . import constants
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
        handle: ReadWriteHandle[_HandleType],
        *,
        _dont_call_this_use_create_or_open_instead: str = "unset",
    ):
        """Internal-only common initializer steps.

        This function is internal; users should open TileDB SOMA objects using
        the :meth:`create` and :meth:`open` factory class methods.
        """
        self._close_stack = ExitStack()
        """An exit stack to manage closing handles owned by this object.

        This is used to manage both our direct handle (in the case of simple
        TileDB objects) and the lifecycle of owned children (in the case of
        Collections).
        """
        # If a user calls something like tiledbsoma.Experiment(something),
        # give them a hint about the right thing to do.
        if _dont_call_this_use_create_or_open_instead != "tiledbsoma-internal-code":
            name = type(self).__name__
            raise TypeError(
                f"{name} objects must be created using a factory function."
                f" To open an existing {name}, use tiledbsoma.open(uri, ...)"
                f" or the {name}.open(uri, ...) class method."
                f" To create a new {name}, use the {name}.create class method."
                f" Directly calling `{name}(...)` is intended for TileDB SOMA"
                f" internal use only."
            )
        self._handle = handle
        self._metadata = MetadataMapping(self._handle)
        self._close_stack.enter_context(self._handle)
        self._closed = False

    _STORAGE_TYPE: StorageType
    _tiledb_type: ClassVar[Type[TDBHandle]]

    @property
    def context(self) -> SOMATileDBContext:
        return self._handle.context

    @property
    def _ctx(self) -> tiledb.Ctx:
        return self.context.tiledb_ctx

    @property
    def metadata(self) -> MetadataMapping:
        # This needs to be implemented as a @property because Python's ABCs
        # require that abstract properties be implemented on the object itself,
        # rather than being a field (i.e., creating self.whatever in __init__).
        return self._metadata

    def __repr__(self) -> str:
        return f'{self.soma_type}(uri="{self.uri}")'

    @property
    def uri(self) -> str:
        """
        Accessor for the object's storage URI
        """
        return self._handle.uri

    def close(self) -> None:
        """
        Release any resources held while the object is open. Closing an already-closed object is a
        no-op.
        """
        self._close_stack.close()
        self._closed = True

    @property
    def closed(self) -> bool:
        """True if the object has been closed. False if it is still open."""
        return self._closed

    @property
    def mode(self) -> options.OpenMode:
        """
        Current open mode: read (r), write (w), or closed (None).
        """
        return self._handle.mode

    @classmethod
    def _set_create_metadata(cls, handle: ReadWriteHandle[TDBHandle]) -> None:
        """Sets the necessary metadata on a newly-created TileDB object."""
        handle.writer.meta.update(
            {
                constants.SOMA_OBJECT_TYPE_METADATA_KEY: cls.soma_type,
                constants.SOMA_ENCODING_VERSION_METADATA_KEY: constants.SOMA_ENCODING_VERSION,
            }
        )
        # HACK: We need this so that the metadata appears on the read handle.
        handle._flush_hack()

    def _check_open_read(self) -> None:
        if self.mode != "r":
            raise ValueError(f"{self} is open for writing, not reading")


AnyTileDBObject = TileDBObject[TDBHandle]
