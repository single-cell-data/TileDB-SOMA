from contextlib import ExitStack
from typing import ClassVar, Generic, Optional, Type, TypeVar

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
        handle: ReadWriteHandle[_HandleType],
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
        self._handle = handle
        self._metadata = MetadataMapping(self._handle)
        self._close_stack = ExitStack()
        self._close_stack.enter_context(self._handle)

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

    # TODO: This is dangerous; two objects with the same URI may be different.
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, TileDBObject):
            return False
        return self.uri == other.uri

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
        my_cs: Optional[ExitStack] = getattr(self, "_close_stack", None)
        if my_cs:
            my_cs.close()

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

    def _ensure_open_read(self) -> None:
        if self.mode != "r":
            raise ValueError(f"{self} is open for writing, not reading")


AnyTileDBObject = TileDBObject[TDBHandle]
