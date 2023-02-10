from contextlib import ExitStack
from typing import Any, Generic, MutableMapping, Optional, Type, TypeVar

import somacore
import tiledb
from somacore import options
from typing_extensions import Self

from . import constants, tdb_handles
from .exception import SOMAError
from .options import SOMATileDBContext

_WrapperType_co = TypeVar(
    "_WrapperType_co", bound=tdb_handles.AnyWrapper, covariant=True
)
"""The type of handle on a backend object that we have.

Covariant because ``_handle`` is read-only.
"""


class TileDBObject(somacore.SOMAObject, Generic[_WrapperType_co]):
    """
    Base class for all TileDB SOMA objects.

    Accepts a SOMATileDBContext, to enable session state to be shared
    across SOMA objects.

    [lifecycle: experimental]
    """

    __slots__ = ("_close_stack", "_handle")

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode = "r",
        *,
        context: Optional[SOMATileDBContext] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> Self:
        """Opens this specific type of SOMA object.

        :param uri: The URI to open.
        :param mode: The mode to open the object in.
            ``r``: Open for reading only (cannot write).
            ``w``: Open for writing only (cannot read).
        """
        del platform_config  # unused
        context = context or SOMATileDBContext()
        handle = cls._wrapper_type.open(uri, mode, context)
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    def __init__(
        self,
        handle: _WrapperType_co,
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
        self._close_stack.enter_context(self._handle)

    _wrapper_type: Type[_WrapperType_co]
    """Class variable of the Wrapper class used to open this object type."""

    @property
    def context(self) -> SOMATileDBContext:
        return self._handle.context

    @property
    def _ctx(self) -> tiledb.Ctx:
        return self.context.tiledb_ctx

    @property
    def metadata(self) -> MutableMapping[str, Any]:
        return self._handle.metadata

    def __repr__(self) -> str:
        return f"<{self._my_repr()}>"

    def _my_repr(self) -> str:
        """``__repr__``, but without the ``<>``."""
        open_str = "CLOSED" if self.closed else "open"
        return f"{type(self).__name__} {self.uri!r} ({open_str} for {self.mode!r})"

    @property
    def uri(self) -> str:
        """
        Accessor for the object's storage URI
        """
        return self._handle.uri

    def close(self) -> None:
        """
        Release any resources held while the object is open.
        Closing an already-closed object is a no-op.
        """
        self._close_stack.close()

    @property
    def closed(self) -> bool:
        """True if the object has been closed. False if it is still open."""
        return self._handle.closed

    @property
    def mode(self) -> options.OpenMode:
        """The mode this object was opened in, either ``r`` or ``w``."""
        return self._handle.mode

    @classmethod
    def exists(cls, uri: str, context: Optional[SOMATileDBContext] = None) -> bool:
        """Finds whether an object of this type exists at the given URI."""
        context = context or SOMATileDBContext()
        try:
            with cls._wrapper_type.open(uri, "r", context) as hdl:
                md_type = hdl.metadata.get(constants.SOMA_OBJECT_TYPE_METADATA_KEY)
                if not isinstance(md_type, str):
                    return False
                return md_type.lower() == cls.soma_type.lower()
        except SOMAError:
            return False

    @classmethod
    def _set_create_metadata(cls, handle: tdb_handles.AnyWrapper) -> None:
        """Sets the necessary metadata on a newly-created TileDB object."""
        handle.writer.meta.update(
            {
                constants.SOMA_OBJECT_TYPE_METADATA_KEY: cls.soma_type,
                constants.SOMA_ENCODING_VERSION_METADATA_KEY: constants.SOMA_ENCODING_VERSION,
            }
        )
        # Semi-hack: flush the metadata immediately upon creation so that the
        # backing storage isn't half-created (i.e., there is a tiledb object
        # on disk, but its type is not stored). This is immutable, so it's fine.
        handle._flush_hack()

    def _check_open_read(self) -> None:
        if self.mode != "r":
            raise ValueError(f"{self} is open for writing, not reading")


AnyTileDBObject = TileDBObject[tdb_handles.AnyWrapper]
