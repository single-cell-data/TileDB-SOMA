# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

import datetime
from contextlib import ExitStack
from typing import Any, Generic, MutableMapping, Optional, Type, TypeVar

import somacore
import tiledb
from somacore import options
from typing_extensions import Self

from . import _constants, _tdb_handles
from ._exception import SOMAError
from ._types import OpenTimestamp
from ._util import check_type, ms_to_datetime
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context

_WrapperType_co = TypeVar(
    "_WrapperType_co", bound=_tdb_handles.AnyWrapper, covariant=True
)
"""The type of handle on a backend object that we have.

Covariant because ``_handle`` is read-only.
"""


class TileDBObject(somacore.SOMAObject, Generic[_WrapperType_co]):
    """Base class for all TileDB SOMA objects.

    Accepts a SOMATileDBContext, to enable session state to be shared
    across SOMA objects.

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_close_stack", "_handle")

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode = "r",
        *,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
        context: Optional[SOMATileDBContext] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> Self:
        """Opens this specific type of SOMA object.

        Args:
            uri:
                The URI to open.
            mode:
                The mode to open the object in.
                - ``r``: Open for reading only (cannot write).
                - ``w``: Open for writing only (cannot read).
            tiledb_timestamp:
                The TileDB timestamp to open this object at,
                measured in milliseconds since the Unix epoch.
                When unset (the default), the current time is used.

        Returns:
            The opened SOMA object.

        Raises:
            DoesNotExistError:
                If the object named by URI can not be accessed.
            SOMAError:
                If the underlying TileDB object specified by ``uri`` is
                not recognized as a SOMA object.
            ValueError:
                If the user-provided ``mode`` is invalid.

        Lifecycle:
            Experimental.
        """
        del platform_config  # unused
        context = _validate_soma_tiledb_context(context)
        handle = cls._wrapper_type.open(uri, mode, context, tiledb_timestamp)
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
        Accessor for the object's storage URI.

        Examples:
            >>> soma_object.uri
            file://tmp/an_object_uri

        Lifecycle:
            Experimental.
        """
        return self._handle.uri

    def close(self) -> None:
        """
        Release any resources held while the object is open.
        Closing an already-closed object is a no-op.

        Examples:
            >>> soma_object.close()

        Lifecycle:
            Experimental.
        """
        self._close_stack.close()

    @property
    def closed(self) -> bool:
        """
        True if the object has been closed. False if it is still open.

        Examples:
            >>> with tiledbsoma.open("an_object") as soma_object:
            ...     print(soma_object.closed)
            ...
            False
            >>> print(soma_object.closed)
            True

        Lifecycle:
            Experimental.
        """
        return self._handle.closed

    @property
    def mode(self) -> options.OpenMode:
        """
        The mode this object was opened in, either ``r`` or ``w``.

        Examples:
            >>> with tiledbsoma.open("an_object") as soma_object:
            ...     print(soma_object.mode)
            ...
            r

        Lifecycle:
            Experimental.
        """
        return self._handle.mode

    @property
    def tiledb_timestamp(self) -> datetime.datetime:
        """The time that this object was opened in UTC."""
        return ms_to_datetime(self.tiledb_timestamp_ms)

    @property
    def tiledb_timestamp_ms(self) -> int:
        """The time this object was opened, as millis since the Unix epoch."""
        return self._handle.timestamp_ms

    @classmethod
    def exists(
        cls,
        uri: str,
        context: Optional[SOMATileDBContext] = None,
        tiledb_timestamp: Optional[OpenTimestamp] = None,
    ) -> bool:
        """
        Finds whether an object of this type exists at the given URI.

        Args:
            uri:
                The URI to open.
            context:
                If provided, the :class:`SOMATileDBContext` to use when creating and
                attempting to access this object.
            tiledb_timestamp:
                The TileDB timestamp to open this object at,
                measured in milliseconds since the Unix epoch.
                When unset (the default), the current time is used.

        Raises:
            TypeError:
                If the ``uri`` is not a string.

        Examples:
            >>> with tiledbsoma.open("a_dataframe") as soma_df:
            ...     print(soma_df.soma_type)
            ...
            SOMADataFrame
            >>> tiledbsoma.DataFrame.exists("./a_dataframe")
            True
            >>> tiledbsoma.SparseNDArray.exists("./a_dataframe")
            False

        Lifecycle:
            Experimental.
        """
        check_type("uri", uri, (str,))
        context = _validate_soma_tiledb_context(context)
        try:
            with cls._wrapper_type.open(uri, "r", context, tiledb_timestamp) as hdl:
                md_type = hdl.metadata.get(_constants.SOMA_OBJECT_TYPE_METADATA_KEY)
                if not isinstance(md_type, str):
                    return False
                return md_type.lower() == cls.soma_type.lower()
        except (SOMAError, tiledb.cc.TileDBError):
            return False

    @classmethod
    def _set_create_metadata(cls, handle: _tdb_handles.AnyWrapper) -> None:
        """Sets the necessary metadata on a newly-created TileDB object."""
        handle.metadata.update(
            {
                _constants.SOMA_OBJECT_TYPE_METADATA_KEY: cls.soma_type,
                _constants.SOMA_ENCODING_VERSION_METADATA_KEY: _constants.SOMA_ENCODING_VERSION,
            }
        )
        # Semi-hack: flush the metadata immediately upon creation so that the
        # backing storage isn't half-created (i.e., there is a tiledb object
        # on disk, but its type is not stored).
        # TODO: We should probably write this metadata at time 0.
        # Doing so would eliminate this last _flush_hack call.
        handle._flush_hack()

    def _check_open_read(self) -> None:
        if self.mode != "r":
            raise ValueError(f"{self} is open for writing, not reading")


AnyTileDBObject = TileDBObject[_tdb_handles.AnyWrapper]
