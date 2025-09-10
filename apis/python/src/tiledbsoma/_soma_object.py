# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import datetime
from collections.abc import MutableMapping
from contextlib import ExitStack
from typing import Any, Generic, TypeVar

import somacore
from somacore import options
from typing_extensions import Self

from . import _constants, _tdb_handles
from ._exception import SOMAError
from ._types import OpenTimestamp
from ._util import check_type, ms_to_datetime
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context

_WrapperType_co = TypeVar("_WrapperType_co", bound=_tdb_handles.AnyWrapper, covariant=True)

"""The type of handle on a backend object that we have.

Covariant because ``_handle_wrapper`` is read-only.
"""


class SOMAObject(somacore.SOMAObject, Generic[_WrapperType_co]):
    """Base class for all TileDB SOMA objects.

    Accepts a SOMATileDBContext, to enable session state to be shared
    across SOMA objects.

    Lifecycle:
        Maturing.
    """

    _wrapper_type: (
        type[_WrapperType_co]
        | type[_tdb_handles.DataFrameWrapper]
        | type[_tdb_handles.DenseNDArrayWrapper]
        | type[_tdb_handles.SparseNDArrayWrapper]
        | type[_tdb_handles.CollectionWrapper]
        | type[_tdb_handles.ExperimentWrapper]
        | type[_tdb_handles.MeasurementWrapper]
        | type[_tdb_handles.SceneWrapper]
        | type[_tdb_handles.MultiscaleImageWrapper]
    )
    """Class variable of the Wrapper class used to open this object type."""

    __slots__ = ("_close_stack", "_handle", "_handle_wrapper")

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode = "r",
        *,
        tiledb_timestamp: OpenTimestamp | None = None,
        context: SOMATileDBContext | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> Self:
        """Opens this specific type of SOMA object.

        Args:
            uri:
                The URI to open.
            mode:
                The mode to open the object in.
                - ``r``: Open for reading only (cannot write or delete).
                - ``w``: Open for writing only (cannot read or delete).
                - ``d``: Open for deleting only (cannot read or write).
            tiledb_timestamp:
                The TileDB timestamp to open this object at,
                either an int representing milliseconds since the Unix epoch
                or a datetime.datetime object.
                When not provided (the default), the current time_wrapper is used.
                A value of zero results in default, i.e., the current time.

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
            Maturing.
        """
        del platform_config  # unused
        context = _validate_soma_tiledb_context(context)
        handle = cls._wrapper_type.open(uri, mode, context, tiledb_timestamp)
        return cls(
            handle,  # type: ignore[arg-type]
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    def __init__(
        self,
        handle: (
            _WrapperType_co
            | _tdb_handles.DataFrameWrapper
            | _tdb_handles.DenseNDArrayWrapper
            | _tdb_handles.SparseNDArrayWrapper
        ),
        *,
        _dont_call_this_use_create_or_open_instead: str = "unset",
    ) -> None:
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
                f" internal use only.",
            )
        self._handle_wrapper = handle
        self._handle = self._handle_wrapper._handle
        self._close_stack.enter_context(self._handle_wrapper)
        self._parse_special_metadata()

    def _parse_special_metadata(self) -> None:
        """Helper function the subclasses can override if they require additional validation or set-up."""
        return

    def reopen(self, mode: options.OpenMode, tiledb_timestamp: OpenTimestamp | None = None) -> Self:
        """Return a new copy of the SOMAObject with the given mode at the current
        Unix timestamp.

        Args:
            mode:
                The mode to open the object in.
                - ``r``: Open for reading only (cannot write or delete).
                - ``w``: Open for writing only (cannot read or delete).
                - ``d``: Open for deleting only (cannot read or write).
            tiledb_timestamp:
                The TileDB timestamp to open this object at, either an int representing milliseconds since the Unix
                epoch or a datetime.datetime object. When not provided (the default), the current time is used.

        Raises:
            ValueError:
                If the user-provided ``mode`` is invalid.
            SOMAError:
                If the object has unwritten metadata.

        Lifecycle:
            Experimental.
        """
        self._handle_wrapper.close()
        self._handle_wrapper = self._wrapper_type.open(
            self._handle_wrapper.uri, mode, self._handle_wrapper.context, tiledb_timestamp
        )  # type: ignore
        self._handle = self._handle_wrapper._handle
        self._close_stack.enter_context(self._handle_wrapper)
        self._parse_special_metadata()
        return self

    @property
    def context(self) -> SOMATileDBContext:
        return self._handle_wrapper.context

    @property
    def metadata(self) -> MutableMapping[str, Any]:
        return self._handle_wrapper.metadata

    def __repr__(self) -> str:
        return f"<{self._my_repr()}>"

    def _my_repr(self) -> str:
        """``__repr__``, but without the ``<>``."""
        open_str = "CLOSED" if self.closed else "open"
        return f"{type(self).__name__} {self.uri!r} ({open_str} for {self.mode!r})"

    @property
    def uri(self) -> str:
        """Accessor for the object's storage URI.

        Examples:
            >>> soma_object.uri
            file://tmp/an_object_uri

        Lifecycle:
            Maturing.
        """
        return self._handle_wrapper.uri

    def close(self) -> None:
        """Release any resources held while the object is open.
        Closing an already-closed object is a no-op.

        Examples:
            >>> soma_object.close()

        Lifecycle:
            Maturing.
        """
        self._close_stack.close()

    @property
    def closed(self) -> bool:
        """True if the object has been closed. False if it is still open.

        Examples:
            >>> with tiledbsoma.open("an_object") as soma_object:
            ...     print(soma_object.closed)
            ...
            False
            >>> print(soma_object.closed)
            True

        Lifecycle:
            Maturing.
        """
        return self._handle_wrapper.closed

    @property
    def mode(self) -> options.OpenMode:
        """The mode this object was opened in, either ``r``, ``w``, or ``d``.

        Examples:
            >>> with tiledbsoma.open("an_object") as soma_object:
            ...     print(soma_object.mode)
            ...
            r

        Lifecycle:
            Maturing.
        """
        return self._handle_wrapper.mode

    def _verify_open_for_deleting(self) -> None:
        """Raises an error if the object is not open for deleting."""
        if self.closed:
            raise SOMAError(f"{self.__class__.__name__} ({self.uri}) must be open for deleting (closed).")
        if self.mode != "d":
            raise SOMAError(f"{self.__class__.__name__} ({self.uri}) must be open for deleting. Mode is '{self.mode}'.")

    def verify_open_for_writing(self) -> None:
        """Raises an error if the object is not open for writing."""
        if self.closed:
            raise SOMAError(f"{self.__class__.__name__} ({self.uri}) must be open for writing (closed)")
        if self.mode != "w":
            raise SOMAError(f"{self.__class__.__name__} ({self.uri}) must be open for writing. Mode is '{self.mode}'.")

    def _verify_open_for_reading(self) -> None:
        """Raises an error if the object is not open for reading."""
        if self.closed:
            raise SOMAError(f"{self.__class__.__name__} ({self.uri}) must be open for reading (closed)")
        if self.mode != "r":
            raise SOMAError(f"{self.__class__.__name__} ({self.uri}) must be open for reading. Mode is '{self.mode}'.")

    @property
    def tiledb_timestamp(self) -> datetime.datetime:
        """The time that this object was opened in UTC."""
        return ms_to_datetime(self.tiledb_timestamp_ms)

    @property
    def tiledb_timestamp_ms(self) -> int:
        """The time this object was opened, as millis since the Unix epoch."""
        return self._handle_wrapper.timestamp_ms

    @classmethod
    def exists(
        cls,
        uri: str,
        context: SOMATileDBContext | None = None,
        tiledb_timestamp: OpenTimestamp | None = None,
    ) -> bool:
        """Finds whether an object of this type exists at the given URI.

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
                A value of zero results in default, i.e., current time.

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
            Maturing.
        """
        check_type("uri", uri, (str,))
        context = _validate_soma_tiledb_context(context)
        try:
            with cls._wrapper_type.open(uri, "r", context, tiledb_timestamp) as hdl:
                md_type = hdl.metadata.get(_constants.SOMA_OBJECT_TYPE_METADATA_KEY)
                if not isinstance(md_type, str):
                    return False
                return md_type.lower() == cls.soma_type.lower()
        except (RuntimeError, SOMAError):
            return False


AnySOMAObject = SOMAObject[_tdb_handles.AnyWrapper]
