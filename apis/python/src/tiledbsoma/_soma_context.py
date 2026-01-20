# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import threading
from concurrent.futures import ThreadPoolExecutor
from typing import Literal

from typing_extensions import Self

from . import pytiledbsoma as clib
from ._types import DataProtocol

_Unset = Literal["__unset__"]
_UNSET: _Unset = "__unset__"


class SOMAContext:
    """Maintains TileDB-specific context for TileDB-SOMA objects.

    The context holds TileDB-specific configuration options. Internally, the context contains the storage manager for
    TileDB operations. The internal TileDB context is lazily created with it is first used.

    Most API functions accept an optional `context` keyword argument, but this is typically only necessary in advanced
    cases with multiple contexts per programs. In most cases, the user should set a global default context before
    calling other TileDB-SOMA operations. If an operation is called with no context and a global context is not already
    set, then a global context will be created with default configuration options.

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_handle", "_lock", "threadpool")

    _default_context: Self | None = None

    @classmethod
    def set_default(
        cls,
        config: dict[str, str | float] | None = None,
        threadpool: ThreadPoolExecutor | None = None,
        replace: bool = False,
    ) -> Self:
        """Initializes and returns a global default context to use in SOMA operations.

        This method should be called once at the beginning of your session before opening any SOMA objects
        if you want to customize the TileDB context parameters that will apply to all subsequent operations.
        Otherwise, a default context will be created automatically with standard parameters when you first open
        a SOMA object.

        If the global context was already set, an error will be raised unless ``replace=True``. Setting a new
        global default context will not change the context for TileDB-SOMA objects that were already created.

        Args:
            config: A dictionary of TileDB configuration options to use, overriding the default configuration.

            threadpool: A threadpool to use for concurrent operations. If not provided, a new ThreadPoolExecutor will
                be created with default settings.

            replace: Allow an existing global default context to be replaced.

        Returns:
                The global default context.

        Lifecycle:
            Experimental
        """
        if not replace and cls._default_context is not None:
            raise RuntimeError(
                "A default context was already created. To replace the default context for new objects call this method"
                " again with `replace=True`."
            )
        cls._default_context = cls.create(config=config, threadpool=threadpool)
        return cls._default_context

    @classmethod
    def get_default(cls) -> Self:
        """Returns the current default context used by TileDB-SOMA operations.

        This function returns the context that was either:

        Raise:
            RuntimeError:
                If no default context is set.

        If no global default context is set, a new context will be created with default configuration options.

        Returns:
            The default

        """
        if cls._default_context is None:
            raise RuntimeError(
                "No default context is set. Call `SOMAContext.set_default(...)` to initialize the context."
            )
        return cls._default_context

    @classmethod
    def has_default(cls) -> bool:
        """Returns if the default context is set."""
        return cls._default_context is not None

    @classmethod
    def create(
        cls,
        config: dict[str, str | float] | None = None,
        threadpool: ThreadPoolExecutor | None = None,
    ) -> Self:
        """Create a new SOMAContext.

        Args:
            config: A dictionary of TileDB configuration options to use, overriding the default configuration.

            threadpool: A threadpool to use for concurrent operations. If not provided, a new ThreadPoolExecutor
                will be created with default settings.

        Lifecycle:
            Maturing.
        """
        config = {} if config is None else {key: str(val) for key, val in config.items()}
        config.setdefault("sm.mem.reader.sparse_global_order.ratio_array_data", "0.3")
        return cls(
            clib.SOMAContext(config),
            threadpool=threadpool,
            _dont_call_this_use_create_instead="tiledbsoma-internal-code",
        )

    def __init__(
        self,
        handle: clib.SOMAContext,
        *,
        threadpool: ThreadPoolExecutor | None = None,
        _dont_call_this_use_create_instead: str = "unset",
    ) -> None:
        """Internal-only initializer steps.

        This function is internal; users should create a TileDB-SOMA context using `meth`:create:.
        """
        if _dont_call_this_use_create_instead != "tiledbsoma-internal-code":
            name = type(self).__name__
            raise RuntimeError(
                f"Directly calling `{name}(...)` is intended for TileDB-SOMA internal use only. Use `name.create(...)` "
                f"create a {name}."
            )

        self._handle = handle
        """Internal clib.SOMAContext."""

        self._lock = threading.Lock()
        """A lock to ensure single initialization of ``_tiledb_ctx``."""

        self.threadpool = threadpool or ThreadPoolExecutor()
        """User specified threadpool. If None, we'll instantiate one ourselves."""

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.config}, {self.threadpool})"

    @property
    def config(self) -> dict[str, str]:
        """The TileDB configuration dictionary for this SOMA context.

        If the internal context is loaded, this will be the configuration settings from the context. Before loading
        it is the input configuration properties.

        Lifecycle:
            Maturing.
        """
        with self._lock:
            return dict(self._handle.config())

    def replace(
        self,
        *,
        config: dict[str, str | float | None] | None = None,
        threadpool: ThreadPoolExecutor | _Unset | None = _UNSET,
    ) -> Self:
        """Create a copy of the context, merging changes.

        Args:
            config: A dictionary of TileDB configuration options. To remove parameters from the existing config set the option
                to have a value of ``None``.
            threadpool:
                A threadpool to replace the current threadpool with.

        Lifecycle:
            Maturing.
        """
        with self._lock:
            new_config: dict[str, str | float] = dict(self._handle.config())
            if config is not None:
                for key, val in config.items():
                    if val is None:
                        new_config.pop(key, None)
                    else:
                        new_config[key] = str(val)
            if threadpool == _UNSET:
                threadpool = self.threadpool

        return type(self).create(config=new_config, threadpool=threadpool)

    def data_protocol(self, uri: str) -> DataProtocol:
        """Return the data protocol in use for this URI and context.

        Return value will be a data model identifier. Currently one of:
        * `tiledbv2` - the legacy data model, supported on all storage platforms except Carrara
        * `tiledbv3` - the new, and currently Carrara-specific, data model.

        Args:
            uri:
                An object URI

        Returns:
            The protocol identifier, currently one of `tiledbv2` or `tiledbv3`

        Lifecycle:
            Experimental.
        """
        with self._lock:
            protocol: DataProtocol = self._handle.data_protocol(uri)
        return protocol

    def is_tiledbv2_uri(self, uri: str) -> bool:
        """Return True if the URI will use `tiledbv2` semantics.

        Lifecycle:
            Experimental.
        """
        return self.data_protocol(uri) == "tiledbv2"

    def is_tiledbv3_uri(self, uri: str) -> bool:
        """Return True if the URI will use `tiledbv3` semantics.

        Lifecycle:
            Experimental.

        """
        return self.data_protocol(uri) == "tiledbv3"
