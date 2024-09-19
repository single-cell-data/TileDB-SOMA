# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

from __future__ import annotations

import datetime
import functools
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, Literal, Mapping, Optional, Union

from somacore import ContextBase
from typing_extensions import Self

import tiledb

from .. import pytiledbsoma as clib
from .._general_utilities import assert_version_before
from .._types import OpenTimestamp
from .._util import ms_to_datetime, to_timestamp_ms


def _warn_ctx_deprecation() -> None:
    assert_version_before(1, 15)
    warnings.warn(
        "tiledb_ctx is now deprecated for removal in 1.15. "
        "Use tiledb_config instead by passing "
        "SOMATileDBContext(tiledb_config=ctx.config().dict()).",
        DeprecationWarning,
        stacklevel=3,
    )


def _default_config(
    override: Mapping[str, Union[str, float]]
) -> Dict[str, Union[str, float]]:
    """Returns a fresh dictionary with TileDB config values.

    These should be reasonable defaults that can be used out-of-the-box.
    ``override`` does exactly what it says: overrides default entries.
    """
    cfg: Dict[str, Union[str, float]] = {
        "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3
    }
    cfg.update(override)
    return cfg


@functools.lru_cache(maxsize=None)
def _default_global_ctx() -> tiledb.Ctx:
    """Lazily builds a default TileDB Context with the default config."""
    return tiledb.Ctx(_default_config({}))


def _maybe_timestamp_ms(input: Optional[OpenTimestamp]) -> Optional[int]:
    if input is None:
        return None
    return to_timestamp_ms(input)


_Unset = Literal["__unset__"]
_UNSET: _Unset = "__unset__"


class SOMATileDBContext(ContextBase):
    """Maintains TileDB-specific context for TileDB-SOMA objects.
    This context can be shared across multiple objects,
    including having a child object inherit it from its parent.

    Treat this as immutable. Only code internal to the ``_soma_tiledb_context``
    module should be modifying this. Use the ``replace`` method to construct
    a new ``SOMATileDBContext`` with new values.

    Lifecycle:
        Maturing.
    """

    def __init__(
        self,
        tiledb_ctx: Optional[tiledb.Ctx] = None,
        tiledb_config: Optional[Dict[str, Union[str, float]]] = None,
        timestamp: Optional[OpenTimestamp] = None,
        threadpool: Optional[ThreadPoolExecutor] = None,
    ) -> None:
        """Initializes a new SOMATileDBContext.

        Either ``tiledb_ctx`` or ``tiledb_config`` may be provided, or both may
        be left at their default. If neither are provided, this will use a
        single shared :class:`tiledb.Ctx` instantiated upon first use.
        If ``tiledb_ctx`` is provided, that exact :class:`tiledb.Ctx` is used.
        If a ``tiledb_config`` is provided (in the form of a ``dict``),
        it is used to construct a new ``Ctx``.

        Args:
            tiledb_ctx: An existing TileDB Context for use as the TileDB Context
                in this configuration.

            tiledb_config: A set of TileDB configuration options to use,
                overriding the default configuration.

            timestamp: The default timestamp for operations on SOMA objects,
                provided either as a ``datetime.datetime`` or a number of
                milliseconds since the Unix epoch.

                WARNING: This should not be set unless you are *absolutely* sure
                you want to use the same timestamp across multiple operations.
                If multiple writes to the same object are performed at the same
                timestamp, they have no defined order. In almost all cases,
                it is better to pass a timestamp to a single ``open`` call,
                or to simply use the default behavior.

                This is used when a timestamp is not provided to
                an ``open`` operation.

                ``None``, the default, sets the timestamp on each root ``open``
                operation. That is, if you ``open`` a collection, and access
                individual members of the collection through indexing or
                ``add_new``, the timestamp of all of those operations will be
                that of the time you called ``open``.

                If a value is passed, that timestamp is used as the timestamp
                to record all operations.

                Set to 0xFFFFFFFFFFFFFFFF (UINT64_MAX) to get the absolute
                latest revision (i.e., including changes that occur "after"
                the current wall time) as of when *each* object is opened.

            threadpool: A threadpool to use for concurrent operations. If not
                provided, a new ThreadPoolExecutor will be created with
                default settings.
        """
        if tiledb_ctx is not None:
            _warn_ctx_deprecation()

        if tiledb_ctx is not None and tiledb_config is not None:
            raise ValueError(
                "only one of tiledb_ctx or tiledb_config"
                " may be set when constructing a SOMATileDBContext"
            )
        self._lock = threading.Lock()
        """A lock to ensure single initialization of ``_tiledb_ctx``."""
        self._initial_config = (
            None if tiledb_config is None else _default_config(tiledb_config)
        )

        """A dictionary of options to override the default TileDB config.

        This includes both the user-provided options and the default options
        that we provide to TileDB. If this is unset, then either we were
        provided with a TileDB Ctx, or we need to use The Default Global Ctx.
        """
        self._tiledb_ctx = tiledb_ctx
        """The TileDB context to use, either provided or lazily constructed."""
        self._timestamp_ms = _maybe_timestamp_ms(timestamp)

        self.threadpool = threadpool or ThreadPoolExecutor()
        """User specified threadpool. If None, we'll instantiate one ourselves."""
        self._native_context: Optional[clib.SOMAContext] = None
        """Lazily construct clib.SOMAContext."""

    @property
    def timestamp_ms(self) -> Optional[int]:
        """
        The default timestamp for SOMA operations, as milliseconds since
        the Unix epoch, or ``None`` if not provided.
        """
        return self._timestamp_ms

    @property
    def timestamp(self) -> Optional[datetime.datetime]:
        """
        The default timestamp for SOMA operations, or ``None`` if not provided.
        """
        if self.timestamp_ms is None:
            return None
        return ms_to_datetime(self.timestamp_ms)

    @property
    def native_context(self) -> clib.SOMAContext:
        """The C++ SOMAContext for this SOMA context."""
        with self._lock:
            if self._native_context is None:
                cfg = self._internal_tiledb_config()
                self._native_context = clib.SOMAContext(
                    {k: str(v) for k, v in cfg.items()}
                )
            return self._native_context

    @property
    def tiledb_ctx(self) -> tiledb.Ctx:
        """The TileDB-Py Context for this SOMA context."""
        _warn_ctx_deprecation()

        with self._lock:
            if self._tiledb_ctx is None:
                if self._initial_config is None:
                    # Special case: we need to use the One Global Default.
                    self._tiledb_ctx = _default_global_ctx()
                else:
                    self._tiledb_ctx = tiledb.Ctx(self._initial_config)
            return self._tiledb_ctx

    @property
    def tiledb_config(self) -> Dict[str, Union[str, float]]:
        """The TileDB configuration dictionary for this SOMA context.

        If this ``SOMATileDBContext`` already has a ``tiledb_ctx``, this will
        return the full set of values from that TileDB Context; otherwise, this
        will only return the values that will be passed into the ``tiledb.Ctx``
        constructor, including both SOMA defaults and the ``tiledb_config``
        parameter passed into this object's constructor.

        This always returns a fresh dictionary.
        """
        with self._lock:
            return self._internal_tiledb_config()

    def _internal_tiledb_config(self) -> Dict[str, Union[str, float]]:
        """Internal function for getting the TileDB Config.

        Returns a new dict with the contents. Caller must hold ``_lock``.
        """
        # We have a clib.SOMAContext. Return its actual config.
        if self._native_context is not None:
            return dict(self._native_context.config())

        # We have TileDB Context. Return its actual config.
        # TODO This block will be deleted once tiledb_ctx is removed in 1.15
        if self._tiledb_ctx is not None:
            return dict(self._tiledb_ctx.config())

        # Our context has not yet been built.
        # We return what will be passed into the context.
        return (
            dict(self._initial_config)
            if self._initial_config is not None
            else _default_config({})
        )

    def replace(
        self,
        *,
        tiledb_config: Optional[Dict[str, Any]] = None,
        tiledb_ctx: Optional[tiledb.Ctx] = None,
        timestamp: Union[None, OpenTimestamp, _Unset] = _UNSET,
        threadpool: Union[None, ThreadPoolExecutor, _Unset] = _UNSET,
    ) -> Self:
        """Create a copy of the context, merging changes.

        Args:
            tiledb_config:
                A dictionary of parameters for `tiledb.Config() <https://tiledb-inc-tiledb.readthedocs-hosted.com/projects/tiledb-py/en/stable/python-api.html#config>`_.
                To remove a parameter from the existing config, provide ``None``
                as the value.
            tiledb_ctx:
                A TileDB Context to replace the current context with.
            timestamp:
                A timestamp to replace the current timestamp with.
                Explicitly passing ``None`` will remove the timestamp.
                For details, see the description of ``timestamp``
                in :meth:`__init__`.
            threadpool:
                A threadpool to replace the current threadpool with.

        Lifecycle:
            Maturing.

        Examples:
            >>> context.replace(timestamp=1_512_658_800_000)  # UNIX millis
            >>> new_region_context = context.replace(
            ...     tiledb_config={"vfs.s3.region": "us-east-2"})
            >>> back_to_default_context = new_region_context.replace(
            ...     tiledb_config={"vfs.s3.region": None})
        """
        with self._lock:
            if tiledb_ctx is not None:
                _warn_ctx_deprecation()

            if tiledb_config is not None:
                if tiledb_ctx:
                    raise ValueError(
                        "Either tiledb_config or tiledb_ctx may be provided"
                        " to replace(), but not both."
                    )
                new_config = self._internal_tiledb_config()
                new_config.update(tiledb_config)
                tiledb_config = {k: v for (k, v) in new_config.items() if v is not None}

            if timestamp == _UNSET:
                # Keep the existing timestamp if not overridden.
                timestamp = self._timestamp_ms
            if threadpool == _UNSET:
                # Keep the existing threadpool if not overridden.
                threadpool = self.threadpool

        assert timestamp is None or isinstance(timestamp, (datetime.datetime, int))
        return type(self)(
            tiledb_config=tiledb_config,
            tiledb_ctx=tiledb_ctx,
            timestamp=timestamp,
            threadpool=threadpool,
        )

    def _open_timestamp_ms(self, in_timestamp: Optional[OpenTimestamp]) -> int:
        """Returns the real timestamp that should be used to open an object."""
        if in_timestamp is not None:
            return to_timestamp_ms(in_timestamp)
        if self.timestamp_ms is not None:
            return self.timestamp_ms
        return int(time.time() * 1000)


def _validate_soma_tiledb_context(context: Any) -> SOMATileDBContext:
    """Returns the argument, as long as it's a ``SOMATileDBContext``, or a new
    one if the argument is ``None``. While we already have static type-checking,
    a few things are extra-important to have runtime validation on.  Since it's
    easy for users to pass a ``tiledb.Ctx`` when a ``SOMATileDBContext`` is
    expected, we should offer a helpful redirect when they do.
    """

    if context is None:
        return SOMATileDBContext()

    if isinstance(context, tiledb.Ctx):
        raise TypeError(
            "context is a tiledb.Ctx, not a SOMATileDBContext -- please wrap it in tiledbsoma.SOMATileDBContext(...)"
        )

    if not isinstance(context, SOMATileDBContext):
        raise TypeError("context is not a SOMATileDBContext")

    return context
