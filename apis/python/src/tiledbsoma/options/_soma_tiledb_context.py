# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

import datetime
import time
from typing import Any, Dict, Optional, Union

import attrs
import tiledb
from typing_extensions import Self

from .._types import OpenTimestamp
from .._util import ms_to_datetime, to_timestamp_ms


def _build_default_tiledb_ctx() -> tiledb.Ctx:
    """Build a TileDB context starting with reasonable defaults,
    and overriding and updating with user-provided config options.
    """

    # Note: Defaults must provide positive out-of-the-box UX!

    cfg: Dict[str, Union[str, float]] = {
        "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3
    }

    return tiledb.Ctx(cfg)


def _maybe_timestamp_ms(input: Optional[OpenTimestamp]) -> Optional[int]:
    if input is None:
        return None
    return to_timestamp_ms(input)


@attrs.define(frozen=True, kw_only=True)
class SOMATileDBContext:
    """Maintains TileDB-specific context for TileDB-SOMA objects.
    This context can be shared across multiple objects,
    including having a child object inherit it from its parent.

    Lifecycle:
        Experimental.
    """

    tiledb_ctx: tiledb.Ctx = _build_default_tiledb_ctx()

    timestamp_ms: Optional[int] = attrs.field(
        default=None, converter=_maybe_timestamp_ms, alias="timestamp"
    )
    """
    Default timestamp for operations on SOMA objects, in milliseconds since the Unix epoch.

    WARNING: This should not be set unless you are *absolutely* sure you want to
    use the same timestamp across multiple operations. If multiple writes to the
    same object are performed at the same timestamp, they have no defined order.
    In most cases, it is better to pass a timestamp to a single ``open`` call,
    or to simply use the default behavior.

    This is used when a timestamp is not provided to an ``open`` operation.

    ``None``, the default, sets the timestamp on each root ``open`` operation.
    That is, if you ``open`` a collection, and access individual members of the
    collection through indexing or ``add_new``, the timestamp of all of those
    operations will be that of the time you called ``open``.

    If a value is passed, that timestamp (representing milliseconds since
    the Unix epoch) is used as the timestamp to record all operations.

    Set to 0xFFFFFFFFFFFFFFFF (UINT64_MAX) to get the absolute latest revision
    (i.e., including changes that occur "after" the current wall time) as of
    when *each* object is opened.
    """

    @property
    def timestamp(self) -> Optional[datetime.datetime]:
        if self.timestamp_ms is None:
            return None
        return ms_to_datetime(self.timestamp_ms)

    def replace(
        self, *, tiledb_config: Optional[Dict[str, Any]] = None, **changes: Any
    ) -> Self:
        """Create a copy of the context, merging changes.

        Args:
            tiledb_config:
                A dictionary of parameters for `tiledb.Config() <https://tiledb-inc-tiledb.readthedocs-hosted.com/projects/tiledb-py/en/stable/python-api.html#config>`_.
            changes:
                Any other parameters will be passed to the class ``__init__``.

        Lifecycle:
            Experimental.

        Examples:
            >>> context.replace(timestamp=0)

            >>> context.replace(tiledb_config={"vfs.s3.region": "us-east-2"})
        """
        if tiledb_config:
            new_config = self.tiledb_ctx.config()
            new_config.update(tiledb_config)
            changes["tiledb_ctx"] = tiledb.Ctx(config=new_config)
        return attrs.evolve(self, **changes)

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
