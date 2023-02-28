from typing import Any, Dict, Optional, Tuple, Union

import attrs
import tiledb
from typing_extensions import Self


def _build_default_tiledb_ctx() -> tiledb.Ctx:
    """
    Build a TileDB context starting with reasonable defaults, and overriding and updating with user-provided config
    options.
    """

    # Note: Defaults must provide positive out-of-the-box UX!

    cfg: Dict[str, Union[str, float]] = {
        "sm.mem.reader.sparse_global_order.ratio_array_data": 0.3
    }

    return tiledb.Ctx(cfg)


@attrs.define(frozen=True, kw_only=True)
class SOMATileDBContext:
    """
    Maintains TileDB-specific context for TileDbObjects. This context can be shared across multiple SOMA objects,
    including having a child object inherit it from its parent.

    [lifecycle: experimental]
    """

    tiledb_ctx: tiledb.Ctx = _build_default_tiledb_ctx()

    timestamp: Optional[int] = attrs.field(default=None)
    """
    Timestamp for operations on SOMA objects, in milliseconds since Unix epoch.

    ``None``, the default, does not provide cross-object timestamp consistency:

    - **Reads** will see data as of the time that the individual TileDB object
      is opened.
    - **Writes** will be recorded as of when the TileDB object is closed.

    If a value is passed, that timestamp (representing milliseconds since
    the Unix epoch) is used as the timestamp to record all operations.

    Set to 0xFFFFFFFFFFFFFFFF (UINT64_MAX) to get the absolute latest revision
    (i.e., including changes that occur "after" the current wall time) as of
    when *each* object is opened.
    """

    timestamp_start: int = 0
    """
    Timestamp range start for operations on SOMA objects. This is usually zero
    except for specific, unusual query requirements. This has no effect on
    timestamps used in writing.
    """

    @timestamp.validator
    def _validate_timestamps(self, _: Any, __: Any) -> None:
        if self.timestamp is None:
            if self.timestamp_start:
                raise ValueError(
                    "SOMATileDBContext: if the current `timestamp` is None,"
                    " the `timestamp_start` cannot be set."
                )
        elif not 0 <= self.timestamp_start <= self.timestamp:
            raise ValueError("SOMATileDBContext: invalid read timestamp range")

    _group_tiledb_ctx_: Optional[tiledb.Ctx] = attrs.field(init=False, default=None)
    """Cache for the context used to open Groups."""

    @property
    def _group_tiledb_ctx(self) -> tiledb.Ctx:
        """Internal-only context specifically for Group operations.

        Unlike Arrays, Groups do not take a ``timestamp`` argument and need
        their timestamps set in their context. This should go away if/when
        :class:`tiledb.Group` starts taking a timestamp argument like
        :func:`tiledb.open`.
        """
        if self._group_tiledb_ctx_:
            return self._group_tiledb_ctx_
        if self.timestamp is None:
            ctx = self.tiledb_ctx
        else:
            ctx = group_timestamp_ctx(
                self.tiledb_ctx,
                timestamp_start=self.timestamp_start,
                timestamp=self.timestamp,
            )
        object.__setattr__(self, "_group_tiledb_ctx_", ctx)
        return ctx

    def _timestamp_arg(self) -> Optional[Tuple[int, int]]:
        """The value to use as the ``timestamp`` argument to ``tiledb.open``."""
        if self.timestamp is None:
            return None
        return (self.timestamp_start, self.timestamp)

    def replace(
        self, *, tiledb_config: Optional[Dict[str, Any]] = None, **changes: Any
    ) -> Self:
        """
        Create a copy of the context, merging changes.

        Parameters
        ----------
        tiledb_config - Dict[str, Any]
            a dictionary of parameters for tiledb.Config()

        changes - Any
            Any other parameters will be passed to the class __init__.

        Examples
        --------
        >>> context.replace(timestamp=0)

        >>> context.replace(tiledb_config={"vfs.s3.region": "us-east-2"})
        """
        if tiledb_config:
            new_config = self.tiledb_ctx.config()
            new_config.update(tiledb_config)
            changes["tiledb_ctx"] = tiledb.Ctx(config=new_config)
        return attrs.evolve(self, **changes)


def group_timestamp_ctx(
    ctx: tiledb.Ctx, *, timestamp_start: int, timestamp: int
) -> tiledb.Ctx:
    """Builds a TileDB context to open groups at the given timestamp."""
    group_cfg = ctx.config().dict()
    group_cfg["sm.group.timestamp_start"] = timestamp_start
    group_cfg["sm.group.timestamp_end"] = timestamp
    return tiledb.Ctx(group_cfg)
