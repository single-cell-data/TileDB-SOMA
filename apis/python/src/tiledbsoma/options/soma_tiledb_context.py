from typing import Dict, Optional, Tuple, Union

import attrs
import tiledb


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

    member_uris_are_relative: Optional[bool] = None
    """Allows "relocatability" for local disk / S3, and correct behavior for TileDB Cloud."""

    read_timestamp: Optional[int] = None
    """
    Timestamp for all array read operations. If unspecified, then always read the latest data.
    The read timestamp must not be set while writing, so that a writer can read its own writes.
    """

    read_timestamp_start: Optional[int] = None
    """
    Timestamp range start for all array read operations.
    This is usually unset (implicitly zero) except for specific, unusual query circumstances.
    """

    def _tiledb_read_timestamp_arg(self) -> Optional[Tuple[int, int]]:
        "(internal) form the read timestamp tuple arg for TileDB methods"
        if self.read_timestamp_start is None and self.read_timestamp is None:
            return None
        start = (
            self.read_timestamp_start if self.read_timestamp_start is not None else 0
        )
        end = (
            self.read_timestamp
            if self.read_timestamp is not None
            else 0xFFFFFFFFFFFFFFFF  # UINT64_MAX
        )
        return (start, end)

    _write_timestamp: Optional[int] = attrs.field(
        alias="_write_timestamp", default=None
    )
    "(internal) override the timestamp applied to all write operations, for testing purposes."
