import os
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

    # This is necessary for smaller tile capacities when querying with a smaller memory budget.

    # Temp workaround pending https://app.shortcut.com/tiledb-inc/story/23827
    region = os.getenv("AWS_DEFAULT_REGION")
    if region is not None:
        cfg["vfs.s3.region"] = region

    return tiledb.Ctx(cfg)


@attrs.define(frozen=True, kw_only=True)
class SOMATileDBContext:
    """
    Maintains TileDB-specific context for TileDbObjects. This context can be shared across multiple SOMA objects,
    including having a child object inherit it from its parent.
    """

    tiledb_ctx: tiledb.Ctx = _build_default_tiledb_ctx()

    member_uris_are_relative: Optional[bool] = None
    """Allows "relocatability" for local disk / S3, and correct behavior for TileDB Cloud."""

    read_timestamp_start: Optional[int] = None
    "Timestamp range start for all array read operations. Usually, implicitly zero."

    read_timestamp_end: Optional[int] = None
    """
    Timestamp range end for all array read operations.
    If unspecified, then always read the latest data.
    When writing SOMA objects, the read timestamp range should usually be left empty, so that the
    writer can read its own writes (and associated metadata).
    """

    write_timestamp: Optional[int] = None
    "Timestamp applied to all array write operations."

    def _tiledb_read_timestamp_arg(self) -> Optional[Tuple[int, int]]:
        "(internal) form the read timestamp tuple arg for TileDB methods"
        if self.read_timestamp_start is None and self.read_timestamp_end is None:
            return None
        start = (
            self.read_timestamp_start if self.read_timestamp_start is not None else 0
        )
        end = (
            self.read_timestamp_end
            if self.read_timestamp_end is not None
            else 0xFFFFFFFFFFFFFFFF  # UINT64_MAX
        )
        return (start, end)
