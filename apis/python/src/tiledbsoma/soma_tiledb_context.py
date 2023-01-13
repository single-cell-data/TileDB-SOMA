import os
from typing import Optional, cast

import attrs
import tiledb


def build_default_tiledb_ctx() -> tiledb.Ctx:
    """
    Build a TileDB context starting with reasonable defaults, and overriding and updating with user-provided config
    options.
    """

    # Note: Defaults must provide positive out-of-the-box UX!

    cfg = {"sm.mem.reader.sparse_global_order.ratio_array_data": 0.3}

    # This is necessary for smaller tile capacities when querying with a smaller memory budget.

    # Temp workaround pending https://app.shortcut.com/tiledb-inc/story/23827
    region = os.getenv("AWS_DEFAULT_REGION")
    if region is not None:
        cfg["vfs.s3.region"] = cast(str, region)  # type: ignore

    return tiledb.Ctx(cfg)


@attrs.define(frozen=True, kw_only=True)
class SOMATileDBContext:
    """
    Maintains TileDB-specific context for TileDbObjects. This context can be shared across multiple SOMA objects,
    including having a child object inherit it from its parent.
    """

    tiledb_ctx: tiledb.Ctx = build_default_tiledb_ctx()

    member_uris_are_relative: Optional[bool] = None
    """Allows "relocatability" for local disk / S3, and correct behavior for TileDB Cloud."""

    # read_timestamp, e.g.
