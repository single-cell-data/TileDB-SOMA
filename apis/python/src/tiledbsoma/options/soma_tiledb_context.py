import time
from typing import Dict, Optional, Union

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

    read_timestamp: int = -1
    """
    Timestamp for operations on SOMA objects open in read mode, in milliseconds since Unix epoch.
    Defaults to the time of context initialization. Set to 0xFFFFFFFFFFFFFFFF (UINT64_MAX) to get
    the latest revision as of when *each* object is opened.

    SOMA objects opened in *write* mode ignore any read timestamp.
    """

    read_timestamp_start: int = 0
    """
    Timestamp range start for operations on SOMA objects opened in read mode. This is usually zero
    except for specific, unusual query requirements.
    """

    write_timestamp: Optional[int] = None
    """
    Timestamp applied to all SOMA object write operations. If unset, each individual write
    operation receives the timestamp as of when the operation executes.

    Caution: overlapping writes (of overlapping array slices, or setting the same collection key)
    should be avoided when write_timestamp is set. Distinct, overlapping write operations given the
    same timestamp may be applied in any order.
    """

    # (internal) tiledb.Ctx specifically for tiledb.Group operations; unlike arrays, tiledb.Group
    # needs its timestamps set in the tiledb.Ctx.
    _group_read_tiledb_ctx: tiledb.Ctx = attrs.field(init=False)
    _group_write_tiledb_ctx: tiledb.Ctx = attrs.field(init=False)

    def __attrs_post_init__(self) -> None:
        "initialization hook invoked by the attrs-generated __init__"

        if self.read_timestamp < 0:
            object.__setattr__(self, "read_timestamp", int(time.time() * 1000))
        assert (
            self.read_timestamp_start >= 0
            and self.read_timestamp >= self.read_timestamp_start
        )
        assert self.write_timestamp is None or self.write_timestamp >= 0

        # prepare the pair of tiledb.Ctx for groups
        group_read_config = self.tiledb_ctx.config().dict()
        group_read_config["sm.group.timestamp_start"] = self.read_timestamp_start
        group_read_config["sm.group.timestamp_end"] = self.read_timestamp
        object.__setattr__(
            self, "_group_read_tiledb_ctx", tiledb.Ctx(group_read_config)
        )
        assert isinstance(self._group_read_tiledb_ctx, tiledb.Ctx)

        if self.write_timestamp is not None:
            group_write_config = self.tiledb_ctx.config().dict()
            group_write_config["sm.group.timestamp_start"] = self.write_timestamp
            group_write_config["sm.group.timestamp_end"] = self.write_timestamp
            object.__setattr__(
                self, "_group_write_tiledb_ctx", tiledb.Ctx(group_write_config)
            )
            assert isinstance(self._group_write_tiledb_ctx, tiledb.Ctx)
        else:
            object.__setattr__(self, "_group_write_tiledb_ctx", self.tiledb_ctx)
