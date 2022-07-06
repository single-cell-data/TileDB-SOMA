from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class SOMAOptions:
    """
    A place to put configuration options various users may wish to change.
    These are mainly TileDB array-schema parameters.
    """

    # TODO: pending further work on
    # https://github.com/single-cell-data/TileDB-SingleCell/issues/27
    obs_extent: int = 256
    var_extent: int = 2048
    X_capacity: int = 100000
    X_tile_order: str = "row-major"
    X_cell_order: str = "row-major"
    # https://github.com/single-cell-data/TileDB-SingleCell/issues/27
    string_dim_zstd_level: int = 22
    write_X_chunked: bool = True
    goal_chunk_nnz: int = 20_000_000
    # Allows relocatability for local disk / S3, and correct behavior for TileDB Cloud
    member_uris_are_relative: Optional[bool] = None
    allows_duplicates: bool = False
