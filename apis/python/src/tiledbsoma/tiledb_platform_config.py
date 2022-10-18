from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class TileDBPlatformConfig:
    """
    A place to put configuration options various users may wish to change.  These are mainly TileDB array-schema parameters.
    """

    # TODO: pending further work on
    # https://github.com/single-cell-data/TileDB-SOMA/issues/27
    obs_extent: int = 256
    var_extent: int = 2048
    X_capacity: int = 100000
    X_tile_order: str = "row-major"
    X_cell_order: str = "row-major"
    # https://github.com/single-cell-data/TileDB-SOMA/issues/27
    string_dim_zstd_level: int = 3
    write_X_chunked: bool = True
    goal_chunk_nnz: int = 200_000_000
    # Allows relocatability for local disk / S3, and correct behavior for TileDB Cloud
    member_uris_are_relative: Optional[bool] = None

    max_thread_pool_workers: int = 8

    # Temporary
    #
    # Feature flag controlling some test features in soma.io. If False,
    # uses original tiledb.from_pandas code to create axis dataframes
    # in soma.io. If True, uses the core SOMA classes to do so.
    #
    # Feature flag in place to allow for additional compatibility testing.
    # TODO: remove when SOMA/Arrow path is known to be working well.
    #
    from_anndata_write_pandas_using_arrow: bool = True
