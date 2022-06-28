class SOMAOptions:
    """
    A place to put configuration options various users may wish to change.
    These are mainly TileDB array-schema parameters.
    """

    # TODO: pending further work on
    # https://github.com/single-cell-data/TileDB-SingleCell/issues/27
    obs_extent: int
    var_extent: int
    X_capacity: int
    X_tile_order: str
    X_cell_order: str
    string_dim_zstd_level: int
    write_X_chunked: bool
    goal_chunk_nnz: int
    member_uris_are_relative: bool
    allows_duplicates: bool

    def __init__(
        self,
        obs_extent=256,
        var_extent=2048,
        X_capacity=100000,
        X_tile_order="row-major",
        X_cell_order="row-major",
        string_dim_zstd_level=22,  # https://github.com/single-cell-data/TileDB-SingleCell/issues/27
        write_X_chunked=True,
        goal_chunk_nnz=20_000_000,
        member_uris_are_relative=None,  # Allows relocatability for local disk / S3, and correct behavior for TileDB Cloud
        allows_duplicates=False,
    ):
        self.obs_extent = obs_extent
        self.var_extent = var_extent
        self.X_capacity = X_capacity
        self.X_tile_order = X_tile_order
        self.X_cell_order = X_cell_order
        self.string_dim_zstd_level = string_dim_zstd_level
        self.write_X_chunked = write_X_chunked
        self.goal_chunk_nnz = goal_chunk_nnz
        self.member_uris_are_relative = member_uris_are_relative
        self.allows_duplicates = allows_duplicates
