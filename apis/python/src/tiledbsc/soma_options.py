from typing import List, Dict

# TODO: when UTF-8 attributes are queryable using TileDB-Py's QueryCondition API we can remove this.
# Context: https://github.com/single-cell-data/TileDB-SingleCell/issues/99.
default_col_names_to_store_as_ascii = {
    "obs": [
        "assay_ontology_term_id",
        "sex_ontology_term_id",
        "organism_ontology_term_id",
        "disease_ontology_term_id",
        "ethnicity_ontology_term_id",
        "development_stage_ontology_term_id",
        "cell_type_ontology_term_id",
        "tissue_ontology_term_id",
        "cell_type",
        "assay",
        "disease",
        "organism",
        "sex",
        "tissue",
        "ethnicity",
        "development_stage",
    ],
    "var": [
        "feature_biotype",
        "feature_name",
        "feature_reference",
    ],
}


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
    col_names_to_store_as_ascii: Dict[str, List[str]]
    member_uris_are_relative: bool

    def __init__(
        self,
        obs_extent=256,
        var_extent=2048,
        X_capacity=100000,
        X_tile_order="row-major",
        X_cell_order="row-major",
        string_dim_zstd_level=22,  # https://github.com/single-cell-data/TileDB-SingleCell/issues/27
        write_X_chunked=True,
        goal_chunk_nnz=10000000,
        col_names_to_store_as_ascii=default_col_names_to_store_as_ascii,
        member_uris_are_relative=None,  # Allows relocatability for local disk / S3, and correct behavior for TileDB Cloud
    ):
        self.obs_extent = obs_extent
        self.var_extent = var_extent
        self.X_capacity = X_capacity
        self.X_tile_order = X_tile_order
        self.X_cell_order = X_cell_order
        self.string_dim_zstd_level = string_dim_zstd_level
        self.write_X_chunked = write_X_chunked
        self.goal_chunk_nnz = goal_chunk_nnz
        self.col_names_to_store_as_ascii = col_names_to_store_as_ascii
        self.member_uris_are_relative = member_uris_are_relative
