mkdir -p doc/tiledbsoma
export PATH=${PATH}:${HOME}/.local/bin
for x in \
    annotation_dataframe \
    annotation_matrix \
    annotation_matrix_group \
    annotation_pairwise_matrix_group \
    assay_matrix \
    assay_matrix_group \
    raw_group \
    soma \
    soma_collection \
    soma_options \
    soma_slice \
    tiledb_array \
    tiledb_group \
    tiledb_object \
    uns_array \
    uns_group \
    util \
    util_ann \
    util_tiledb \
    \
    tiledbsoma/general_utilities \
    tiledbsoma/io \
    tiledbsoma/logging \
    tiledbsoma/soma_collection \
    tiledbsoma/soma_dataframe \
    tiledbsoma/soma_dense_nd_array \
    tiledbsoma/soma_experiment \
    tiledbsoma/soma_indexed_dataframe \
    tiledbsoma/soma_measurement \
    tiledbsoma/soma_metadata_mapping \
    tiledbsoma/soma_sparse_nd_array \
    tiledbsoma/test_general_utilities \
    tiledbsoma/tiledb_array \
    tiledbsoma/tiledb_object \
    tiledbsoma/tiledb_platform_config \
    tiledbsoma/util \
    tiledbsoma/util_ann \
    tiledbsoma/util_arrow \
    tiledbsoma/util_pandas \
    tiledbsoma/util_tiledb
do
    echo $x
    pydoc-markdown -I src -m tiledbsc.$x > doc/$x.md
done
