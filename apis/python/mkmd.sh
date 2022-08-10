mkdir -p doc/v1
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
    v1/general_utilities \
    v1/io \
    v1/logging \
    v1/soma_collection \
    v1/soma_dataframe \
    v1/soma_dense_nd_array \
    v1/soma_experiment \
    v1/soma_indexed_dataframe \
    v1/soma_measurement \
    v1/soma_metadata_mapping \
    v1/soma_sparse_nd_array \
    v1/test_general_utilities \
    v1/tiledb_array \
    v1/tiledb_object \
    v1/tiledb_platform_config \
    v1/util \
    v1/util_ann \
    v1/util_arrow \
    v1/util_pandas \
    v1/util_tiledb
do
    echo $x
    pydoc-markdown -I src -m tiledbsc.$x > doc/$x.md
done
