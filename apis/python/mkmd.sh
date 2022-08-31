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
    util_tiledb
do
    echo $x
    pydoc-markdown -I src -m tiledbsc.$x > doc/$x.md
done
