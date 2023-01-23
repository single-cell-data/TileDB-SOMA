#!/usr/bin/env python

"""
Batch-process the mean of X/data grouping by obs['cell_type_ontology_term_id']
"""

import sys

import pandas as pd
import tiledb

import tiledbsoma
import tiledbsoma.util

if len(sys.argv) == 2:
    soco_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
    sys.exit(1)

soco = tiledbsoma.SOMACollection(soco_path)

# per-column buffer size
ctx = tiledb.Ctx({"py.init_buffer_bytes": 4 * 1024**3})

var_ids_column = []
ctot_ids_column = []
means_column = []

ctot_ids = soco.find_unique_obs_values("cell_type_ontology_term_id")
n = len(ctot_ids)
print("cell_type_ontology_term_id count =", n)
for i, ctot_id in enumerate(ctot_ids):
    soma_slices = soco.query(
        obs_attrs=["cell_type_ontology_term_id"],
        obs_query_string=f'cell_type_ontology_term_id == {ctot_id!r}',
    )
    if soma_slices == []:
        continue

    result_soma_slice = tiledbsoma.SOMASlice.concat(soma_slices)

    if result_soma_slice is not None:

        slice_means = result_soma_slice.X["data"].mean(axis=0)
        j = 0
        for var_id in result_soma_slice.var.index:
            value = slice_means[0, j]
            var_ids_column.append(var_id)
            ctot_ids_column.append(ctot_id)
            means_column.append(value)
            j += 1
        print("... %6.2f%% done %s" % (100 * i / n, ctot_id))

sparse_means_df = pd.DataFrame(
    {
        "var_id": var_ids_column,
        "cell_type_ontology_term_id": ctot_ids_column,
        "value": means_column,
    }
)

sparse_means_df.set_index(["var_id", "cell_type_ontology_term_id"], inplace=True)
print(sparse_means_df)

# This is in 'triples' format with three columns var_id, cell_type_ontology_term_id, and value.
# Convert it to a dataframe with var_id row labels, cell_type_ontology_term_id column labels, and value data.
print()
print("As dense with zero-fill:")
dense_means_df = tiledbsoma.util.triples_to_dense_df(sparse_means_df)
print(dense_means_df)

# Remove all-zeroes rows
trimmed_dense_means_df = dense_means_df.loc[(dense_means_df != 0).any(axis=1)]
print()
print("Trimmed of all-zeroes rows:")
print(trimmed_dense_means_df)
