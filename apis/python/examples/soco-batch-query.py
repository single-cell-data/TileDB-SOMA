#!/usr/bin/env python

"""
Batch-process the mean of X/data grouping by obs['cell_type_ontology_term_id']
"""

import tiledbsc
import tiledbsc.util
import sys

import pandas as pd  # so we can type it either way

if len(sys.argv) == 2:
    soco_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
    sys.exit(1)

soco = tiledbsc.SOMACollection(soco_path)

var_ids_column = []
ctot_ids_column = []
means_column = []

ctot_ids = soco.find_unique_obs_values("cell_type_ontology_term_id")
n = len(ctot_ids)
for i, ctot_id in enumerate(ctot_ids):
    soma_slice = soco.attribute_filter(
        obs_attr_names=["cell_type_ontology_term_id"],
        obs_query_string=f'cell_type_ontology_term_id == "{ctot_id}"',
    )
    if soma_slice is None:
        continue

    slice_means = soma_slice.X.mean(axis=0)
    j = 0
    for var_id in soma_slice.var.index:
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
dense_means_df = tiledbsc.util.triples_to_dense_df(sparse_means_df)
print(dense_means_df)

# Remove all-zeroes rows
trimmed_dense_means_df = dense_means_df.loc[(dense_means_df != 0).any(axis=1)]
print()
print("Trimmed of all-zeroes rows:")
print(trimmed_dense_means_df)
