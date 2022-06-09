#!/usr/bin/env python
#!/usr/bin/env python -i

import tiledb
import tiledbsc
import sys, os

import anndata
import anndata as ad  # so we can type it either way

import pandas
import pandas as pd  # so we can type it either way
import numpy
import numpy as np  # so we can type it either way
import scipy

from typing import List, Dict

if len(sys.argv) == 1:
    soco_path = "soma-collection"
elif len(sys.argv) == 2:
    soco_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
    sys.exit(1)

soco = tiledbsc.SOMACollection(soco_path)

# Interact at the Python prompt now

# ----------------------------------------------------------------
def time_soco_attribute_filter_prototype(
    soco: tiledbsc.SOMACollection,
    obs_attr_names: List[str],
    obs_query_string: str,
    var_attr_names: List[str],
    var_query_string: str,
) -> None:

    print()
    print(
        "================================================================ TIME ATTRIBUTE-FILTER QUERY"
    )
    s = tiledbsc.util.get_start_stamp()

    soma_slices = []
    for soma in soco:

        # E.g. querying for 'cell_type == "blood"' but this SOMA doesn'tiledbsc have a cell_type column in
        # its obs at all.
        if not soma.obs.has_attr_names(obs_attr_names):
            continue
        # E.g. querying for 'feature_name == "MT-CO3"' but this SOMA doesn'tiledbsc have a feature_name
        # column in its var at all.
        if not soma.var.has_attr_names(var_attr_names):
            continue

        soma_slice = soma.attribute_filter(obs_query_string, var_query_string)
        if soma_slice != None:
            print()
            print("----------------------------------------------------------------")
            print(soma.uri)
            # soma_slice.describe()
            # xxx temp health check: a = soma_slice.to_anndata()
            # xxx temp as another health check try: soma.from_soma_slice()
            # print(soma_slice)
            # print()

            soma_slices.append(soma_slice)

    print(tiledbsc.util.format_elapsed(s, f"attribute-filter prototype query"))

    result_soma_slice = tiledbsc.SOMASlice.concat(soma_slices)
    if result_soma_slice is None:
        print("Empty slice")
    else:
        output_file_name = "first-light.h5ad"
        a = result_soma_slice.to_anndata()
        a.write_h5ad(output_file_name)
        print("Wrote", output_file_name)

        output_soma_path = "first-light"
        soma = tiledbsc.SOMA.from_soma_slice(result_soma_slice, output_soma_path)
        print("Wrote", output_soma_path)


# ================================================================
time_soco_attribute_filter_prototype(
    soco,
    obs_attr_names=["tissue"],
    obs_query_string='tissue == "blood"',
    var_attr_names=["feature_name"],
    var_query_string='feature_name == "MT-CO3"',
)
