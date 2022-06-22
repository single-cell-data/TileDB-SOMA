#!/usr/bin/env python

import tiledbsc
import pandas as pd
import sys

if len(sys.argv) == 2:
    soco_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
    sys.exit(1)
soco = tiledbsc.SOMACollection(soco_path)

obs_labels = ["cell_type", "tissue", "cell_type_ontology_term_id"]
var_labels = ["feature_name"]

print("================================================================")
for obs_label in obs_labels:
    print()
    print("----------------------------------------------------------------")
    print("Per-SOMA values for", obs_label)
    for soma in soco:
        print()
        print("SOMA", soma.name)
        print(soma.obs.df().groupby(obs_label).size())

print("================================================================")
for var_label in var_labels:
    print()
    print("----------------------------------------------------------------")
    print("Per-SOMA values for", var_label)
    for soma in soco:
        print()
        print("SOMA", soma.name)
        print(soma.var.df().groupby(var_label).size())

print()
print("================================================================")
for obs_label in obs_labels:
    print()
    print("----------------------------------------------------------------")
    print("Counts of SOMAs having", obs_label)
    print()
    print("obs_label", obs_label)
    unique_values_to_counts = {}
    for soma in soco:
        if obs_label in soma.obs.keys():
            unique_values = soma.obs.df().groupby(obs_label).size().index
            for unique_value in unique_values:
                if unique_value in unique_values_to_counts:
                    unique_values_to_counts[unique_value] += 1
                else:
                    unique_values_to_counts[unique_value] = 1
    df = pd.DataFrame.from_dict(
        {
            "obs_label_value": unique_values_to_counts.keys(),
            "count": unique_values_to_counts.values(),
        }
    )
    df.set_index("obs_label_value", inplace=True)
    print(df)


print()
print("================================================================")
for obs_label in obs_labels:
    print()
    print("----------------------------------------------------------------")
    print("Collection-wide counts of values of", obs_label)
    print()
    print("obs_label", obs_label)
    dfs = []
    for soma in soco:
        if obs_label in soma.obs.keys():
            df = soma.get_obs_value_counts(obs_label)
            dfs.append(df)
    df = pd.concat(dfs)
    print(df)
