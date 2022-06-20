#!/usr/bin/env python

import tiledbsc
import sys

if len(sys.argv) == 2:
    soco_path = sys.argv[1]
else:
    print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
    sys.exit(1)
soco = tiledbsc.SOMACollection(soco_path)

obs_labels = ["cell_type", "tissue", "cell_type_ontology_term_id"]

print("================================================================")
for obs_label in obs_labels:
    print()
    print("----------------------------------------------------------------")
    print("Per-SOMA values for", obs_label)
    for soma in soco:
        print()
        print("SOMA", soma.name)
        print(soma.get_obs_value_counts(obs_label))

print()
print("================================================================")
for obs_label in obs_labels:
    print()
    print("----------------------------------------------------------------")
    print("Counts of SOMAs having", obs_label)
    print()
    print("obs_label", obs_label)
    df = soco.get_obs_value_counts(obs_label, False)
    print(df)

print()
print("================================================================")
for obs_label in obs_labels:
    print()
    print("----------------------------------------------------------------")
    print("Collection-wide counts of values of", obs_label)
    print()
    print("obs_label", obs_label)
    df = soco.get_obs_value_counts(obs_label, True)
    print(df)
    print()
    print("TOTAL", df.sum())
