#!/usr/bin/env python

import tiledbsc as t
import pandas as pd
import sys
from typing import List, Dict

# ================================================================
def main():
    if len(sys.argv) == 2:
        soco_path = sys.argv[1]
    else:
        print(f"{sys.argv[0]}: need just one soma-collection path.", file=sys.stderr)
        sys.exit(1)

    soco = tiledbsc.SOMACollection(soco_path)

    #    print()
    #    print("================================================================")
    #    print("NAMES AND URIS")
    #    print_names_and_uris(soco)
    #
    #    print()
    #    print("================================================================")
    #    print("OBS NAMES")
    #    show_obs_names(soco)
    #
    #    print()
    #    print("================================================================")
    #    print("VAR NAMES")
    #    show_var_names(soco)
    #
    #    print()
    #    print("================================================================")
    #    print("SOMAS HAVING ALL THREE")
    #    show_somas_with_all_three(soco)
    #
    #    print()
    #    print("================================================================")
    #    print("OBS_ID COUNTS")
    #    show_obs_id_counts(soco)
    #
    #    print()
    #    print("================================================================")
    #    print("VAR_ID COUNTS")
    #    show_var_id_counts(soco)
    #
    #    print()
    #    print("================================================================")
    #    print("OBS UNIQUE VALUES FOR CELL_TYPE")
    #    show_obs_column_unique_values(soco, "cell_type")
    #
    #    print()
    #    print("================================================================")
    #    print("OBS UNIQUE VALUES FOR FEATURE_NAME")
    #    show_var_column_unique_values(soco, "feature_name")
    #
    #    print()
    #    print("================================================================")
    #    print("OBS VALUE COUNTS FOR CELL_TYPE AND TISSUE")
    #    show_obs_value_counts(soco, ["cell_type", "tissue"])
    #
    #    print()
    #    print("================================================================")
    #    print("VAR VALUE COUNTS FOR CELL_TYPE AND FEATURE_NAME")
    #    show_var_value_counts(soco, ["feature_name"])

    print()
    print("================================================================")
    print("SHOW SOMAS HAVING")
    show_somas_having(
        soco,
        {"cell_type": ["B cell", "T cell"], "tissue": ["blood", "lung"]},
        {"feature_name": ["MT-CO3"]},
    )


# ----------------------------------------------------------------
def print_names_and_uris(soco: tiledbsc.SOMACollection) -> None:
    for soma in soco:
        print("%-40s %s" % (soma.name, soma.uri))


# ----------------------------------------------------------------
def show_obs_names(soco: tiledbsc.SOMACollection) -> None:
    for soma in soco:
        print(soma.uri)
        for attr_name in soma.obs.keys():
            print("  obs", attr_name)


# ----------------------------------------------------------------
def show_var_names(soco: tiledbsc.SOMACollection) -> None:
    for soma in soco:
        print(soma.uri)
        for attr_name in soma.var.keys():
            print("  var", attr_name)


# ----------------------------------------------------------------
def show_somas_with_all_three(soco: tiledbsc.SOMACollection) -> None:
    for soma in soco:
        if "cell_type" in soma.obs.attr_names():
            if "tissue" in soma.obs.attr_names():
                if "feature_name" in soma.var.attr_names():
                    print(soma.uri)


# ----------------------------------------------------------------
def show_obs_id_counts(soco: tiledbsc.SOMACollection) -> None:
    counts = {}
    for soma in soco:
        for oid in soma.obs.ids():
            if oid in counts:
                counts[oid] += 1
            else:
                counts[oid] = 1
    df = pd.DataFrame.from_dict({"obs_id": counts.keys(), "counts": counts.values()})
    # print(df.head())
    print(df)


# ----------------------------------------------------------------
def show_var_id_counts(soco: tiledbsc.SOMACollection) -> None:
    counts = {}
    for soma in soco:
        for oid in soma.var.ids():
            if oid in counts:
                counts[oid] += 1
            else:
                counts[oid] = 1
    df = pd.DataFrame.from_dict({"var_id": counts.keys(), "counts": counts.values()})
    # print(df.head())
    print(df)


# ----------------------------------------------------------------
def show_obs_column_unique_values(soco: tiledbsc.SOMACollection, col_name: str) -> None:
    for soma in soco:
        print()
        print(soma.uri)
        if col_name in soma.obs.keys():
            print(soma.obs.df()[col_name].unique())


# ----------------------------------------------------------------
def show_var_column_unique_values(soco: tiledbsc.SOMACollection, col_name: str) -> None:
    for soma in soco:
        print()
        print(soma.uri)
        if col_name in soma.var.keys():
            print(soma.var.df()[col_name].unique())


# ----------------------------------------------------------------
def show_obs_value_counts(soco: tiledbsc.SOMACollection, obs_labels: List[str]) -> None:

    for obs_label in obs_labels:
        counts = {}

        for soma in soco:
            print("...", soma.name)
            # print("\n".join(sorted(soma.obs.attr_names())))
            obs = soma.obs.df()

            if not obs_label in obs:
                continue

            obs_label_values = sorted(list(set(obs[obs_label])))
            for obs_label_value in obs_label_values:
                if obs_label_value in counts:
                    counts[obs_label_value] += 1
                else:
                    counts[obs_label_value] = 1

        print(
            "----------------------------------------------------------------",
            obs_label,
        )
        for k, v in dict(sorted(counts.items(), key=lambda item: item[1])).items():
            print(k, v)


# ----------------------------------------------------------------
def show_var_value_counts(soco: tiledbsc.SOMACollection, var_labels: List[str]) -> None:

    for var_label in var_labels:
        counts = {}

        for soma in soco:
            print("...", soma.name)
            # print("\n".join(sorted(soma.var.attr_names())))
            var = soma.var.df()

            if not var_label in var:
                continue

            var_label_values = sorted(list(set(var[var_label])))
            for var_label_value in var_label_values:
                if var_label_value in counts:
                    counts[var_label_value] += 1
                else:
                    counts[var_label_value] = 1

        print(
            "----------------------------------------------------------------",
            var_label,
        )
        for k, v in dict(sorted(counts.items(), key=lambda item: item[1])).items():
            print(k, v)


# ----------------------------------------------------------------
def show_somas_having(
    soco: tiledbsc.SOMACollection,
    obs_labels_to_values: Dict[str, List],
    var_labels_to_values: Dict[str, List],
) -> None:

    for soma in soco:
        print(soma.uri)

        obs = soma.obs.df()
        for obs_label in obs_labels_to_values:
            if not obs_label in obs:
                print("out1")
                continue
            soma_obs_label_values = sorted(list(set(obs[obs_label])))
            for sought_obs_label_value in obs_labels_to_values[obs_label]:
                if sought_obs_label_value in soma_obs_label_values:
                    print("  found obs", sought_obs_label_value)

        var = soma.var.df()
        for var_label in var_labels_to_values:
            if not var_label in var:
                print("out2")
                continue
            soma_var_label_values = sorted(list(set(var[var_label])))
            for sought_var_label_value in var_labels_to_values[var_label]:
                if sought_var_label_value in soma_var_label_values:
                    print("  found var", sought_var_label_value)


# ================================================================
if __name__ == "__main__":
    main()
