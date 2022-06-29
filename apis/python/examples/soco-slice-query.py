#!/usr/bin/env python

import os
import shutil
import sys
from typing import List, Optional

import tiledbsc


# ----------------------------------------------------------------
def soco_query_and_store(
    soco: tiledbsc.SOMACollection,
    output_h5ad_path: str,
    output_soma_path: str,
    obs_attrs: Optional[List[str]] = None,
    obs_query_string: str = None,
    var_attrs: Optional[List[str]] = None,
    var_query_string: str = None,
) -> None:

    result_soma_slice = soco.query(
        obs_attrs=obs_attrs,
        obs_query_string=obs_query_string,
        var_attrs=var_attrs,
        var_query_string=var_query_string,
    )

    if result_soma_slice is None:
        print("Empty slice")
        return

    a = result_soma_slice.to_anndata()
    a.write_h5ad(output_h5ad_path)
    print("Wrote", output_h5ad_path, a.X.shape)

    if os.path.exists(output_soma_path):
        shutil.rmtree(output_soma_path)
    soma = tiledbsc.SOMA.from_soma_slice(result_soma_slice, output_soma_path)
    print("Wrote", output_soma_path, soma.X.data.shape())


# ----------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Need option")
        sys.exit(1)

    if sys.argv[1] == "1":
        print()
        print("TWO-SIDED QUERY")
        soco_query_and_store(
            soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
            output_h5ad_path="mini-atlas-two-sided.h5ad",
            output_soma_path="mini-atlas-two-sided",
            obs_attrs=["cell_type"],
            obs_query_string='cell_type == "B cell"',
            var_attrs=["feature_name"],
            var_query_string='feature_name == "MT-CO3"',
        )

    if sys.argv[1] == "2":
        print()
        print("OBS-ONLY QUERY")
        soco_query_and_store(
            soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
            output_h5ad_path="mini-atlas-obs-sided.h5ad",
            output_soma_path="mini-atlas-obs-sided",
            obs_attrs=["cell_type"],
            obs_query_string='cell_type == "B cell"',
        )

    if sys.argv[1] == "3":
        print()
        print("VAR-ONLY QUERY")
        soco_query_and_store(
            soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
            output_h5ad_path="mini-atlas-var-sided.h5ad",
            output_soma_path="mini-atlas-var-sided",
            var_attrs=["feature_name"],
            var_query_string='feature_name == "MT-CO3"',
        )

    if sys.argv[1] == "4":
        print()
        print("OBS-ONLY QUERY")
        soco_query_and_store(
            soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
            output_h5ad_path="cell-ontology-236.h5ad",
            output_soma_path="cell-ontology-236",
            obs_attrs=["cell_type_ontology_term_id"],
            obs_query_string='cell_type_ontology_term_id == "CL:0000236"',
        )

    if sys.argv[1] == "5":
        print()
        print("OBS-ONLY QUERY")
        soco_query_and_store(
            soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
            output_h5ad_path="kidney.h5ad",
            output_soma_path="kidney",
            obs_attrs=["tissue"],
            obs_query_string='tissue == "kidney"',
        )

    if sys.argv[1] == "6":
        print()
        print("OBS-ONLY QUERY")
        soco_query_and_store(
            soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
            output_h5ad_path="platelet.h5ad",
            output_soma_path="platelet",
            obs_attrs=["cell_type"],
            obs_query_string='cell_type == "platelet"',
        )

    if sys.argv[1] == "7":
        print()
        print("OBS-ONLY QUERY")
        soco_query_and_store(
            soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
            output_h5ad_path="platelet.h5ad",
            output_soma_path="platelet",
            obs_attrs=["cell_type", "tissue"],
            obs_query_string='cell_type == "B cell" and tissue == "blood"',
        )

    if sys.argv[1] == "8":
        print()
        print("OBS-ONLY QUERY")
        soco_query_and_store(
            soco=tiledbsc.SOMACollection("/Users/johnkerl/mini-corpus/atlas"),
            output_h5ad_path="platelet.h5ad",
            output_soma_path="platelet",
            obs_attrs=["cell_type", "tissue"],
            obs_query_string='cell_type == "B cell" or cell_type == "T cell"',
        )
