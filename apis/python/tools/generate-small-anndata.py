#!/usr/bin/env python

# This is a little script to demonstrate how to write arbitrary test data into an AnnData .h5ad file.

import anndata as ad
import numpy as np
import pandas as pd

import random

obs_ids = [
    "AAATTCGAATCACG",
    "AATGTTGACAGTCA",
    "AGAGATGATCTCGC",
    "CATGGCCTGTGCAT",
    "CCCAACTGCAATCG",
    "CTAACGGAACCGAT",
    "GAACCTGATGAACC",
    "GCGTAAACACGGTT",
    "TTACCATGAATCGC",
    "TTACGTACGTTCAG",
]

var_ids = [
    "AKR1C3",
    "CA2",
    "CD1C",
    "GNLY",
    "HLA-DPB1",
    "HLA-DQA1",
    "IGLL5",
    "MYL9",
    "PARVB",
    "PF4",
    "PGRMC1",
    "PPBP",
    "RP11-290F20.3",
    "S100A9",
    "SDPR",
    "TREML1",
]

n_obs = len(obs_ids)
n_var = len(var_ids)

cell_types = ["blööd" if obs_id[1] == "A" else "lung" for obs_id in obs_ids]
feature_names = [
    "ENSG00000999999" if var_id[1] < "M" else "ENSG00000123456" for var_id in var_ids
]

# AnnData requires string indices for obs/var
obs = pd.DataFrame(
    data={
        "obs_id": np.asarray(obs_ids),
        "cell_type": np.asarray(cell_types),
        "is_primary_data": np.asarray([True] * n_obs),
    },
    index=np.arange(n_obs).astype(str),
)
obs.set_index("obs_id", inplace=True)
var = pd.DataFrame(
    data={
        "var_id": np.asarray(var_ids),
        "feature_name": np.asarray(feature_names),
    },
    index=np.arange(n_var).astype(str),
)
var.set_index("var_id", inplace=True)

X = np.zeros((n_obs, n_var))
for i in range(n_obs):
    for j in range(n_var):
        if random.uniform(0, 1) < 0.3:
            # if i == j and i != 1:
            X[i, j] = i * j

ann = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

output_file_name = "fake-small.h5ad"
ann.write_h5ad(output_file_name)
print("Wrote", output_file_name)
