#!/usr/bin/env python

"""
Given one of the outputs from cartorapher.py, splits it into a few little pieces for demo purposes.
"""

import tiledbsc
import tiledbsc.io

import numpy as np
import os, shutil

# ----------------------------------------------------------------
def write_subset(soma, obs_indices, output_path):
    if os.path.exists(output_path):
        shutil.rmtree(output_path)
    subset_soma = tiledbsc.SOMA(output_path)
    subset_soma._create()

    obs_ids = list(soma.obs.df().index[obs_indices])
    var_ids = list(soma.var.df().index)

    subset_obs = soma.obs.df(obs_ids)
    subset_var = soma.var.df()
    subset_X_data = soma.X.data.csr(obs_ids, None)

    subset_obs["is_primary_data"] = np.asarray([True] * len(obs_ids))

    subset_soma.obs.from_dataframe(subset_obs, extent=2048)
    subset_soma.var.from_dataframe(subset_var, extent=2048)

    print("S OBS", subset_soma.obs.shape())
    print("S VAR", subset_soma.var.shape())
    print("S XDA", subset_X_data.shape)

    subset_soma.X.add_layer_from_matrix_and_dim_values(subset_X_data, obs_ids, var_ids)

    tiledbsc.io.to_h5ad(subset_soma, output_path + ".h5ad")


# ----------------------------------------------------------------
input_soma = tiledbsc.SOMA("atlas/4056cbab-2a32-4c9e-a55f-c930bc793fb6")

write_subset(input_soma, range(0, 100), "subset-soma-01")
write_subset(input_soma, range(100, 200), "subset-soma-02")
write_subset(input_soma, range(200, 300), "subset-soma-03")
write_subset(input_soma, range(300, 400), "subset-soma-04")
