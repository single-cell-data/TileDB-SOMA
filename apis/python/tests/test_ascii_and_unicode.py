import anndata as ad
import tiledb
from tiledbsc import SOMA
import tiledbsc.io as io
import pandas as pd
import numpy as np
from scipy import sparse

from pathlib import Path

import pytest


def test_readback(tmp_path):
    """
    Validate correct encode/decode of non-ASCII attribute text.
    """
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
        "HLA-DPB1",
        "IGLL5",
        "PARVB",
        "PGRMC1",
        "RP11-290F20.3",
        "SDPR",
    ]

    n_obs = len(obs_ids)
    n_var = len(var_ids)

    cell_types = ["blööd" if obs_id[1] == "A" else "lung" for obs_id in obs_ids]
    feature_names = [
        "ENSG00000999999" if var_id[1] < "M" else "ENSG00000123456"
        for var_id in var_ids
    ]

    obs = pd.DataFrame(
        data={
            "obs_id": np.asarray(obs_ids),
            "cell_type": np.asarray(cell_types),
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

    X = np.eye(n_obs, n_var)

    adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    soma = SOMA(tmp_path.as_posix())

    io.from_anndata(soma, adata)

    assert (
        len(soma.obs.attribute_filter('cell_type=="lung"', ["cell_type"])["cell_type"])
        == 6
    )
    assert (
        len(soma.obs.attribute_filter('cell_type=="blööd"', ["cell_type"])["cell_type"])
        == 4
    )
