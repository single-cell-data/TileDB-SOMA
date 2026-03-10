"""
Test soma.io
"""

from __future__ import annotations

import pathlib
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
import pandas.testing
import scipy.sparse as sp

import tiledbsoma as soma
import tiledbsoma.io

from ._util import carrara_cleanup_asset
from .conftest import BASE_URI


def array_eq(a1, a2) -> bool:
    if isinstance(a1, np.ndarray) and isinstance(a2, np.ndarray):
        return np.array_equal(a1, a2)

    if isinstance(a1, sp.spmatrix) and isinstance(a2, sp.spmatrix):
        return (a1.tocsr() != a2.tocsr()).nnz == 0

    print(f"Oops, unsupported types: {type(a1)}, {type(a2)}")

    return False


def test_soma_io_roundtrip(small_pbmc: ad.AnnData, carrara_group_path: str, carrara_context: soma.SOMAContext) -> None:
    soma.io.from_anndata(carrara_group_path, small_pbmc, measurement_name="RNA", context=carrara_context)

    with soma.open(carrara_group_path, context=carrara_context) as exp:
        assert exp.obs.count == len(small_pbmc.obs)
        assert "RNA" in exp.ms
        assert exp.ms["RNA"].var.count == len(small_pbmc.var)

        assert set(exp.ms["RNA"].X) == {"data"} | set(small_pbmc.layers)
        assert set(exp.ms["RNA"].obsm) == set(small_pbmc.obsm)
        assert set(exp.ms["RNA"].varm) == set(small_pbmc.varm)
        assert set(exp.ms["RNA"].obsp) == set(small_pbmc.obsp)

        assert exp.ms["raw"].var.count == len(small_pbmc.raw.var)

        # spot checks

        obs_df = exp.obs.read().concat().to_pandas()[["obs_id", *list(small_pbmc.obs.keys())]].set_index("obs_id")
        obs_df.index.name = small_pbmc.obs.index.name
        pd.testing.assert_frame_equal(obs_df, small_pbmc.obs)

        var_df = (
            exp.ms["RNA"].var.read().concat().to_pandas()[["var_id", *list(small_pbmc.var.keys())]].set_index("var_id")
        )
        var_df.index.name = small_pbmc.var.index.name
        pd.testing.assert_frame_equal(var_df, small_pbmc.var)

        assert exp.ms["RNA"].X["data"].shape == small_pbmc.X.shape
        assert exp.ms["raw"].X["data"].shape == small_pbmc.raw.X.shape

        for slot in ("obsm", "varm", "obsp", "varp"):
            if slot in exp.ms["RNA"]:
                for k in exp.ms["RNA"][slot]:
                    assert exp.ms["RNA"][slot][k].shape == getattr(small_pbmc, slot)[k].shape

        # --- outgest and verify round-trip

        adata = soma.io.to_anndata(exp, measurement_name="RNA")
        assert adata.shape == (exp.obs.count, exp.ms["RNA"].var.count)
        adata.obs.index.name = small_pbmc.obs.index.name
        adata.var.index.name = small_pbmc.var.index.name
        pd.testing.assert_frame_equal(adata.obs, small_pbmc.obs)
        pd.testing.assert_frame_equal(adata.var, small_pbmc.var)

        assert (adata.X != small_pbmc.X).nnz == 0
        for slot in ("obsm", "varm", "obsp", "varp"):
            for k in getattr(small_pbmc, slot):
                assert array_eq(getattr(adata, slot)[k], getattr(small_pbmc, slot)[k]), f"{slot}[{k}] not EQ"


def test_soma_io_from_h5ad(
    tmp_path: pathlib.Path, small_pbmc: ad.AnnData, carrara_group_path: str, carrara_context: soma.SOMAContext
) -> None:
    """Test ability to ingest an H5AD sourced from a Carrara asset URL."""
    import tiledb.client

    local_small_pbmc_path = tmp_path / "small_pbmc.h5ad"
    small_pbmc.write(local_small_pbmc_path)

    with carrara_cleanup_asset(f"{BASE_URI}/{uuid4()}.h5ad") as h5ad_uri, open(local_small_pbmc_path, "rb") as fs:
        # stage H5AD to Carrara as a file asset
        tiledb.client.files.upload_file(fs, h5ad_uri)
        # ingest from Carrara asset to Carrara SOMA experiment
        soma.io.from_h5ad(carrara_group_path, h5ad_uri, measurement_name="RNA", context=carrara_context)
        # verify
        with soma.open(carrara_group_path) as exp:
            assert exp.obs.count == small_pbmc.n_obs
            assert exp.ms["RNA"].var.count == small_pbmc.n_vars
