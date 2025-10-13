"""
Test soma.io
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pandas.testing

import tiledbsoma as soma
import tiledbsoma.io


def test_soma_io_import(small_pbmc: ad.AnnData, group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    soma.io.from_anndata(group_path, small_pbmc, measurement_name="RNA", context=carrara_context)

    with soma.open(group_path, context=carrara_context) as exp:
        print(small_pbmc)
        print(exp)
        print(exp.ms["RNA"])
        for k in exp.ms["RNA"]:
            print(exp.ms["RNA"][k])
        for k in exp.ms["raw"]:
            print(exp.ms["raw"][k])
        print(exp)

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
        pd.testing.assert_frame_equal(adata.raw.var, small_pbmc.raw.var)
        assert (adata.X != small_pbmc.X).nnz == 0
        assert (adata.raw.X != small_pbmc.raw.X).nnz == 0
        for slot in ("obsm", "varm", "obsp", "varp"):
            for k in getattr(small_pbmc, slot):
                if isinstance(getattr(adata, slot)[k], np.ndarray):
                    assert getattr(adata, slot)[k] == getattr(small_pbmc, slot)[k], f"{slot}[{k}] not EQ"
                else:
                    assert (getattr(adata, slot)[k] != getattr(small_pbmc, slot)[k]).nnz == 0, f"{slot}[{k}] not EQ"
