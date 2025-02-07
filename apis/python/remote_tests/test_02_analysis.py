from __future__ import annotations

import os

import pandas as pd
import pytest
import scanpy as sc

import tiledbsoma
import tiledbsoma.io
import tiledbsoma.logging

from .util import util_make_uri, util_tear_down_uri

if os.getenv("TILEDB_REST_UNITTEST_TOKEN") is None:
    pytest.skip(
        reason="$TILEDB_REST_UNITTEST_TOKEN is not set", allow_module_level=True
    )


def test_write_with_updates(
    conftest_context, conftest_namespace, conftest_default_s3_path
):
    (creation_uri, readback_uri) = util_make_uri(
        "soma-prod-ephemeral-data",
        "ephemeral_analysis",
        conftest_namespace,
        conftest_default_s3_path,
    )

    adata = sc.datasets.pbmc3k()

    tiledbsoma.logging.info()
    tiledbsoma.io.from_anndata(
        creation_uri,
        adata,
        measurement_name="RNA",
        context=conftest_context,
    )

    with tiledbsoma.Experiment.open(readback_uri, context=conftest_context) as exp:
        assert "RNA" in exp.ms

        assert exp.metadata.get("dataset_type") == "soma"
        assert exp.metadata.get("soma_object_type") == "SOMAExperiment"
        assert exp.obs.metadata.get("soma_object_type") == "SOMADataFrame"
        assert exp.ms["RNA"].var.metadata.get("soma_object_type") == "SOMADataFrame"
        assert "data" in exp.ms["RNA"].X
        assert (
            exp.ms["RNA"].X["data"].metadata.get("soma_object_type")
            == "SOMASparseNDArray"
        )

        assert exp.obs.count == adata.obs.shape[0]
        assert exp.ms["RNA"].var.count == adata.var.shape[0]

        obs_arrow = exp.obs.read().concat()
        obs_pandas = obs_arrow.to_pandas()
        assert obs_pandas.shape[0] == adata.obs.shape[0]

    # Here we augment that with some on-the-fly computed data. This imitates a common customer workflow.
    # Add a categorical column
    parity = [["even", "odd"][e % 2] for e in range(len(adata.obs))]
    adata.obs["parity"] = pd.Categorical(parity)
    with tiledbsoma.Experiment.open(creation_uri, "w", context=conftest_context) as exp:
        tiledbsoma.io.update_obs(exp, adata.obs, context=conftest_context)

    with tiledbsoma.Experiment.open(readback_uri, context=conftest_context) as exp:
        obs_arrow = exp.obs.read().concat()
        obs_pandas = obs_arrow.to_pandas()
        assert obs_pandas.shape[0] == adata.obs.shape[0]

    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata, copy=False)

    with tiledbsoma.open(creation_uri, "w", context=conftest_context) as exp:
        tiledbsoma.io.add_X_layer(
            exp,
            measurement_name="RNA",
            X_layer_name="logcounts",
            X_layer_data=adata.X,
            context=conftest_context,
        )

    with tiledbsoma.open(readback_uri, "w", context=conftest_context) as exp:
        assert sorted(list(exp.ms["RNA"].X.keys())) == ["data", "logcounts"]

    # Add dimensional-reduction results
    sc.pp.highly_variable_genes(adata, inplace=True)
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata)
    sc.tl.pca(adata, use_highly_variable=True, n_comps=5)

    with tiledbsoma.open(creation_uri, "w", context=conftest_context) as exp:
        tiledbsoma.io.add_matrix_to_collection(
            exp=exp,
            measurement_name="RNA",
            collection_name="obsm",
            matrix_name="logcounts_pca",
            matrix_data=adata.obsm["X_pca"],
            context=conftest_context,
        )

    with tiledbsoma.open(readback_uri, "w", context=conftest_context) as exp:
        assert sorted(list(exp.ms["RNA"].obsm.keys())) == ["logcounts_pca"]

    with tiledbsoma.open(creation_uri, "w", context=conftest_context) as exp:
        tiledbsoma.io.add_matrix_to_collection(
            exp=exp,
            measurement_name="RNA",
            collection_name="varm",
            matrix_name="logcounts_pcs",
            matrix_data=adata.varm["PCs"],
            context=conftest_context,
        )
    with tiledbsoma.open(exp.uri, context=conftest_context) as exp:
        assert sorted(list(exp.ms["RNA"].varm.keys())) == ["logcounts_pcs"]

    util_tear_down_uri(readback_uri)
