from __future__ import annotations

import os

import pytest

import tiledbsoma
import tiledbsoma.io

from .util import util_pbmc3k_unprocessed_versions

if os.getenv("TILEDB_REST_UNITTEST_TOKEN") is None:
    pytest.skip(
        reason="$TILEDB_REST_UNITTEST_TOKEN is not set", allow_module_level=True
    )


@pytest.mark.parametrize(
    "uri_and_info",
    util_pbmc3k_unprocessed_versions(),
)
def test_basic_readback(conftest_context, uri_and_info):
    uri, info = uri_and_info
    with tiledbsoma.Experiment.open(uri, context=conftest_context) as exp:

        md = dict(exp.metadata)
        assert md["dataset_type"] == "soma"
        assert md["soma_object_type"] == "SOMAExperiment"

        md = dict(exp.obs.metadata)
        assert md["soma_object_type"] == "SOMADataFrame"

        md = dict(exp.ms["RNA"].var.metadata)
        assert md["soma_object_type"] == "SOMADataFrame"

        md = dict(exp.ms["RNA"].X["data"].metadata)
        assert md["soma_object_type"] == "SOMASparseNDArray"

        obs_table = exp.obs.read().concat()
        assert len(obs_table) == 2700
        obs_df = obs_table.to_pandas()
        assert obs_df.shape == (2700, 6)

        var_table = exp.ms["RNA"].var.read().concat()
        assert len(var_table) == 13714
        var_df = var_table.to_pandas()
        assert var_df.shape == (13714, 2)

        X_coo = exp.ms["RNA"].X["data"].read().coos().concat()
        if info["shape"] == "old":
            assert X_coo.shape == (2147483646, 2147483646)
        else:
            assert X_coo.shape == (2700, 13714)

        # Implicitly checking for no throw
        adata = tiledbsoma.io.to_anndata(exp, "RNA")

        assert adata.obs.shape == (2700, 4)
        assert adata.var.shape == (13714, 0)
        assert adata.X.shape == (2700, 13714)


@pytest.mark.parametrize(
    "uri_and_info",
    util_pbmc3k_unprocessed_versions(),
)
def test_dataframe_queries(conftest_context, uri_and_info):
    uri, info = uri_and_info
    with tiledbsoma.Experiment.open(uri, context=conftest_context) as exp:

        qobs = (
            exp.obs.read(
                coords=[slice(0, 99)],
                value_filter="nFeature_RNA > 1000",
                column_names=["soma_joinid", "obs_id", "nFeature_RNA"],
            )
            .concat()
            .to_pandas()
        )
        assert qobs.shape == (22, 3)

        qvar = (
            exp.ms["RNA"]
            .var.read(
                value_filter="var_id in ['ANXA1', 'IFI44', 'IFI44L', 'OAS1']",
            )
            .concat()
            .to_pandas()
        )
        assert qvar.shape == (4, 2)


@pytest.mark.parametrize(
    "uri_and_info",
    util_pbmc3k_unprocessed_versions(),
)
def test_experiment_queries(conftest_context, uri_and_info):
    uri, info = uri_and_info
    with tiledbsoma.Experiment.open(uri, context=conftest_context) as exp:

        query = tiledbsoma.ExperimentAxisQuery(
            experiment=exp,
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(
                value_filter="nFeature_RNA > 1000",
            ),
            var_query=tiledbsoma.AxisQuery(
                value_filter="var_id in  ['ANXA1', 'IFI44', 'IFI44L', 'OAS1']",
            ),
        )

        assert (query.n_obs, query.n_vars) == (530, 4)
