from __future__ import annotations

import os

import pytest

import tiledbsoma
import tiledbsoma.io
import scanpy as sc

from .util import util_make_uri

if os.getenv("TILEDB_REST_UNITTEST_TOKEN") is None:
    pytest.skip(
        reason="$TILEDB_REST_UNITTEST_TOKEN is not set", allow_module_level=True
    )


def test_basic_append(conftest_context, conftest_namespace, conftest_default_s3_path):
    (creation_uri, readback_uri) = util_make_uri(
        "soma-prod-ephemeral-data",
        "ephemeral_basic_append",
        conftest_namespace,
        conftest_default_s3_path,
    )

    measurement_name = "RNA"

    adata1 = sc.datasets.pbmc3k()
    adata1 .obs["when"] = ["Monday"] * len(adata1 .obs)
    tiledbsoma.io.from_anndata(creation_uri, adata1 , measurement_name=measurement_name)

    with tiledbsoma.Experiment.open(readback_uri) as exp:
        assert exp.obs.count == 2700
        assert exp.ms["RNA"].var.count == 32738
        assert exp.ms["RNA"].X["data"].shape == (2700, 32738)

    adata2 = sc.datasets.pbmc3k()
    adata2.obs.index = [e.replace("-1", "-2") for e in adata1 .obs.index]
    adata2.obs["when"] = ["Tuesday"] * len(adata2.obs)
    adata2.X *= 10

    rd = tiledbsoma.io.register_anndatas(
        readback_uri,
        [adata2],
        measurement_name=measurement_name,
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    tiledbsoma.io.resize_experiment(
        creation_uri,
        nobs=rd.get_obs_shape(),
        nvars=rd.get_var_shapes(),
    )

    tiledbsoma.io.from_anndata(
        creation_uri,
        adata2,
        measurement_name=measurement_name,
        registration_mapping=rd,
    )

    with tiledbsoma.Experiment.open(readback_uri) as exp:
        assert exp.obs.count == 5400
        assert exp.ms["RNA"].var.count == 32738
        assert exp.ms["RNA"].X["data"].shape == (5400, 32738)
