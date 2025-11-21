# These are test that need to run first to check basic functionality, before we
# go on to test other, more complex things.

from __future__ import annotations

import os
import sys

import pytest
import scanpy as sc

import tiledbsoma
import tiledbsoma.io
import tiledb.cloud

from .util import util_make_uri, util_tear_down_uri

# Nominally this is the 'unittest' SaaS user. What we require is:
#
# * The user can _read_ data in the 'unittest' namespace.
# * For data _written_, the namespace and default_s3_path are taken from the
#   cloud profile.
#
# For CI, this environment variable is a GitHub Actions secret, propagated in
# the CI YAML.
if os.getenv("TILEDB_REST_UNITTEST_TOKEN") is None:
    pytest.skip(reason="$TILEDB_REST_UNITTEST_TOKEN is not set", allow_module_level=True)


def test_skipping_correctly():
    assert os.getenv("TILEDB_REST_UNITTEST_TOKEN") is not None


def test_basic_read(conftest_context):
    uri = "tiledb://unittest/pbmc3k_unprocessed_1_15_7"
    assert tiledbsoma.Experiment.exists(uri, context=conftest_context)
    with tiledbsoma.Experiment.open(uri, context=conftest_context) as exp:
        assert exp.obs.count == 2700
        assert "RNA" in exp.ms
        assert exp.ms["RNA"].var.count == 13714


def test_basic_write(conftest_context, conftest_namespace, conftest_default_s3_path):
    (creation_uri, readback_uri) = util_make_uri(
        "soma-prod-ephemeral-data",
        "ephemeral_basic_write",
        conftest_namespace,
        conftest_default_s3_path,
    )

    adata = sc.datasets.pbmc3k()

    tiledbsoma.io.from_anndata(
        creation_uri,
        adata,
        measurement_name="RNA",
        context=conftest_context,
    )

    with tiledbsoma.Experiment.open(readback_uri, context=conftest_context) as exp:
        assert exp.obs.count == 2700
        assert "RNA" in exp.ms
        assert exp.ms["RNA"].var.count == 32738

    util_tear_down_uri(readback_uri)


@pytest.mark.skipif(
    (sys.version_info.major, sys.version_info.minor) != (3, 9),
    reason="As of 2025-02-05 UDFs require Python 3.9",
)
def test_remote_version(conftest_tiledb_cloud_login):
    def remote_version():
        import tiledbsoma

        return {"tiledbsoma": tiledbsoma.__version__}

    output = tiledb.cloud.udf.exec(remote_version)
    assert "tiledbsoma" in output
    assert output["tiledbsoma"].startswith("2.")
