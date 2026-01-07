from __future__ import annotations

import io
import os

import pytest

import tiledbsoma
import tiledbsoma.io

from .util import util_pbmc3k_unprocessed_versions

if os.getenv("TILEDB_REST_UNITTEST_TOKEN") is None:
    pytest.skip(reason="$TILEDB_REST_UNITTEST_TOKEN is not set", allow_module_level=True)


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
    uri, _ = uri_and_info
    with tiledbsoma.Experiment.open(uri, context=conftest_context) as exp:
        qobs = (
            exp.obs
            .read(
                coords=[slice(0, 99)],
                value_filter="nFeature_RNA > 1000",
                column_names=["soma_joinid", "obs_id", "nFeature_RNA"],
            )
            .concat()
            .to_pandas()
        )
        assert qobs.shape == (22, 3)

        qvar = (
            exp
            .ms["RNA"]
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
    uri, _ = uri_and_info
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


@pytest.mark.parametrize(
    "uri_and_info",
    util_pbmc3k_unprocessed_versions(),
)
def test_upgrade_experiment_shapes(conftest_context, uri_and_info):
    uri, info = uri_and_info

    handle = io.StringIO()
    upgradeable = tiledbsoma.io.upgrade_experiment_shapes(
        uri,
        check_only=True,
        context=conftest_context,
        output_handle=handle,
    )
    handle.seek(0)
    lines = handle.readlines()
    handle.close()
    body = "\n".join(lines)

    assert "Dry run" in body
    if info["shape"] == "old":
        assert upgradeable
    else:
        assert not upgradeable
        assert "dataframe already has its domain set" in body


@pytest.mark.parametrize(
    "uri_and_info",
    util_pbmc3k_unprocessed_versions(),
)
def test_resize_experiment_too_small(conftest_context, uri_and_info):
    uri, _ = uri_and_info

    handle = io.StringIO()
    ok = tiledbsoma.io.resize_experiment(
        uri,
        nobs=10,
        nvars={"RNA": 20},
        check_only=True,
        context=conftest_context,
        output_handle=handle,
    )

    handle.seek(0)
    lines = handle.readlines()
    handle.close()
    body = "\n".join(lines)

    assert "Dry run" in body
    assert not ok


@pytest.mark.parametrize(
    "uri_and_info",
    util_pbmc3k_unprocessed_versions(),
)
def test_resize_experiment_ok(conftest_context, uri_and_info):
    uri, info = uri_and_info

    handle = io.StringIO()
    ok = tiledbsoma.io.resize_experiment(
        uri,
        nobs=100_000,
        nvars={"RNA": 200_000},
        check_only=True,
        context=conftest_context,
        output_handle=handle,
    )

    handle.seek(0)
    lines = handle.readlines()
    handle.close()
    body = "\n".join(lines)

    if info["shape"] == "old":
        assert not ok
        assert "dataframe currently has no domain set" in body
    else:
        assert ok


@pytest.mark.parametrize(
    "uri_and_info",
    util_pbmc3k_unprocessed_versions(),
)
def test_get_experiment_shapes(conftest_context, uri_and_info):
    uri, info = uri_and_info

    dict_output = tiledbsoma.io.get_experiment_shapes(uri, context=conftest_context)

    # Mask out the URIs for which we needn't compare the exact UUID values.
    dict_output["obs"]["uri"] = "test"
    dict_output["ms"]["RNA"]["var"]["uri"] = "test"
    dict_output["ms"]["RNA"]["X"]["data"]["uri"] = "test"

    if info["shape"] == "old":
        expect = {
            "obs": {
                "uri": "test",
                "type": "DataFrame",
                "count": 2700,
                "non_empty_domain": ((0, 2699),),
                "domain": ((0, 2147483646),),
                "maxdomain": ((0, 2147483646),),
                "upgraded": False,
            },
            "ms": {
                "RNA": {
                    "var": {
                        "uri": "test",
                        "type": "DataFrame",
                        "count": 13714,
                        "non_empty_domain": ((0, 13713),),
                        "domain": ((0, 2147483646),),
                        "maxdomain": ((0, 2147483646),),
                        "upgraded": False,
                    },
                    "X": {
                        "data": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2699), (0, 13713)),
                            "shape": (2147483646, 2147483646),
                            "maxshape": (2147483646, 2147483646),
                            "upgraded": False,
                        },
                    },
                },
            },
        }
        assert dict_output == expect

    else:
        expect = {
            "obs": {
                "uri": "test",
                "type": "DataFrame",
                "count": 2700,
                "non_empty_domain": ((0, 2699),),
                "domain": ((0, 2699),),
                "maxdomain": ((0, 9223372036854773758),),
                "upgraded": True,
            },
            "ms": {
                "RNA": {
                    "var": {
                        "uri": "test",
                        "type": "DataFrame",
                        "count": 13714,
                        "non_empty_domain": ((0, 13713),),
                        "domain": ((0, 13713),),
                        "maxdomain": ((0, 9223372036854773758),),
                        "upgraded": True,
                    },
                    "X": {
                        "data": {
                            "uri": "test",
                            "type": "SparseNDArray",
                            "non_empty_domain": ((0, 2699), (0, 13713)),
                            "shape": (2700, 13714),
                            "maxshape": (9223372036854773759, 9223372036854773759),
                            "upgraded": True,
                        },
                    },
                },
            },
        }
        assert dict_output == expect
