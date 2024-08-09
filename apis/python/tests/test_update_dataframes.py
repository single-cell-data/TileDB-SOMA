import tempfile
from dataclasses import dataclass
from typing import Optional, Type

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from anndata import AnnData
from pandas._testing import assert_frame_equal
from pyarrow import Schema

import tiledbsoma
import tiledbsoma.io

from ._util import (
    assert_adata_equal,
    maybe_raises,
)


@dataclass
class TestCase:
    """Group of related objects, used by test cases in this file, that share a creation code-path and are exposed as
    `pytest.fixture`s below."""

    experiment_path: str
    old_anndata: AnnData
    new_anndata: AnnData
    new_obs: pd.DataFrame
    new_var: pd.DataFrame
    obs_schema: Schema
    var_schema: Schema


@pytest.fixture
def multiple_fixtures_with_readback(request, conftest_pbmc_small) -> TestCase:
    """
    Ingest an `AnnData`to a SOMA `Experiment`, yielding a `TestCase` with the old and new AnnData objects.

    * The input AnnData is always `conftest_pbmc_small` (from conftest.py), not specifiable here as an argument
    * `request.param` (a.k.a. `use_readback` below) is a boolean, populated by pytest, that specifies whether:
      * `False`: the returned `new_obs` and `new_var` come directly from the input `conftest_pbmc_small`, or
      * `True`: `new_obs` and `new_var` are returned after round-tripping `conftest_pbmc_small` through SOMA and back to
        AnnData/Pandas.
    * Each `TestCase` member is also exposed directly as its own `fixture` below.
    """
    with tempfile.TemporaryDirectory() as experiment_path:
        old_anndata = conftest_pbmc_small.copy()
        tiledbsoma.io.from_anndata(
            experiment_path, conftest_pbmc_small, measurement_name="RNA"
        )

        # Check that the anndata-to-soma ingestion didn't modify the old_anndata object (which is
        # passed by reference to the ingestor) while it was doing the ingest
        assert_adata_equal(old_anndata, conftest_pbmc_small)

        use_readback = request.param

        with tiledbsoma.Experiment.open(experiment_path) as exp:
            obs_schema = exp.obs.schema
            var_schema = exp.ms["RNA"].var.schema
            if use_readback:
                new_obs = exp.obs.read().concat().to_pandas()
                new_var = exp.ms["RNA"].var.read().concat().to_pandas()
            else:
                new_obs = conftest_pbmc_small.obs
                new_var = conftest_pbmc_small.var

            yield TestCase(
                experiment_path=experiment_path,
                old_anndata=old_anndata,
                new_anndata=conftest_pbmc_small,
                new_obs=new_obs,
                new_var=new_var,
                obs_schema=obs_schema,
                var_schema=var_schema,
            )


# Expose each field of the TestCase object as a pytest.fixture, to reduce boilerplate in test functions.
@pytest.fixture
def experiment_path(multiple_fixtures_with_readback):
    return multiple_fixtures_with_readback.experiment_path


@pytest.fixture
def old_anndata(multiple_fixtures_with_readback):
    return multiple_fixtures_with_readback.old_anndata


@pytest.fixture
def new_anndata(multiple_fixtures_with_readback):
    return multiple_fixtures_with_readback.new_anndata


@pytest.fixture
def new_obs(multiple_fixtures_with_readback):
    return multiple_fixtures_with_readback.new_obs


@pytest.fixture
def new_var(multiple_fixtures_with_readback):
    return multiple_fixtures_with_readback.new_var


@pytest.fixture
def obs_schema(multiple_fixtures_with_readback):
    return multiple_fixtures_with_readback.obs_schema


@pytest.fixture
def var_schema(multiple_fixtures_with_readback):
    return multiple_fixtures_with_readback.var_schema


def verify_schemas(experiment_path, obs_schema, var_schema):
    """Read {obs,var} schemas, verify they match initial versions."""
    with tiledbsoma.Experiment.open(experiment_path) as exp:
        other_obs_schema = exp.obs.schema
        other_var_schema = exp.ms["RNA"].var.schema
    assert obs_schema == other_obs_schema
    assert var_schema == other_var_schema


def verify_updates(
    experiment_path,
    obs,
    var,
    exc: Optional[Type[ValueError]] = None,
):
    """
    Calls `update_obs` and `update_var` on the experiment. Also verifies that the
    updater code didn't inadvertently modify the `obs` and `var` objects, which are
    passed into the updater by reference.
    """

    obs0 = obs.copy()
    var0 = var.copy()

    with tiledbsoma.Experiment.open(experiment_path, "w") as exp:
        with maybe_raises(exc):
            tiledbsoma.io.update_obs(exp, obs)
        with maybe_raises(exc):
            tiledbsoma.io.update_var(exp, var, "RNA")

    assert_frame_equal(obs0, obs)
    assert_frame_equal(var0, var)


# `pytest.mark.parametrize` wrapper for running a test twice:
# 1. `readback=False`: `new_obs` and `new_var` come directly from the input `conftest_pbmc_small` object
# 2. `readback=True`: `new_obs` and `new_var` are ingested to SOMA, then exported back to pandas DataFrames.
with_and_without_soma_roundtrip = pytest.mark.parametrize(
    "multiple_fixtures_with_readback", [False, True], indirect=True
)


@with_and_without_soma_roundtrip
def test_no_change(
    experiment_path, old_anndata, new_anndata, new_obs, new_var, obs_schema, var_schema
):
    verify_updates(experiment_path, new_obs, new_var)
    verify_schemas(experiment_path, obs_schema, var_schema)
    assert_adata_equal(old_anndata, new_anndata)


@with_and_without_soma_roundtrip
def test_add(experiment_path, new_obs, new_var):
    # boolean
    new_obs["is_g1"] = new_obs["groups"] == "g1"
    # int
    new_obs["seq"] = np.arange(new_obs.shape[0], dtype=np.int32)
    # categorical of string
    new_obs["parity"] = pd.Categorical(
        np.asarray([["even", "odd"][e % 2] for e in range(len(new_obs))])
    )

    new_var["vst.mean.sq"] = new_var["vst.mean"] ** 2

    verify_updates(experiment_path, new_obs, new_var)

    with tiledbsoma.Experiment.open(experiment_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema
        obs = exp.obs.read().concat().to_pandas()

    assert o2.field("is_g1").type == pa.bool_()
    assert o2.field("seq").type == pa.int32()
    # tiledbsoma.io upgrades int8 and int16 to int32 for appendability
    assert o2.field("parity").type == pa.dictionary(
        index_type=pa.int32(), value_type=pa.string(), ordered=False
    )
    assert obs["parity"][0] == "even"
    assert obs["parity"][1] == "odd"
    assert v2.field("vst.mean.sq").type == pa.float64()


@with_and_without_soma_roundtrip
def test_drop(experiment_path, new_obs, new_var):
    del new_obs["groups"]
    del new_var["vst.mean"]

    verify_updates(experiment_path, new_obs, new_var)

    with tiledbsoma.Experiment.open(experiment_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    with pytest.raises(KeyError):
        o2.field("groups")
    with pytest.raises(KeyError):
        v2.field("vst.mean")


@with_and_without_soma_roundtrip
def test_change(experiment_path, new_obs, new_var, obs_schema, var_schema):
    new_obs["groups"] = np.arange(new_obs.shape[0], dtype=np.int16)
    new_var["vst.mean"] = np.arange(new_var.shape[0], dtype=np.int32)
    verify_updates(experiment_path, new_obs, new_var, exc=ValueError)
    verify_schemas(experiment_path, obs_schema, var_schema)


@with_and_without_soma_roundtrip
@pytest.mark.parametrize("shift_and_exc", [[0, None], [1, ValueError]])
def test_change_counts(
    experiment_path,
    old_anndata,
    new_anndata,
    new_obs,
    new_var,
    obs_schema,
    var_schema,
    shift_and_exc,
):
    shift, exc = shift_and_exc

    new_nobs = len(new_obs)
    new_nvar = len(new_var)

    new_nobs2 = new_nobs + shift
    new_nvar2 = new_nvar + shift

    new_obs2 = pd.DataFrame(
        data={
            "somebool": np.asarray([True] * new_nobs2),
        },
        index=np.arange(new_nobs2).astype(str),
    )
    new_var2 = pd.DataFrame(
        data={
            "somebool": np.asarray([True] * new_nvar2),
        },
        index=np.arange(new_nvar2).astype(str),
    )

    if exc is None:
        verify_updates(experiment_path, new_obs2, new_var2)
    else:
        verify_updates(experiment_path, new_obs2, new_var2, exc=ValueError)
        assert_adata_equal(old_anndata, new_anndata)
        verify_schemas(experiment_path, obs_schema, var_schema)


@pytest.mark.parametrize("separate_ingest", [False, True])
def test_update_non_null_to_null(tmp_path, conftest_pbmc3k_adata, separate_ingest):
    uri = tmp_path.as_uri()

    # Two ways to test:
    #
    # One way:
    # * Create columns with non-nulls before from_anndata
    # * Ingest, having those new columns with non-nulls
    # * Call update_obs to set the columns to have nulls
    #
    # Other way:
    # * Ingest, without any new columns
    # * Call update_obs once to add the new columns with non-null values
    # * Call update_obs again to add the new columns with null values

    if separate_ingest:
        tiledbsoma.io.from_anndata(
            uri,
            conftest_pbmc3k_adata,
            measurement_name="RNA",
            uns_keys=[],
        )

        conftest_pbmc3k_adata.obs["batch_id"] = "testing"
        conftest_pbmc3k_adata.obs["myfloat"] = 12.34
        verify_updates(uri, conftest_pbmc3k_adata.obs, conftest_pbmc3k_adata.var)

    else:
        conftest_pbmc3k_adata.obs["batch_id"] = "testing"
        conftest_pbmc3k_adata.obs["myfloat"] = 12.34

        tiledbsoma.io.from_anndata(
            uri,
            conftest_pbmc3k_adata,
            measurement_name="RNA",
            uns_keys=[],
        )

    conftest_pbmc3k_adata.obs["batch_id"] = pd.NA
    conftest_pbmc3k_adata.obs["myfloat"] = np.nan
    # We need nan_safe since pd.NA != pd.NA
    verify_updates(uri, conftest_pbmc3k_adata.obs, conftest_pbmc3k_adata.var)
