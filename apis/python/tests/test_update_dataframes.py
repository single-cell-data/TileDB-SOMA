import tempfile
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa

# ruff: noqa
import pyarrow_hotfix
import pytest

import tiledbsoma
import tiledbsoma.io

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    input_path = HERE.parent / "testdata/pbmc-small.h5ad"
    # input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)


@pytest.mark.parametrize("readback", [False, True])
def test_no_change(adata, readback):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o1 = exp.obs.schema
        v1 = exp.ms["RNA"].var.schema

        if readback:
            new_obs = exp.obs.read().concat().to_pandas()
            new_var = exp.ms["RNA"].var.read().concat().to_pandas()
        else:
            new_obs = adata.obs
            new_var = adata.var

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_obs(exp, new_obs)
        tiledbsoma.io.update_var(exp, new_var, "RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    assert o1 == o2
    assert v1 == v2


@pytest.mark.parametrize("readback", [False, True])
def test_add(adata, readback):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        exp.ms["RNA"].var.schema

        if readback:
            new_obs = exp.obs.read().concat().to_pandas()
            new_var = exp.ms["RNA"].var.read().concat().to_pandas()
        else:
            new_obs = adata.obs
            new_var = adata.var

    # boolean
    new_obs["is_g1"] = new_obs["groups"] == "g1"
    # int
    new_obs["seq"] = np.arange(new_obs.shape[0], dtype=np.int32)
    # categorical of string
    new_obs["parity"] = pd.Categorical(
        np.asarray([["even", "odd"][e % 2] for e in range(len(new_obs))])
    )

    new_var["vst.mean.sq"] = new_var["vst.mean"] ** 2

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_obs(exp, new_obs)
        tiledbsoma.io.update_var(exp, new_var, "RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema
        obs = exp.obs.read().concat().to_pandas()

    assert o2.field("is_g1").type == pa.bool_()
    assert o2.field("seq").type == pa.int32()
    assert o2.field("parity").type == pa.dictionary(
        index_type=pa.int8(), value_type=pa.string(), ordered=False
    )
    assert obs["parity"][0] == "even"
    assert obs["parity"][1] == "odd"
    assert v2.field("vst.mean.sq").type == pa.float64()


@pytest.mark.parametrize("readback", [False, True])
def test_drop(adata, readback):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        exp.ms["RNA"].var.schema

        if readback:
            new_obs = exp.obs.read().concat().to_pandas()
            new_var = exp.ms["RNA"].var.read().concat().to_pandas()
        else:
            new_obs = adata.obs
            new_var = adata.var

    del new_obs["groups"]
    del new_var["vst.mean"]

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_obs(exp, new_obs)
        tiledbsoma.io.update_var(exp, new_var, "RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    with pytest.raises(KeyError):
        o2.field("groups")
    with pytest.raises(KeyError):
        v2.field("vst.mean")


@pytest.mark.parametrize("readback", [False, True])
def test_change(adata, readback):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o1 = exp.obs.schema
        v1 = exp.ms["RNA"].var.schema

        if readback:
            new_obs = exp.obs.read().concat().to_pandas()
            new_var = exp.ms["RNA"].var.read().concat().to_pandas()
        else:
            new_obs = adata.obs
            new_var = adata.var

    new_obs["groups"] = np.arange(new_obs.shape[0], dtype=np.int16)
    new_var["vst.mean"] = np.arange(new_var.shape[0], dtype=np.int32)

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        with pytest.raises(ValueError):
            tiledbsoma.io.update_obs(exp, new_obs)
        with pytest.raises(ValueError):
            tiledbsoma.io.update_var(exp, new_var, "RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    assert o1 == o2
    assert v1 == v2


@pytest.mark.parametrize("readback", [False, True])
@pytest.mark.parametrize("shift_and_exc", [[0, None], [1, ValueError]])
def test_change_counts(adata, readback, shift_and_exc):
    shift, exc = shift_and_exc
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o1 = exp.obs.schema
        v1 = exp.ms["RNA"].var.schema

        if readback:
            old_obs = exp.obs.read().concat().to_pandas()
            old_var = exp.ms["RNA"].var.read().concat().to_pandas()
        else:
            old_obs = adata.obs
            old_var = adata.var

    old_nobs = len(old_obs)
    old_nvar = len(old_var)

    new_nobs = old_nobs + shift
    new_nvar = old_nvar + shift

    new_obs = pd.DataFrame(
        data={
            "somebool": np.asarray([True] * new_nobs),
        },
        index=np.arange(new_nobs).astype(str),
    )
    new_var = pd.DataFrame(
        data={
            "somebool": np.asarray([True] * new_nvar),
        },
        index=np.arange(new_nvar).astype(str),
    )

    if exc is None:
        with tiledbsoma.Experiment.open(output_path, "w") as exp:
            tiledbsoma.io.update_obs(exp, new_obs)
            tiledbsoma.io.update_var(exp, new_var, measurement_name="RNA")

    else:
        with tiledbsoma.Experiment.open(output_path, "w") as exp:
            with pytest.raises(exc):
                tiledbsoma.io.update_obs(exp, new_obs)
            with pytest.raises(exc):
                tiledbsoma.io.update_var(exp, new_var, measurement_name="RNA")
        with tiledbsoma.Experiment.open(output_path) as exp:
            o2 = exp.obs.schema
            v2 = exp.ms["RNA"].var.schema
            assert o1 == o2
            assert v1 == v2
