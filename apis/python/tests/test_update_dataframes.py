import tempfile
from pathlib import Path

import anndata
import numpy as np
import pandas as pd
import pyarrow as pa
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


def test_no_change(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o1 = exp.obs.schema
        v1 = exp.ms["RNA"].var.schema

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_obs(exp, adata.obs)
        tiledbsoma.io.update_var(exp, adata.var, "RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    assert o1 == o2
    assert v1 == v2


def test_add(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        exp.ms["RNA"].var.schema

    new_obs = adata.obs
    new_var = adata.var

    new_obs["is_g1"] = new_obs["groups"] == "g1"
    new_obs["seq"] = np.arange(new_obs.shape[0], dtype=np.int32)

    new_var["vst.mean.sq"] = new_var["vst.mean"] ** 2

    with tiledbsoma.Experiment.open(output_path, "w") as exp:
        tiledbsoma.io.update_obs(exp, new_obs)
        tiledbsoma.io.update_var(exp, new_var, "RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    assert o2.field("is_g1").type == pa.bool_()
    assert o2.field("seq").type == pa.int32()
    assert v2.field("vst.mean.sq").type == pa.float64()


def test_drop(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        exp.ms["RNA"].var.schema

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


def test_change(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o1 = exp.obs.schema
        v1 = exp.ms["RNA"].var.schema

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


@pytest.mark.parametrize("shift_and_exc", [[0, None], [1, ValueError]])
def test_change_counts(adata, shift_and_exc):
    shift, exc = shift_and_exc
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    tiledbsoma.io.from_anndata(output_path, adata, measurement_name="RNA")

    with tiledbsoma.Experiment.open(output_path) as exp:
        o1 = exp.obs.schema
        v1 = exp.ms["RNA"].var.schema

    old_nobs = len(adata.obs)
    old_nvar = len(adata.var)

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
