import tempfile
from dataclasses import dataclass, fields

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from anndata import AnnData
from pyarrow import Schema

import tiledbsoma
import tiledbsoma.io
from tiledbsoma._util import anndata_dataframe_unmodified


@dataclass
class TestCase:
    exp_path: str
    original: AnnData
    new: AnnData
    new_obs: pd.DataFrame
    new_var: pd.DataFrame
    o1: Schema
    v1: Schema


@pytest.fixture
def case(request, adata) -> TestCase:
    with tempfile.TemporaryDirectory() as exp_path:
        original = adata.copy()
        tiledbsoma.io.from_anndata(exp_path, adata, measurement_name="RNA")
        assert anndata_dataframe_unmodified(original.obs, adata.obs)
        assert anndata_dataframe_unmodified(original.var, adata.var)
        readback = request.param
        with tiledbsoma.Experiment.open(exp_path) as exp:
            o1 = exp.obs.schema
            v1 = exp.ms["RNA"].var.schema
            if readback:
                new_obs = exp.obs.read().concat().to_pandas()
                new_var = exp.ms["RNA"].var.read().concat().to_pandas()
            else:
                new_obs = adata.obs
                new_var = adata.var

            yield TestCase(
                exp_path=exp_path,
                original=original,
                new=adata,
                new_obs=new_obs,
                new_var=new_var,
                o1=o1,
                v1=v1,
            )


# Dynamically create a fixture for each field in the TestCase dataclass
def create_member_fixture(name):
    @pytest.fixture
    def member_fixture(case):
        return getattr(case, name)

    return member_fixture


# Register the dynamically created fixtures with pytest
for field in fields(TestCase):
    globals()[field.name] = create_member_fixture(field.name)


@pytest.mark.parametrize("case", [False, True], indirect=True)
def test_no_change(exp_path, original, new, new_obs, new_var, o1, v1):
    with tiledbsoma.Experiment.open(exp_path, "w") as exp:
        tiledbsoma.io.update_obs(exp, new_obs)
        tiledbsoma.io.update_var(exp, new_var, "RNA")
    assert anndata_dataframe_unmodified(original.obs, new.obs)
    assert anndata_dataframe_unmodified(original.var, new.var)

    with tiledbsoma.Experiment.open(exp_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    assert o1 == o2
    assert v1 == v2


@pytest.mark.parametrize("case", [False, True], indirect=True)
def test_add(exp_path, new_obs, new_var):
    # boolean
    new_obs["is_g1"] = new_obs["groups"] == "g1"
    # int
    new_obs["seq"] = np.arange(new_obs.shape[0], dtype=np.int32)
    # categorical of string
    new_obs["parity"] = pd.Categorical(
        np.asarray([["even", "odd"][e % 2] for e in range(len(new_obs))])
    )

    new_var["vst.mean.sq"] = new_var["vst.mean"] ** 2

    new_obs_save = new_obs.copy()
    new_var_save = new_var.copy()
    with tiledbsoma.Experiment.open(exp_path, "w") as exp:
        tiledbsoma.io.update_obs(exp, new_obs)
        tiledbsoma.io.update_var(exp, new_var, "RNA")
    assert anndata_dataframe_unmodified(new_obs, new_obs_save)
    assert anndata_dataframe_unmodified(new_var, new_var_save)

    with tiledbsoma.Experiment.open(exp_path) as exp:
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


@pytest.mark.parametrize("case", [False, True], indirect=True)
def test_drop(exp_path, new_obs, new_var):
    del new_obs["groups"]
    del new_var["vst.mean"]

    new_obs_save = new_obs.copy()
    new_var_save = new_var.copy()
    with tiledbsoma.Experiment.open(exp_path, "w") as exp:
        tiledbsoma.io.update_obs(exp, new_obs)
        tiledbsoma.io.update_var(exp, new_var, "RNA")
    assert anndata_dataframe_unmodified(new_obs, new_obs_save)
    assert anndata_dataframe_unmodified(new_var, new_var_save)

    with tiledbsoma.Experiment.open(exp_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    with pytest.raises(KeyError):
        o2.field("groups")
    with pytest.raises(KeyError):
        v2.field("vst.mean")


@pytest.mark.parametrize("case", [False, True], indirect=True)
def test_change(exp_path, new_obs, new_var, o1, v1):
    new_obs["groups"] = np.arange(new_obs.shape[0], dtype=np.int16)
    new_var["vst.mean"] = np.arange(new_var.shape[0], dtype=np.int32)

    new_obs_save = new_obs.copy()
    new_var_save = new_var.copy()
    with tiledbsoma.Experiment.open(exp_path, "w") as exp:
        with pytest.raises(ValueError):
            tiledbsoma.io.update_obs(exp, new_obs)
        with pytest.raises(ValueError):
            tiledbsoma.io.update_var(exp, new_var, "RNA")
    assert anndata_dataframe_unmodified(new_obs, new_obs_save)
    assert anndata_dataframe_unmodified(new_var, new_var_save)

    with tiledbsoma.Experiment.open(exp_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema

    assert o1 == o2
    assert v1 == v2


@pytest.mark.parametrize("case", [False, True], indirect=True)
@pytest.mark.parametrize("shift_and_exc", [[0, None], [1, ValueError]])
def test_change_counts(
    exp_path, original, new, new_obs, new_var, shift_and_exc, o1, v1
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
        new_obs2_save = new_obs2.copy()
        new_var2_save = new_var2.copy()
        with tiledbsoma.Experiment.open(exp_path, "w") as exp:
            tiledbsoma.io.update_obs(exp, new_obs2)
            tiledbsoma.io.update_var(exp, new_var2, measurement_name="RNA")

        assert anndata_dataframe_unmodified(new_obs2, new_obs2_save)
        assert anndata_dataframe_unmodified(new_var2, new_var2_save)

    else:
        with tiledbsoma.Experiment.open(exp_path, "w") as exp:
            with pytest.raises(exc):
                tiledbsoma.io.update_obs(exp, new_obs2)
            with pytest.raises(exc):
                tiledbsoma.io.update_var(exp, new_var2, measurement_name="RNA")

        assert anndata_dataframe_unmodified(original.obs, new.obs)
        assert anndata_dataframe_unmodified(original.var, new.var)

        with tiledbsoma.Experiment.open(exp_path) as exp:
            o2 = exp.obs.schema
            v2 = exp.ms["RNA"].var.schema
            assert o1 == o2
            assert v1 == v2
