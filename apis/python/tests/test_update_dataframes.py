import tempfile
from contextlib import nullcontext
from dataclasses import dataclass, fields

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from anndata import AnnData
from pyarrow import Schema

import tiledbsoma
import tiledbsoma.io
from tiledbsoma._util import anndata_dataframe_unmodified, verify_obs_var


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
    """Ingest an ``AnnData``, yield a ``TestCase`` with the original and new AnnData objects."""
    with tempfile.TemporaryDirectory() as exp_path:
        original = adata.copy()
        tiledbsoma.io.from_anndata(exp_path, adata, measurement_name="RNA")
        verify_obs_var(original, adata)
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


def create_member_fixture(name):
    """Create a ``pytest.fixture`` for a ``TestCase`` field."""

    @pytest.fixture
    def member_fixture(case):
        return getattr(case, name)

    return member_fixture


for field in fields(TestCase):
    """Create ``pytest.fixture``s for each ``TestCase`` field."""
    globals()[field.name] = create_member_fixture(field.name)


def verify_schemas(exp_path, o1, v1):
    """Read {obs,var} schemas, verify they match initial versions."""
    with tiledbsoma.Experiment.open(exp_path) as exp:
        o2 = exp.obs.schema
        v2 = exp.ms["RNA"].var.schema
    assert o1 == o2
    assert v1 == v2


def verify_updates(exp_path, obs, var, exc=False):
    obs0 = obs.copy()
    var0 = var.copy()

    def ctx():
        return pytest.raises(ValueError) if exc else nullcontext()

    with tiledbsoma.Experiment.open(exp_path, "w") as exp:
        with ctx():
            tiledbsoma.io.update_obs(exp, obs)
        with ctx():
            tiledbsoma.io.update_var(exp, var, "RNA")

    assert anndata_dataframe_unmodified(obs0, obs)
    assert anndata_dataframe_unmodified(var0, var)


@pytest.mark.parametrize("case", [False, True], indirect=True)
def test_no_change(exp_path, original, new, new_obs, new_var, o1, v1):
    verify_updates(exp_path, new_obs, new_var)
    verify_schemas(exp_path, o1, v1)
    verify_obs_var(original, new)


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

    verify_updates(exp_path, new_obs, new_var)

    with tiledbsoma.Experiment.open(exp_path) as exp:
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


@pytest.mark.parametrize("case", [False, True], indirect=True)
def test_drop(exp_path, new_obs, new_var):
    del new_obs["groups"]
    del new_var["vst.mean"]

    verify_updates(exp_path, new_obs, new_var)

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
    verify_updates(exp_path, new_obs, new_var, exc=True)
    verify_schemas(exp_path, o1, v1)


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
        verify_updates(exp_path, new_obs2, new_var2)
    else:
        verify_updates(exp_path, new_obs2, new_var2, exc=True)
        verify_obs_var(original, new)
        verify_schemas(exp_path, o1, v1)
