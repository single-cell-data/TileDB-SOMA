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
from tiledbsoma._util import anndata_dataframe_unmodified, verify_obs_and_var_same


@dataclass
class TestCase:
    """Holds pieces which would otherwise be multiple fixtures in pytest.mark.parametrize."""

    experiment_path: str
    old_anndata: AnnData
    new_anndata: AnnData
    new_obs: pd.DataFrame
    new_var: pd.DataFrame
    obs_schema: Schema
    var_schema: Schema


# Magical conftest.py fixture: conftest_adata
# Also magical: request
@pytest.fixture
def multiple_fixtures_with_readback(request, conftest_adata) -> TestCase:
    """
    Ingest an ``AnnData``to a SOMA ``Experiment``, yielding a ``TestCase`` with the old and new AnnData objects.

    * The input AnnData is from a hard-coded fixture in conftest.py, not specifiable here as an argument
    * The `request` is nominally for tests to give a boolean for (False) whether they want the
      `new_obs` and `new_var` to be the same as the input conftest_adata, or (True) whether they want
      the `new_obs` and `new_var` to be gotten from writing the conftest_adata object via
      tiledbsoma.io.from_anndata to a SOMA experiment and then read back via the SOMA API and
      to_pandas() calls.
    * The slots in `TestCase` will be available to tests which use this multi-fixture:
      e.g. if the TestCase class has a `foo` slot then a unit-test function in this file can be
      passed a `foo` argument
    """
    with tempfile.TemporaryDirectory() as experiment_path:
        old_anndata = conftest_adata.copy()
        tiledbsoma.io.from_anndata(
            experiment_path, conftest_adata, measurement_name="RNA"
        )

        # Check that the anndata-to-soma ingestion didn't modify the old_anndata object (which is
        # passed by reference to the ingestor) while it was doing the ingest
        verify_obs_and_var_same(old_anndata, conftest_adata)

        use_readback = request.param

        with tiledbsoma.Experiment.open(experiment_path) as exp:
            obs_schema = exp.obs.schema
            var_schema = exp.ms["RNA"].var.schema
            if use_readback:
                new_obs = exp.obs.read().concat().to_pandas()
                new_var = exp.ms["RNA"].var.read().concat().to_pandas()
            else:
                new_obs = conftest_adata.obs
                new_var = conftest_adata.var

            yield TestCase(
                experiment_path=experiment_path,
                old_anndata=old_anndata,
                new_anndata=conftest_adata,
                new_obs=new_obs,
                new_var=new_var,
                obs_schema=obs_schema,
                var_schema=var_schema,
            )


for field in fields(TestCase):

    def create_member_fixture(name):
        """Creates a ``pytest.fixture`` for a ``TestCase`` field."""

        @pytest.fixture
        def member_fixture(multiple_fixtures_with_readback):
            """
            Given a `TestCase` (multi-fixture) object, returns a getter callback to access one of that
            object's attributes. So if the `TestCase` class has a slat called `foo`, this is what will
            let tests in this file ask for `foo` as an argument passed to them.
            """
            return getattr(multiple_fixtures_with_readback, name)

        return member_fixture

    """
    Create ``pytest.fixture``s for each ``TestCase`` field. So if the `TestCase` class has slots
    `foo` and `bar` then tests in this file can ask for arguments named `foo` and/or `bar`.
    """
    globals()[field.name] = create_member_fixture(field.name)


def verify_schemas(experiment_path, obs_schema, var_schema):
    """Read {obs,var} schemas, verify they match initial versions."""
    with tiledbsoma.Experiment.open(experiment_path) as exp:
        other_obs_schema = exp.obs.schema
        other_var_schema = exp.ms["RNA"].var.schema
    assert obs_schema == other_obs_schema
    assert var_schema == other_var_schema


def verify_updates(experiment_path, obs, var, exc=False):
    """
    Calls `update_obs` and `update_var` on the experiment. Also verifies that the
    updater code didn't inadvertently modify the `obs` and `var` objects, which are
    passed into the updater by reference.
    """

    obs0 = obs.copy()
    var0 = var.copy()

    # It'd be lovely if we could do `pytest.raises(None)`, and lovelier indeed to use
    # `pytest.mark.parametrize` to make should-not-throw and should-throw cases with things like
    # `exc=[None, ValueError]`. Alas, `pytest.raises()` doesn't accept `None`. So here we
    # do a little roll-our-own keystroke-saver.
    def ctx():
        return pytest.raises(ValueError) if exc else nullcontext()

    with tiledbsoma.Experiment.open(experiment_path, "w") as exp:
        with ctx():
            tiledbsoma.io.update_obs(exp, obs)
        with ctx():
            tiledbsoma.io.update_var(exp, var, "RNA")

    assert anndata_dataframe_unmodified(obs0, obs)
    assert anndata_dataframe_unmodified(var0, var)


@pytest.mark.parametrize(
    "multiple_fixtures_with_readback", [False, True], indirect=True
)
def test_no_change(
    experiment_path, old_anndata, new_anndata, new_obs, new_var, obs_schema, var_schema
):
    verify_updates(experiment_path, new_obs, new_var)
    verify_schemas(experiment_path, obs_schema, var_schema)
    verify_obs_and_var_same(old_anndata, new_anndata)


@pytest.mark.parametrize(
    "multiple_fixtures_with_readback", [False, True], indirect=True
)
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


@pytest.mark.parametrize(
    "multiple_fixtures_with_readback", [False, True], indirect=True
)
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


@pytest.mark.parametrize(
    "multiple_fixtures_with_readback", [False, True], indirect=True
)
def test_change(experiment_path, new_obs, new_var, obs_schema, var_schema):
    new_obs["groups"] = np.arange(new_obs.shape[0], dtype=np.int16)
    new_var["vst.mean"] = np.arange(new_var.shape[0], dtype=np.int32)
    verify_updates(experiment_path, new_obs, new_var, exc=True)
    verify_schemas(experiment_path, obs_schema, var_schema)


@pytest.mark.parametrize(
    "multiple_fixtures_with_readback", [False, True], indirect=True
)
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
        verify_updates(experiment_path, new_obs2, new_var2, exc=True)
        verify_obs_and_var_same(old_anndata, new_anndata)
        verify_schemas(experiment_path, obs_schema, var_schema)
