"""
Test join-id registrations for ingesting multiple AnnData objects into a single SOMA Experiment.
"""

import math
from typing import Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import tiledbsoma.io
import tiledbsoma.io._registration as registration


def _create_anndata(
    *,
    obs_ids: Sequence[str],
    var_ids: Sequence[str],
    X_base: int,
    raw_var_ids: Optional[Sequence[str]] = None,
    measurement_name: str = "RNA",
    X_density: float = 0.3,
):
    n_obs = len(obs_ids)
    n_var = len(var_ids)

    cell_types = [["B cell", "T cell"][e % 2] for e in range(n_obs)]
    obs = pd.DataFrame(
        data={
            "obs_id": np.asarray(obs_ids),
            "cell_type": np.asarray(cell_types),
            "is_primary_data": np.asarray([True] * n_obs),
        },
        index=np.arange(n_obs).astype(str),
    )
    obs.set_index("obs_id", inplace=True)

    def _fake_mean(var_id):
        return math.sqrt(sum(ord(c) for c in var_id))

    def _make_var(arg_var_ids):
        means = [_fake_mean(var_id) for var_id in arg_var_ids]
        var = pd.DataFrame(
            data={
                "var_id": np.asarray(arg_var_ids),
                "means": np.asarray(means, dtype=np.float32),
            },
            index=np.arange(len(arg_var_ids)).astype(str),
        )
        var.set_index("var_id", inplace=True)
        return var

    def _make_X(n_obs, n_var, X_base):
        X = np.zeros((n_obs, n_var))
        for i in range(n_obs):
            for j in range(n_var):
                if (i + j) % 2 == 1:
                    X[i, j] = X_base + 10 * i + j
        return X

    var = _make_var(var_ids)
    X = _make_X(n_obs, n_var, X_base)

    adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    if raw_var_ids is not None:
        raw_var = _make_var(raw_var_ids)
        raw_X = _make_X(n_obs, len(raw_var_ids), X_base)
        raw = ad.Raw(adata, var=raw_var, X=raw_X)
        adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype, raw=raw)

    return adata


def create_h5ad(adata, path):
    adata.write_h5ad(path)
    return path


@pytest.fixture
def anndata1():
    return _create_anndata(
        obs_ids=["AAAT", "ACTG", "AGAG"],
        var_ids=["AKT1", "APOE", "ESR1", "TP53", "VEGFA"],
        raw_var_ids=["AKT1", "APOE", "ESR1", "TP53", "VEGFA", "RAW1", "RAW2"],
        X_base=100,
    )


@pytest.fixture
def anndata2():
    return _create_anndata(
        obs_ids=["CAAT", "CCTG", "CGAG"],
        var_ids=["APOE", "ESR1", "TP53", "VEGFA"],
        raw_var_ids=["APOE", "ESR1", "TP53", "VEGFA"],
        X_base=200,
    )


@pytest.fixture
def anndata3():
    return _create_anndata(
        obs_ids=["GAAT", "GCTG", "GGAG"],
        var_ids=["APOE", "EGFR", "ESR1", "TP53", "VEGFA"],
        raw_var_ids=["APOE", "EGFR", "ESR1", "TP53", "VEGFA", "RAW1", "RAW3"],
        X_base=300,
    )


@pytest.fixture
def anndata4():
    return _create_anndata(
        obs_ids=["TAAT", "TCTG", "TGAG"],
        var_ids=["AKT1", "APOE", "ESR1", "TP53", "VEGFA", "ZZZ3"],
        raw_var_ids=[
            "AKT1",
            "APOE",
            "ESR1",
            "TP53",
            "VEGFA",
            "ZZZ3",
            "RAW1",
            "RAW3",
            "RAW2",
        ],
        X_base=400,
    )


@pytest.fixture
def h5ad1(tmp_path, anndata1):
    return create_h5ad(anndata1, (tmp_path / "1.h5ad").as_posix())


@pytest.fixture
def h5ad2(tmp_path, anndata2):
    return create_h5ad(anndata2, (tmp_path / "2.h5ad").as_posix())


@pytest.fixture
def h5ad3(tmp_path, anndata3):
    return create_h5ad(anndata3, (tmp_path / "3.h5ad").as_posix())


@pytest.fixture
def h5ad4(tmp_path, anndata4):
    return create_h5ad(anndata4, (tmp_path / "4.h5ad").as_posix())


@pytest.fixture
def soma1(tmp_path, h5ad1):
    uri = (tmp_path / "soma1").as_posix()
    tiledbsoma.io.from_h5ad(uri, h5ad1, "RNA")
    return uri


def test_axis_mappings(anndata1):
    mapping = registration.AxisIDMapping.identity(10)
    assert mapping.data == tuple(range(10))

    dictionary = registration.AxisAmbientLabelMapping(
        data={"a": 10, "b": 20, "c": 30}, field_name="obs_id"
    )
    assert dictionary.id_mapping_from_values(["a", "b", "c"]).data == (10, 20, 30)
    assert dictionary.id_mapping_from_values(["c", "a"]).data == (30, 10)
    assert dictionary.id_mapping_from_values([]).data == ()

    d = registration.AxisAmbientLabelMapping.from_isolated_dataframe(
        anndata1.obs,
        index_field_name="obs_id",
    )
    assert d.id_mapping_from_values([]).data == ()
    assert d.id_mapping_from_values(["AAAT", "AGAG"]).data == (0, 2)
    keys = list(anndata1.obs.index)
    assert d.id_mapping_from_values(keys).data == tuple(range(len(keys)))


def test_isolated_anndata_mappings(anndata1):
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_anndata(
        anndata1, measurement_name="RNA"
    )
    assert rd.obs_axis.id_mapping_from_values([]).data == ()
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == (2, 1)
    assert rd.var_axes["RNA"].id_mapping_from_values(["TP53", "VEGFA"]).data == (3, 4)
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == (6, 3, 4)


def test_isolated_h5ad_mappings(h5ad1):
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_h5ad(
        h5ad1,
        measurement_name="RNA",
    )
    assert rd.obs_axis.id_mapping_from_values([]).data == ()
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == (2, 1)
    assert rd.var_axes["RNA"].id_mapping_from_values(["TP53", "VEGFA"]).data == (3, 4)
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == (6, 3, 4)


def test_isolated_soma_experiment_mappings(soma1):
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_soma_experiment(soma1)
    assert rd.obs_axis.id_mapping_from_values([]).data == ()
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == (2, 1)
    assert rd.var_axes["RNA"].id_mapping_from_values(["TP53", "VEGFA"]).data == (3, 4)
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == (6, 3, 4)


@pytest.mark.parametrize("permutation", [[0, 1, 2, 3], [2, 3, 0, 1], [3, 2, 1, 0]])
@pytest.mark.parametrize("solo_experiment_first", [True, False])
def test_multiples_without_experiment(
    tmp_path,
    h5ad1,
    h5ad2,
    h5ad3,
    h5ad4,
    permutation,
    solo_experiment_first,
):
    experiment_uri = (tmp_path / "exp").as_posix()
    h5ad_file_names = [h5ad1, h5ad2, h5ad3, h5ad4]

    if solo_experiment_first:
        # Write the first H5AD as a solo experiment. Then append the rest.
        tiledbsoma.io.from_h5ad(
            experiment_uri,
            h5ad_file_names[0],
            measurement_name="RNA",
            ingest_mode="write",
        )
        rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
            experiment_uri=experiment_uri,
            h5ad_file_names=h5ad_file_names,
            measurement_name="RNA",
            obs_field_name="obs_id",
            var_field_name="var_id",
        )

    else:
        # "Append" all the H5ADs where no experiment exists yet.
        rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
            experiment_uri=None,
            h5ad_file_names=h5ad_file_names,
            measurement_name="RNA",
            obs_field_name="obs_id",
            var_field_name="var_id",
        )

    assert rd.obs_axis.id_mapping_from_values(["AGAG", "GGAG"]).data == (2, 8)
    assert rd.var_axes["RNA"].id_mapping_from_values(["ESR1", "VEGFA"]).data == (2, 4)
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["ZZZ3", "RAW2", "TP53", "VEGFA"]
    ).data == (9, 6, 3, 4)

    assert rd.obs_axis.data == {
        "AAAT": 0,
        "ACTG": 1,
        "AGAG": 2,
        "CAAT": 3,
        "CCTG": 4,
        "CGAG": 5,
        "GAAT": 6,
        "GCTG": 7,
        "GGAG": 8,
        "TAAT": 9,
        "TCTG": 10,
        "TGAG": 11,
    }

    assert rd.var_axes["RNA"].data == {
        "AKT1": 0,
        "APOE": 1,
        "EGFR": 5,
        "ESR1": 2,
        "TP53": 3,
        "VEGFA": 4,
        "ZZZ3": 6,
    }

    assert rd.var_axes["raw"].data == {
        "AKT1": 0,
        "APOE": 1,
        "ESR1": 2,
        "TP53": 3,
        "VEGFA": 4,
        "RAW1": 5,
        "RAW2": 6,
        "EGFR": 7,
        "RAW3": 8,
        "ZZZ3": 9,
    }

    # Now do the ingestion per se.  Note that once registration is done sequentially, ingest order
    # mustn't matter, and in fact, can be done in parallel. This is why we test various permutations
    # of the ordering of the h5ad file names.
    for h5ad_file_name in [
        h5ad_file_names[permutation[0]],
        h5ad_file_names[permutation[1]],
        h5ad_file_names[permutation[2]],
        h5ad_file_names[permutation[3]],
    ]:
        tiledbsoma.io.from_h5ad(
            experiment_uri,
            h5ad_file_name,
            measurement_name="RNA",
            ingest_mode="write",
            registration_mapping=rd,
        )

    expect_obs_soma_joinids = list(range(12))
    expect_var_soma_joinids = list(range(7))

    expect_obs_obs_ids = [
        "AAAT",
        "ACTG",
        "AGAG",
        "CAAT",
        "CCTG",
        "CGAG",
        "GAAT",
        "GCTG",
        "GGAG",
        "TAAT",
        "TCTG",
        "TGAG",
    ]

    expect_var_var_ids = [
        "AKT1",
        "APOE",
        "ESR1",
        "TP53",
        "VEGFA",
        "EGFR",
        "ZZZ3",
    ]

    expect_X = pd.DataFrame(
        {
            "soma_dim_0": np.asarray(
                [
                    0,
                    0,
                    1,
                    1,
                    1,
                    2,
                    2,
                    3,
                    3,
                    4,
                    4,
                    5,
                    5,
                    6,
                    6,
                    7,
                    7,
                    7,
                    8,
                    8,
                    9,
                    9,
                    9,
                    10,
                    10,
                    10,
                    11,
                    11,
                    11,
                ],
                dtype=np.int64,
            ),
            "soma_dim_1": np.asarray(
                [
                    1,
                    3,
                    0,
                    2,
                    4,
                    1,
                    3,
                    2,
                    4,
                    1,
                    3,
                    2,
                    4,
                    3,
                    5,
                    1,
                    2,
                    4,
                    3,
                    5,
                    1,
                    3,
                    6,
                    0,
                    2,
                    4,
                    1,
                    3,
                    6,
                ],
                dtype=np.int64,
            ),
            "soma_data": np.asarray(
                [
                    101.0,
                    103.0,
                    110.0,
                    112.0,
                    114.0,
                    121.0,
                    123.0,
                    201.0,
                    203.0,
                    210.0,
                    212.0,
                    221.0,
                    223.0,
                    303.0,
                    301.0,
                    310.0,
                    312.0,
                    314.0,
                    323.0,
                    321.0,
                    401.0,
                    403.0,
                    405.0,
                    410.0,
                    412.0,
                    414.0,
                    421.0,
                    423.0,
                    425.0,
                ],
                dtype=np.float64,
            ),
        }
    )

    with tiledbsoma.Experiment.open(experiment_uri) as exp:
        obs = exp.obs.read().concat()
        var = exp.ms["RNA"].var.read().concat()

        actual_obs_soma_joinids = obs["soma_joinid"].to_pylist()
        actual_obs_obs_ids = obs["obs_id"].to_pylist()

        actual_var_soma_joinids = var["soma_joinid"].to_pylist()
        actual_var_var_ids = var["var_id"].to_pylist()

        actual_X = exp.ms["RNA"].X["data"].read().tables().concat().to_pandas()

        assert actual_obs_soma_joinids == expect_obs_soma_joinids
        assert actual_var_soma_joinids == expect_var_soma_joinids
        assert actual_obs_obs_ids == expect_obs_obs_ids
        assert actual_var_var_ids == expect_var_var_ids

        # In the happy case a simple, single-line `assert all(X == expect)` covers all of this. But
        # if the lengths don't match, the error message is useless -- "ValueError: Can only compare
        # identically-labeled DataFrame objects" -- and so for mercy to anyone debugging future
        # unit-test failures, we split out some sub-asserts for clarity.

        assert len(actual_X["soma_dim_0"].values) == len(expect_X["soma_dim_0"].values)
        assert len(actual_X["soma_dim_1"].values) == len(expect_X["soma_dim_1"].values)
        assert len(actual_X["soma_data"].values) == len(expect_X["soma_data"].values)
        assert all(actual_X.dtypes == expect_X.dtypes)
        assert all(actual_X == expect_X)


def test_multiples_with_experiment(soma1, h5ad2, h5ad3, h5ad4):
    rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
        experiment_uri=soma1,
        h5ad_file_names=[h5ad2, h5ad3, h5ad4],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "GGAG"]).data == (2, 8)
    assert rd.var_axes["RNA"].id_mapping_from_values(["ESR1", "VEGFA"]).data == (2, 4)
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["ZZZ3", "RAW2", "TP53", "VEGFA"]
    ).data == (9, 6, 3, 4)

    assert rd.obs_axis.data == {
        "AAAT": 0,
        "ACTG": 1,
        "AGAG": 2,
        "CAAT": 3,
        "CCTG": 4,
        "CGAG": 5,
        "GAAT": 6,
        "GCTG": 7,
        "GGAG": 8,
        "TAAT": 9,
        "TCTG": 10,
        "TGAG": 11,
    }

    assert rd.var_axes["RNA"].data == {
        "AKT1": 0,
        "APOE": 1,
        "EGFR": 5,
        "ESR1": 2,
        "TP53": 3,
        "VEGFA": 4,
        "ZZZ3": 6,
    }

    assert rd.var_axes["raw"].data == {
        "AKT1": 0,
        "APOE": 1,
        "ESR1": 2,
        "TP53": 3,
        "VEGFA": 4,
        "RAW1": 5,
        "RAW2": 6,
        "EGFR": 7,
        "RAW3": 8,
        "ZZZ3": 9,
    }


def test_append_dataframes_only_with_experiment(soma1, h5ad2):
    rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
        experiment_uri=soma1,
        h5ad_file_names=[h5ad2],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    adata2 = ad.read_h5ad(h5ad2)

    with tiledbsoma.Experiment.open(soma1, "w") as exp1:
        tiledbsoma.io.append_obs(
            exp1.obs.uri,
            adata2.obs,
            registration_mapping=rd,
        )

    expect_obs_soma_joinids = list(range(6))
    expect_var_soma_joinids = list(range(5))

    expect_obs_obs_ids = [
        "AAAT",
        "ACTG",
        "AGAG",
        "CAAT",
        "CCTG",
        "CGAG",
    ]

    expect_var_var_ids = [
        "AKT1",
        "APOE",
        "ESR1",
        "TP53",
        "VEGFA",
    ]

    with tiledbsoma.Experiment.open(soma1) as exp:
        obs = exp.obs.read().concat()
        var = exp.ms["RNA"].var.read().concat()

        actual_obs_soma_joinids = obs["soma_joinid"].to_pylist()
        actual_obs_obs_ids = obs["obs_id"].to_pylist()

        actual_var_soma_joinids = var["soma_joinid"].to_pylist()
        actual_var_var_ids = var["var_id"].to_pylist()

        assert actual_obs_soma_joinids == expect_obs_soma_joinids
        assert actual_var_soma_joinids == expect_var_soma_joinids
        assert actual_obs_obs_ids == expect_obs_obs_ids
        assert actual_var_var_ids == expect_var_var_ids
