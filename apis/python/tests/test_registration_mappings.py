import math
import random
from typing import Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import tiledbsoma.io
import tiledbsoma.io.registration as registration

"""
Test join-id registrations for ingesting multiple AnnData objects into a single SOMA Experiment.
"""


def _create_anndata(
    obs_ids: Sequence[str],
    var_ids: Sequence[str],
    *,
    raw_var_ids: Optional[Sequence[str]] = None,
    measurement_name: str = "RNA",
    X_density: float = 0.3,
):
    n_obs = len(obs_ids)
    n_var = len(var_ids)

    cell_types = [random.choice(["B cell", "T cell"]) for e in range(n_obs)]
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

    def _make_X(n_obs, n_var):
        X = np.zeros((n_obs, n_var))
        for i in range(n_obs):
            for j in range(n_var):
                if random.uniform(0, 1) < X_density:
                    X[i, j] = int(random.gauss(3, 1) ** 2)
        return X

    var = _make_var(var_ids)
    X = _make_X(n_obs, n_var)

    adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    if raw_var_ids is not None:
        raw_var = _make_var(raw_var_ids)
        raw_X = _make_X(n_obs, len(raw_var_ids))
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
    )


@pytest.fixture
def anndata2():
    return _create_anndata(
        obs_ids=["CAAT", "CCTG", "CGAG"],
        var_ids=["APOE", "ESR1", "TP53", "VEGFA"],
        raw_var_ids=["APOE", "ESR1", "TP53", "VEGFA"],
    )


@pytest.fixture
def anndata3():
    return _create_anndata(
        obs_ids=["GAAT", "GCTG", "GGAG"],
        var_ids=["APOE", "EGFR", "ESR1", "TP53", "VEGFA"],
        raw_var_ids=["APOE", "EGFR", "ESR1", "TP53", "VEGFA", "RAW1", "RAW3"],
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
    assert mapping.data == list(range(10))

    dictionary = registration.AxisAmbientLabelMapping(
        {"a": 10, "b": 20, "c": 30}, "obs_id"
    )
    assert dictionary.id_mapping_from_values(["a", "b", "c"]).data == [10, 20, 30]
    assert dictionary.id_mapping_from_values(["c", "a"]).data == [30, 10]
    assert dictionary.id_mapping_from_values([]).data == []

    d = registration.AxisAmbientLabelMapping.from_isolated_dataframe(
        anndata1.obs,
        index_field_name="obs_id",
    )
    assert d.id_mapping_from_values([]).data == []
    assert d.id_mapping_from_values(["AAAT", "AGAG"]).data == [0, 2]
    keys = list(anndata1.obs.index)
    assert d.id_mapping_from_values(keys).data == list(range(len(keys)))


def test_isolated_anndata_mappings(anndata1):
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_anndata(
        anndata1, measurement_name="RNA"
    )
    assert rd.obs_axis.id_mapping_from_values([]).data == []
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == [2, 1]
    assert rd.var_axes["RNA"].id_mapping_from_values(["TP53", "VEGFA"]).data == [3, 4]
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == [6, 3, 4]


def test_isolated_h5ad_mappings(h5ad1):
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_h5ad(
        h5ad1,
        measurement_name="RNA",
    )
    assert rd.obs_axis.id_mapping_from_values([]).data == []
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == [2, 1]
    assert rd.var_axes["RNA"].id_mapping_from_values(["TP53", "VEGFA"]).data == [3, 4]
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == [6, 3, 4]


def test_isolated_soma_experiment_mappings(soma1):
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_soma_experiment(soma1)
    assert rd.obs_axis.id_mapping_from_values([]).data == []
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == [2, 1]
    assert rd.var_axes["RNA"].id_mapping_from_values(["TP53", "VEGFA"]).data == [3, 4]
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == [6, 3, 4]


def test_multiples_without_experiment(h5ad1, h5ad2, h5ad3, h5ad4):
    rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
        experiment_uri=None,
        h5ad_file_names=[h5ad1, h5ad2, h5ad3, h5ad4],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "GGAG"]).data == [2, 8]
    assert rd.var_axes["RNA"].id_mapping_from_values(["ESR1", "VEGFA"]).data == [2, 4]
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["ZZZ3", "RAW2", "TP53", "VEGFA"]
    ).data == [9, 6, 3, 4]

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


def test_multiples_with_experiment(soma1, h5ad2, h5ad3, h5ad4):
    rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
        experiment_uri=soma1,
        h5ad_file_names=[h5ad2, h5ad3, h5ad4],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "GGAG"]).data == [2, 8]
    assert rd.var_axes["RNA"].id_mapping_from_values(["ESR1", "VEGFA"]).data == [2, 4]
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["ZZZ3", "RAW2", "TP53", "VEGFA"]
    ).data == [9, 6, 3, 4]

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
