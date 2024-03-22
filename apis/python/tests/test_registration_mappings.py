"""
Test join-id registrations for ingesting multiple AnnData objects into a single SOMA Experiment.
"""

import math
import tempfile
from typing import Optional, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import tiledbsoma.io
import tiledbsoma.io._registration as registration
from tiledbsoma._util import verify_obs_var


def _create_anndata(
    *,
    obs_ids: Sequence[str],
    var_ids: Sequence[str],
    obs_field_name: str,
    var_field_name: str,
    X_value_base: int,
    measurement_name: str,
    raw_var_ids: Optional[Sequence[str]] = None,
    X_density: float = 0.3,
):
    n_obs = len(obs_ids)
    n_var = len(var_ids)

    cell_types = [["B cell", "T cell"][e % 2] for e in range(n_obs)]
    obs = pd.DataFrame(
        data={
            obs_field_name: np.asarray(obs_ids),
            "cell_type": np.asarray(cell_types),
            "is_primary_data": np.asarray([True] * n_obs),
        },
        index=np.arange(n_obs).astype(str),
    )
    obs.set_index(obs_field_name, inplace=True)

    def _fake_mean(var_id):
        return math.sqrt(sum(ord(c) for c in var_id))

    def _make_var(arg_var_ids):
        means = [_fake_mean(var_id) for var_id in arg_var_ids]
        var = pd.DataFrame(
            data={
                var_field_name: np.asarray(arg_var_ids),
                "means": np.asarray(means, dtype=np.float32),
            },
            index=np.arange(len(arg_var_ids)).astype(str),
        )
        var.set_index(var_field_name, inplace=True)
        return var

    def _make_X(n_obs, n_var, X_value_base):
        X = np.zeros((n_obs, n_var))
        for i in range(n_obs):
            for j in range(n_var):
                if (i + j) % 2 == 1:
                    X[i, j] = X_value_base + 10 * i + j
        return X

    var = _make_var(var_ids)
    X = _make_X(n_obs, n_var, X_value_base)

    adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype)

    if raw_var_ids is not None:
        raw_var = _make_var(raw_var_ids)
        raw_X = _make_X(n_obs, len(raw_var_ids), X_value_base)
        raw = ad.Raw(adata, var=raw_var, X=raw_X)
        adata = ad.AnnData(X=X, obs=obs, var=var, dtype=X.dtype, raw=raw)

    return adata


def create_h5ad(adata, path):
    adata.write_h5ad(path)
    return path


# This is the central data generator for this file. It makes one of four
# datasets, varying cell and gene IDs.
def create_anndata_canned(which: int, obs_field_name: str, var_field_name: str):
    if which == 1:
        obs_ids = ["AAAT", "ACTG", "AGAG"]
        var_ids = ["AKT1", "APOE", "ESR1", "TP53", "VEGFA"]
        raw_var_ids = ["AKT1", "APOE", "ESR1", "TP53", "VEGFA", "RAW1", "RAW2"]
        X_value_base = 100

    elif which == 2:
        obs_ids = ["CAAT", "CCTG", "CGAG"]
        var_ids = ["APOE", "ESR1", "TP53", "VEGFA"]
        raw_var_ids = ["APOE", "ESR1", "TP53", "VEGFA"]
        X_value_base = 200

    elif which == 3:
        obs_ids = ["GAAT", "GCTG", "GGAG"]
        var_ids = ["APOE", "EGFR", "ESR1", "TP53", "VEGFA"]
        raw_var_ids = ["APOE", "EGFR", "ESR1", "TP53", "VEGFA", "RAW1", "RAW3"]
        X_value_base = 300

    elif which == 4:
        obs_ids = ["TAAT", "TCTG", "TGAG"]
        var_ids = ["AKT1", "APOE", "ESR1", "TP53", "VEGFA", "ZZZ3"]
        raw_var_ids = [
            "AKT1",
            "APOE",
            "ESR1",
            "TP53",
            "VEGFA",
            "ZZZ3",
            "RAW1",
            "RAW3",
            "RAW2",
        ]
        X_value_base = 400

    else:
        raise Exception(f"create_anndata_canned takes 1..4; got {which}")

    return _create_anndata(
        obs_ids=obs_ids,
        var_ids=var_ids,
        raw_var_ids=raw_var_ids,
        X_value_base=X_value_base,
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )


def create_h5ad_canned(which: int, obs_field_name: str, var_field_name: str):
    tmp_path = tempfile.TemporaryDirectory()
    anndata = create_anndata_canned(which, obs_field_name, var_field_name)
    return create_h5ad(
        anndata,
        (tmp_path.name + f"{which}.h5ad"),
    )


def create_soma_canned(which: int, obs_field_name, var_field_name):
    tmp_path = tempfile.TemporaryDirectory()
    h5ad = create_h5ad_canned(which, obs_field_name, var_field_name)
    uri = tmp_path.name + f"soma{which}"
    tiledbsoma.io.from_h5ad(uri, h5ad, "measname")
    return uri


@pytest.fixture
def anndata_larger():
    return _create_anndata(
        obs_ids=["id_%08d" % e for e in range(1000)],
        var_ids=["AKT1", "APOE", "ESR1", "TP53", "VEGFA", "ZZZ3"],
        X_value_base=0,
        measurement_name="measname",
        obs_field_name="cell_id",
        var_field_name="gene_id",
    )


@pytest.fixture
def soma_larger(anndata_larger):
    tmp_path = tempfile.TemporaryDirectory()
    uri = tmp_path.name + "soma-larger"
    tiledbsoma.io.from_anndata(uri, anndata_larger, "measname")
    return uri


@pytest.mark.parametrize(
    "args",
    [
        # SOMA ID column is to be obs_id, and it is the Pandas index named "obs_id"
        {
            "do_set_index": True,
            "index_name_to_set": "obs_id",
            "do_rename_axis": False,
            "axis_name_to_set": None,
            "registration_index_column_name": "obs_id",
            "expected_signature": {"obs_id": "string", "alt_id": "string"},
        },
        # SOMA ID column is to be obs_id, and it is the Pandas index named "index"
        {
            "do_set_index": True,
            "index_name_to_set": "obs_id",
            "do_rename_axis": True,
            "axis_name_to_set": "index",
            "registration_index_column_name": "obs_id",
            "expected_signature": {"obs_id": "string", "alt_id": "string"},
        },
        # SOMA ID column is to be obs_id, and it is the Pandas unnamed index
        {
            "do_set_index": True,
            "index_name_to_set": "obs_id",
            "do_rename_axis": True,
            "axis_name_to_set": None,
            "registration_index_column_name": "obs_id",
            "expected_signature": {"obs_id": "string", "alt_id": "string"},
        },
        # SOMA ID column is to be obs_id, and the Pandas index is named something else
        {
            "do_set_index": True,
            "index_name_to_set": "alt_id",
            "do_rename_axis": False,
            "axis_name_to_set": None,
            "registration_index_column_name": "obs_id",
            "expected_signature": {"alt_id": "string", "obs_id": "string"},
        },
        # SOMA ID column is to be obs_id, and the Pandas index is unnamed
        {
            "do_set_index": True,
            "index_name_to_set": "alt_id",
            "do_rename_axis": True,
            "axis_name_to_set": None,
            "registration_index_column_name": "obs_id",
            "expected_signature": {"obs_id": "string"},
        },
        # SOMA ID column is to be obs_id, and the Pandas index is named "index"
        {
            "do_set_index": True,
            "index_name_to_set": "alt_id",
            "do_rename_axis": True,
            "axis_name_to_set": "index",
            "registration_index_column_name": "obs_id",
            "expected_signature": {"obs_id": "string"},
        },
        # SOMA ID column is to be obs_id, and the Pandas index is implicitized integers
        {
            "do_set_index": False,
            "index_name_to_set": None,
            "do_rename_axis": False,
            "axis_name_to_set": None,
            "registration_index_column_name": "obs_id",
            "expected_signature": {"alt_id": "string", "obs_id": "string"},
        },
        # SOMA ID column is to be alt_id, and it is the Pandas index named "alt_id"
        {
            "do_set_index": True,
            "index_name_to_set": "alt_id",
            "do_rename_axis": False,
            "axis_name_to_set": None,
            "registration_index_column_name": "alt_id",
            "expected_signature": {"alt_id": "string", "obs_id": "string"},
        },
        # SOMA ID column is to be alt_id, and it is the Pandas index named "index"
        {
            "do_set_index": True,
            "index_name_to_set": "alt_id",
            "do_rename_axis": True,
            "axis_name_to_set": "index",
            "registration_index_column_name": "alt_id",
            "expected_signature": {"alt_id": "string", "obs_id": "string"},
        },
        # SOMA ID column is to be alt_id, and it is the Pandas unnamed index
        {
            "do_set_index": True,
            "index_name_to_set": "alt_id",
            "do_rename_axis": True,
            "axis_name_to_set": None,
            "registration_index_column_name": "alt_id",
            "expected_signature": {"alt_id": "string", "obs_id": "string"},
        },
        # SOMA ID column is to be alt_id, and the Pandas index is named something else
        {
            "do_set_index": True,
            "index_name_to_set": "obs_id",
            "do_rename_axis": False,
            "axis_name_to_set": None,
            "registration_index_column_name": "alt_id",
            "expected_signature": {"obs_id": "string", "alt_id": "string"},
        },
        # SOMA ID column is to be alt_id, and the Pandas index is unnamed
        {
            "do_set_index": True,
            "index_name_to_set": "obs_id",
            "do_rename_axis": True,
            "axis_name_to_set": None,
            "registration_index_column_name": "alt_id",
            "expected_signature": {"alt_id": "string"},
        },
        # SOMA ID column is to be alt_id, and the Pandas index is named "index"
        {
            "do_set_index": True,
            "index_name_to_set": "obs_id",
            "do_rename_axis": True,
            "axis_name_to_set": "index",
            "registration_index_column_name": "alt_id",
            "expected_signature": {"alt_id": "string"},
        },
        # SOMA ID column is to be alt_id, and the Pandas index is implicitized integers
        {
            "do_set_index": False,
            "index_name_to_set": None,
            "do_rename_axis": False,
            "axis_name_to_set": None,
            "registration_index_column_name": "alt_id",
            "expected_signature": {"alt_id": "string", "obs_id": "string"},
        },
    ],
)
def test_pandas_indexing(args):
    """
    The index-column name for registration can take a variety of forms.
    This test exercises all of them.
    """

    df = pd.DataFrame(
        data={
            "soma_joinid": np.arange(3, dtype=np.int64),
            "alt_id": ["A", "C", "G"],
            "obs_id": ["AT", "CT", "GT"],
        }
    )
    if args["do_set_index"]:
        df.set_index(args["index_name_to_set"], inplace=True)
    if args["do_rename_axis"]:
        df.rename_axis(args["axis_name_to_set"], inplace=True)

    actual_signature = registration.signatures._string_dict_from_pandas_dataframe(
        df,
        args["registration_index_column_name"],
    )
    assert actual_signature == args["expected_signature"]


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
def test_axis_mappings(obs_field_name, var_field_name):
    anndata1 = create_anndata_canned(1, obs_field_name, var_field_name)
    mapping = registration.AxisIDMapping.identity(10)
    assert mapping.data == tuple(range(10))

    dictionary = registration.AxisAmbientLabelMapping(
        data={"a": 10, "b": 20, "c": 30},
        field_name=obs_field_name,
    )
    assert dictionary.id_mapping_from_values(["a", "b", "c"]).data == (10, 20, 30)
    assert dictionary.id_mapping_from_values(["c", "a"]).data == (30, 10)
    assert dictionary.id_mapping_from_values([]).data == ()

    d = registration.AxisAmbientLabelMapping.from_isolated_dataframe(
        anndata1.obs,
        index_field_name=obs_field_name,
    )
    assert d.id_mapping_from_values([]).data == ()
    assert d.id_mapping_from_values(["AAAT", "AGAG"]).data == (0, 2)
    keys = list(anndata1.obs.index)
    assert d.id_mapping_from_values(keys).data == tuple(range(len(keys)))


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
def test_isolated_anndata_mappings(obs_field_name, var_field_name):
    anndata1 = create_anndata_canned(1, obs_field_name, var_field_name)
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_anndata(
        anndata1, measurement_name="measname"
    )
    assert rd.obs_axis.id_mapping_from_values([]).data == ()
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == (2, 1)
    assert rd.var_axes["measname"].id_mapping_from_values(["TP53", "VEGFA"]).data == (
        3,
        4,
    )
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == (6, 3, 4)


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
def test_isolated_h5ad_mappings(obs_field_name, var_field_name):
    h5ad1 = create_h5ad_canned(1, obs_field_name, var_field_name)
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_h5ad(
        h5ad1,
        measurement_name="measname",
    )
    assert rd.obs_axis.id_mapping_from_values([]).data == ()
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == (2, 1)
    assert rd.var_axes["measname"].id_mapping_from_values(["TP53", "VEGFA"]).data == (
        3,
        4,
    )
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == (6, 3, 4)


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
def test_isolated_soma_experiment_mappings(obs_field_name, var_field_name):
    soma1 = create_soma_canned(1, obs_field_name, var_field_name)
    rd = registration.ExperimentAmbientLabelMapping.from_isolated_soma_experiment(
        soma1, obs_field_name=obs_field_name, var_field_name=var_field_name
    )
    assert rd.obs_axis.id_mapping_from_values([]).data == ()
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == (2, 1)
    assert rd.var_axes["measname"].id_mapping_from_values(["TP53", "VEGFA"]).data == (
        3,
        4,
    )
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == (6, 3, 4)


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
@pytest.mark.parametrize("permutation", [[0, 1, 2, 3], [2, 3, 0, 1], [3, 2, 1, 0]])
@pytest.mark.parametrize("solo_experiment_first", [True, False])
def test_multiples_without_experiment(
    tmp_path,
    obs_field_name,
    var_field_name,
    permutation,
    solo_experiment_first,
):
    h5ad1 = create_h5ad_canned(1, obs_field_name, var_field_name)
    h5ad2 = create_h5ad_canned(2, obs_field_name, var_field_name)
    h5ad3 = create_h5ad_canned(3, obs_field_name, var_field_name)
    h5ad4 = create_h5ad_canned(4, obs_field_name, var_field_name)

    experiment_uri = (tmp_path / "exp").as_posix()
    h5ad_file_names = [h5ad1, h5ad2, h5ad3, h5ad4]

    if solo_experiment_first:
        # Write the first H5AD as a solo experiment. Then append the rest.
        tiledbsoma.io.from_h5ad(
            experiment_uri,
            h5ad_file_names[0],
            measurement_name="measname",
            ingest_mode="write",
        )
        rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
            experiment_uri=experiment_uri,
            h5ad_file_names=h5ad_file_names,
            measurement_name="measname",
            obs_field_name=obs_field_name,
            var_field_name=var_field_name,
        )

    else:
        # "Append" all the H5ADs where no experiment exists yet.
        rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
            experiment_uri=None,
            h5ad_file_names=h5ad_file_names,
            measurement_name="measname",
            obs_field_name=obs_field_name,
            var_field_name=var_field_name,
        )

    assert rd.obs_axis.id_mapping_from_values(["AGAG", "GGAG"]).data == (2, 8)
    assert rd.var_axes["measname"].id_mapping_from_values(["ESR1", "VEGFA"]).data == (
        2,
        4,
    )
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

    assert rd.var_axes["measname"].data == {
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
            measurement_name="measname",
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
        var = exp.ms["measname"].var.read().concat()

        actual_obs_soma_joinids = obs["soma_joinid"].to_pylist()
        actual_obs_obs_ids = obs[obs_field_name].to_pylist()

        actual_var_soma_joinids = var["soma_joinid"].to_pylist()
        actual_var_var_ids = var[var_field_name].to_pylist()

        actual_X = exp.ms["measname"].X["data"].read().tables().concat().to_pandas()

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

        X = exp.ms["measname"].X["data"]
        assert X.used_shape() == ((0, 11), (0, 6))
        assert X.non_empty_domain() == ((0, 11), (0, 6))


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
def test_multiples_with_experiment(obs_field_name, var_field_name):
    soma1 = create_soma_canned(1, obs_field_name, var_field_name)
    h5ad2 = create_h5ad_canned(2, obs_field_name, var_field_name)
    h5ad3 = create_h5ad_canned(3, obs_field_name, var_field_name)
    h5ad4 = create_h5ad_canned(4, obs_field_name, var_field_name)

    rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
        experiment_uri=soma1,
        h5ad_file_names=[h5ad2, h5ad3, h5ad4],
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "GGAG"]).data == (2, 8)
    assert rd.var_axes["measname"].id_mapping_from_values(["ESR1", "VEGFA"]).data == (
        2,
        4,
    )
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

    assert rd.var_axes["measname"].data == {
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


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
def test_append_items_with_experiment(obs_field_name, var_field_name):
    soma1 = create_soma_canned(1, obs_field_name, var_field_name)
    h5ad2 = create_h5ad_canned(2, obs_field_name, var_field_name)
    rd = registration.ExperimentAmbientLabelMapping.from_h5ad_appends_on_experiment(
        experiment_uri=soma1,
        h5ad_file_names=[h5ad2],
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )

    adata2 = ad.read_h5ad(h5ad2)

    original = adata2.copy()

    with tiledbsoma.Experiment.open(soma1, "w") as exp1:
        tiledbsoma.io.append_obs(
            exp1,
            adata2.obs,
            registration_mapping=rd,
        )

        tiledbsoma.io.append_var(
            exp1,
            adata2.var,
            measurement_name="measname",
            registration_mapping=rd,
        )

        tiledbsoma.io.append_X(
            exp1,
            adata2.X,
            measurement_name="measname",
            X_layer_name="data",
            obs_ids=list(adata2.obs.index),
            var_ids=list(adata2.var.index),
            registration_mapping=rd,
        )

    verify_obs_var(original, adata2)

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
        var = exp.ms["measname"].var.read().concat()
        exp.ms["measname"].X["data"].read().tables().concat()

        actual_obs_soma_joinids = obs["soma_joinid"].to_pylist()
        actual_obs_obs_ids = obs[obs_field_name].to_pylist()

        actual_var_soma_joinids = var["soma_joinid"].to_pylist()
        actual_var_var_ids = var[var_field_name].to_pylist()

        assert actual_obs_soma_joinids == expect_obs_soma_joinids
        assert actual_var_soma_joinids == expect_var_soma_joinids
        assert actual_obs_obs_ids == expect_obs_obs_ids
        assert actual_var_var_ids == expect_var_var_ids

        actual_X = exp.ms["measname"].X["data"].read().tables().concat().to_pandas()

        expect_X = pd.DataFrame(
            {
                "soma_dim_0": np.asarray(
                    [0, 0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5], dtype=np.int64
                ),
                "soma_dim_1": np.asarray(
                    [1, 3, 0, 2, 4, 1, 3, 2, 4, 1, 3, 2, 4], dtype=np.int64
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
                    ],
                    dtype=np.float64,
                ),
            }
        )

        assert all(actual_X == expect_X)


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
@pytest.mark.parametrize("use_same_cells", [True, False])
def test_append_with_disjoint_measurements(
    tmp_path, obs_field_name, var_field_name, use_same_cells
):
    anndata1 = create_anndata_canned(1, obs_field_name, var_field_name)
    anndata4 = create_anndata_canned(4, obs_field_name, var_field_name)
    soma_uri = tmp_path.as_posix()

    tiledbsoma.io.from_anndata(soma_uri, anndata1, measurement_name="one")

    with tiledbsoma.open(soma_uri, "w") as exp:
        exp.ms.add_new_collection("two", kind=tiledbsoma.Measurement)

    anndata2 = anndata1 if use_same_cells else anndata4

    original = anndata2.copy()

    rd = tiledbsoma.io.register_anndatas(
        soma_uri,
        [anndata2],
        measurement_name="two",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )

    tiledbsoma.io.from_anndata(
        soma_uri,
        anndata2,
        measurement_name="two",
        registration_mapping=rd,
    )

    verify_obs_var(original, anndata2)

    # exp/obs, use_same_cells=True:                       exp/obs, use_same_cells=False:
    #    soma_joinid obs_id cell_type  is_primary_data       soma_joinid obs_id cell_type  is_primary_data
    # 0            0   AAAT    B cell                1    0            0   AAAT    B cell                1
    # 1            1   ACTG    T cell                1    1            1   ACTG    T cell                1
    # 2            2   AGAG    B cell                1    2            2   AGAG    B cell                1
    #                                                     3            3   TAAT    B cell                1
    #                                                     4            4   TCTG    T cell                1
    #                                                     5            5   TGAG    B cell                1
    #
    # exp/ms/one/var, use_same_cells=True:                exp/ms/one/var, use_same_cells=False:
    #    soma_joinid var_id      means                       soma_joinid var_id      means
    # 0            0   AKT1  16.522711                    0            0   AKT1  16.522711
    # 1            1   APOE  17.117243                    1            1   APOE  17.117243
    # 2            2   ESR1  16.822603                    2            2   ESR1  16.822603
    # 3            3   TP53  16.370705                    3            3   TP53  16.370705
    # 4            4  VEGFA  19.000000                    4            4  VEGFA  19.000000
    #
    # exp/ms/two/var, use_same_cells=True:                exp/ms/two/var, use_same_cells=False:
    #    soma_joinid var_id      means                       soma_joinid var_id      means
    # 0            0   AKT1  16.522711                    0            0   AKT1  16.522711
    # 1            1   APOE  17.117243                    1            1   APOE  17.117243
    # 2            2   ESR1  16.822603                    2            2   ESR1  16.822603
    # 3            3   TP53  16.370705                    3            3   TP53  16.370705
    # 4            4  VEGFA  19.000000                    4            4  VEGFA  19.000000
    #                                                     5            5   ZZZ3  17.916473
    #
    # exp/ms/one/X/data, use_same_cells=True:             exp/ms/one/X/data, use_same_cells=False:
    #    soma_dim_0  soma_dim_1  soma_data                   soma_dim_0  soma_dim_1  soma_data
    # 0           0           1      101.0                0           0           1      101.0
    # 1           0           3      103.0                1           0           3      103.0
    # 2           1           0      110.0                2           1           0      110.0
    # 3           1           2      112.0                3           1           2      112.0
    # 4           1           4      114.0                4           1           4      114.0
    # 5           2           1      121.0                5           2           1      121.0
    # 6           2           3      123.0                6           2           3      123.0
    #
    # exp/ms/two/X/data, use_same_cells=True:             exp/ms/two/X/data, use_same_cells=False:
    #    soma_dim_0  soma_dim_1  soma_data                   soma_dim_0  soma_dim_1  soma_data
    # 0           0           1      101.0                0           3           1      401.0
    # 1           0           3      103.0                1           3           3      403.0
    # 2           1           0      110.0                2           3           5      405.0
    # 3           1           2      112.0                3           4           0      410.0
    # 4           1           4      114.0                4           4           2      412.0
    # 5           2           1      121.0                5           4           4      414.0
    # 6           2           3      123.0                6           5           1      421.0
    #                                                     7           5           3      423.0
    #                                                     8           5           5      425.0

    if use_same_cells:
        # Data for the same cells, ingested into two different measurements:
        # obs should not extend.
        expect_obs_soma_joinids = list(range(3))
        expect_obs_obs_ids = [
            "AAAT",
            "ACTG",
            "AGAG",
        ]
    else:
        # Data for different cells, ingested into two different measurements:
        # obs should extend.
        expect_obs_soma_joinids = list(range(6))
        expect_obs_obs_ids = [
            "AAAT",
            "ACTG",
            "AGAG",
            "TAAT",
            "TCTG",
            "TGAG",
        ]

    # anndata1 cell IDs go into measurement one's var.
    expect_var_one_soma_joinids = list(range(5))
    expect_var_one_var_ids = [
        "AKT1",
        "APOE",
        "ESR1",
        "TP53",
        "VEGFA",
    ]

    # anndata2 cell IDs go into measurement two's var.
    if use_same_cells:
        expect_var_two_soma_joinids = list(range(5))
        expect_var_two_var_ids = [
            "AKT1",
            "APOE",
            "ESR1",
            "TP53",
            "VEGFA",
        ]
    else:
        expect_var_two_soma_joinids = list(range(6))
        expect_var_two_var_ids = [
            "AKT1",
            "APOE",
            "ESR1",
            "TP53",
            "VEGFA",
            "ZZZ3",
        ]

    with tiledbsoma.Experiment.open(soma_uri) as exp:
        obs = exp.obs.read().concat()
        var_one = exp.ms["one"].var.read().concat()
        var_two = exp.ms["two"].var.read().concat()

        actual_obs_soma_joinids = obs["soma_joinid"].to_pylist()
        actual_obs_obs_ids = obs[obs_field_name].to_pylist()

        actual_var_one_soma_joinids = var_one["soma_joinid"].to_pylist()
        actual_var_one_var_ids = var_one[var_field_name].to_pylist()

        actual_var_two_soma_joinids = var_two["soma_joinid"].to_pylist()
        actual_var_two_var_ids = var_two[var_field_name].to_pylist()

        assert actual_obs_soma_joinids == expect_obs_soma_joinids
        assert actual_obs_obs_ids == expect_obs_obs_ids

        assert actual_var_one_soma_joinids == expect_var_one_soma_joinids
        assert actual_var_one_var_ids == expect_var_one_var_ids

        assert actual_var_two_soma_joinids == expect_var_two_soma_joinids
        assert actual_var_two_var_ids == expect_var_two_var_ids

        actual_X_one = exp.ms["one"].X["data"].read().tables().concat().to_pandas()
        actual_X_two = exp.ms["two"].X["data"].read().tables().concat().to_pandas()

        expect_X_one = pd.DataFrame(
            {
                "soma_dim_0": np.asarray([0, 0, 1, 1, 1, 2, 2], dtype=np.int64),
                "soma_dim_1": np.asarray([1, 3, 0, 2, 4, 1, 3], dtype=np.int64),
                "soma_data": np.asarray(
                    [
                        101.0,
                        103.0,
                        110.0,
                        112.0,
                        114.0,
                        121.0,
                        123.0,
                    ],
                    dtype=np.float64,
                ),
            }
        )

        if use_same_cells:
            expect_X_two = pd.DataFrame(
                {
                    "soma_dim_0": np.asarray([0, 0, 1, 1, 1, 2, 2], dtype=np.int64),
                    "soma_dim_1": np.asarray([1, 3, 0, 2, 4, 1, 3], dtype=np.int64),
                    "soma_data": np.asarray(
                        [
                            101.0,
                            103.0,
                            110.0,
                            112.0,
                            114.0,
                            121.0,
                            123.0,
                        ],
                        dtype=np.float64,
                    ),
                }
            )
        else:
            expect_X_two = pd.DataFrame(
                {
                    "soma_dim_0": np.asarray(
                        [3, 3, 3, 4, 4, 4, 5, 5, 5], dtype=np.int64
                    ),
                    "soma_dim_1": np.asarray(
                        [1, 3, 5, 0, 2, 4, 1, 3, 5], dtype=np.int64
                    ),
                    "soma_data": np.asarray(
                        [
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

        assert all(actual_X_one == expect_X_one)
        assert all(actual_X_two == expect_X_two)


@pytest.mark.parametrize("use_small_buffer", [False, True])
def test_registration_with_batched_reads(tmp_path, soma_larger, use_small_buffer):
    # This tests https://github.com/single-cell-data/TileDB-SOMA/pull/2112.
    # We check that it takes more than one batch iteration to read all of
    # obs, and then we check that the registration got all the obs IDs.

    context = None
    if use_small_buffer:
        context = tiledbsoma.SOMATileDBContext(
            tiledb_config={
                "soma.init_buffer_bytes": 2048,
            }
        )

    with tiledbsoma.Experiment.open(soma_larger, context=context) as exp:
        assert exp.obs.count == 1000

        nbatch = 0
        for batch in exp.obs.read():
            nbatch += 1
        if use_small_buffer:
            assert nbatch > 1

    rd = registration.ExperimentAmbientLabelMapping.from_isolated_soma_experiment(
        soma_larger,
        context=context,
        obs_field_name="cell_id",
        var_field_name="gene_id",
    )

    assert len(rd.obs_axis.data) == 1000


def test_ealm_expose():
    """Checks that this is exported from tiledbsoma.io._registration"""
    # All we want to check is that the import doesn't throw. Job done. Period.
    # However, the pre-commit hook will strip out this import statement as "unused".
    # So, assert something.
    assert tiledbsoma.io.ExperimentAmbientLabelMapping is not None
