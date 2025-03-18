"""
Test join-id registrations for ingesting multiple AnnData objects into a single SOMA Experiment.
"""

from __future__ import annotations

import math
import tempfile
from contextlib import nullcontext
from typing import List, Sequence, Tuple, Union

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import tiledbsoma.io
import tiledbsoma.io._registration as registration
from tiledbsoma.io import conversions

from ._util import assert_adata_equal


def _create_anndata(
    *,
    obs_ids: Sequence[str],
    var_ids: Sequence[str],
    obs_field_name: str,
    var_field_name: str,
    X_value_base: int,
    raw_var_ids: Sequence[str] | None = None,
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


def create_h5ad(conftest_pbmc_small, path):
    conftest_pbmc_small.write_h5ad(path)
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

    elif which == 8:
        obs_ids = ["TAAT", "TCTG", "TGAG", "DUP", "DUP"]
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
        X_value_base = 800

    elif which == 9:
        obs_ids = ["TAAT", "TCTG", "TGAG"]
        var_ids = ["AKT1", "DUP", "ESR1", "TP53", "DUP", "ZZZ3"]
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
        X_value_base = 900

    else:
        raise Exception(f"create_anndata_canned got unrecognized which={which}")

    return _create_anndata(
        obs_ids=obs_ids,
        var_ids=var_ids,
        raw_var_ids=raw_var_ids,
        X_value_base=X_value_base,
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )


def create_h5ad_canned(which: int, obs_field_name: str, var_field_name: str):
    tmp_path = tempfile.TemporaryDirectory(prefix="create_h5ad_canned_")
    anndata = create_anndata_canned(which, obs_field_name, var_field_name)
    return create_h5ad(
        anndata,
        (tmp_path.name + f"{which}.h5ad"),
    )


def create_soma_canned(which: int, obs_field_name, var_field_name):
    tmp_path = tempfile.TemporaryDirectory(prefix="create_soma_canned_")
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
        obs_field_name="cell_id",
        var_field_name="gene_id",
    )


@pytest.fixture
def soma_larger(anndata_larger):
    tmp_path = tempfile.TemporaryDirectory(prefix="soma_larger_")
    uri = tmp_path.name + "soma-larger"
    tiledbsoma.io.from_anndata(uri, anndata_larger, "measname")
    return uri


# fmt: off
PANDAS_INDEXING_TEST_DF = pd.DataFrame(
    data={
        "soma_joinid": np.arange(3, dtype=np.int64),
        "alt_id": ["A", "C", "G"],
        "obs_id": ["AT", "CT", "GT"],
    }
)
@pytest.mark.parametrize(
    [          "index_col_and_name"      ,  "default_index_name"  ,  "signature_col_names"  ],
    [   # |   Set this   |  If present,  |     signatures.py      |        Expected:        |
        # |  col as idx  | rename index  |  `default_index_name`  |   signature col names   |

        # `default_index_name` matches column that was made index ⇒ both columns present in signature
        [ (   "obs_id"   ,               ),       "obs_id"        , ( "obs_id" , "alt_id" ) ],
        [ (   "obs_id"   ,    "index"    ),       "obs_id"        , ( "obs_id" , "alt_id" ) ],
        [ (   "obs_id"   ,      None     ),       "obs_id"        , ( "obs_id" , "alt_id" ) ],
        [ (   "alt_id"   ,               ),       "alt_id"        , ( "alt_id" , "obs_id" ) ],
        [ (   "alt_id"   ,    "index"    ),       "alt_id"        , ( "alt_id" , "obs_id" ) ],
        [ (   "alt_id"   ,      None     ),       "alt_id"        , ( "alt_id" , "obs_id" ) ],

        # `default_index_name` is the column that was not made index ⇒ index dropped when named "index" or None
        [ (   "alt_id"   ,               ),       "obs_id"        , ( "obs_id" , "alt_id" ) ],
        [ (   "alt_id"   ,    "index"    ),       "obs_id"        , ( "obs_id" ,          ) ],
        [ (   "alt_id"   ,      None     ),       "obs_id"        , ( "obs_id" ,          ) ],
        [ (   "obs_id"   ,               ),       "alt_id"        , ( "alt_id" , "obs_id" ) ],
        [ (   "obs_id"   ,    "index"    ),       "alt_id"        , ( "alt_id" ,          ) ],
        [ (   "obs_id"   ,      None     ),       "alt_id"        , ( "alt_id" ,          ) ],

        # default RangeIndex ⇒ columns are preserved
        [ (      None    ,               ),       "obs_id"        , ( "obs_id" , "alt_id" ) ],
        [ (      None    ,               ),       "alt_id"        , ( "alt_id" , "obs_id" ) ],
    ]
)
def test_pandas_indexing(
    index_col_and_name: Union[Tuple[str | None], Tuple[str, str | None]],
    default_index_name: str,
    signature_col_names: List[Union[str, Tuple[str, str]]],
):
    """
    The `default_index_name` for registration can interact with column- and
    index-names in a variety of ways; this test exercises several of them.
    """
    df = PANDAS_INDEXING_TEST_DF.copy()
    index_col = index_col_and_name[0]
    if index_col is not None:
        df.set_index(index_col, inplace=True)
        if len(index_col_and_name) == 2:
            df.index.name = index_col_and_name[1]

    arrow_schema = conversions.df_to_arrow_schema(df, default_index_name)
    actual_signature = conversions._string_dict_from_arrow_schema(arrow_schema)
    expected_signature = dict((col, "string") for col in signature_col_names)
    assert actual_signature == expected_signature
# fmt: on


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

    for reload in [False, True]:
        if reload:
            dictionary = registration.AxisAmbientLabelMapping.from_json(
                dictionary.to_json()
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

    for reload in [False, True]:
        if reload:
            rd = registration.ExperimentAmbientLabelMapping.from_json(rd.to_json())

    assert rd.obs_axis.id_mapping_from_values([]).data == ()
    assert rd.obs_axis.id_mapping_from_values(["AGAG", "ACTG"]).data == (2, 1)
    assert rd.var_axes["measname"].id_mapping_from_values(["TP53", "VEGFA"]).data == (
        3,
        4,
    )
    assert rd.var_axes["raw"].id_mapping_from_values(
        ["RAW2", "TP53", "VEGFA"]
    ).data == (6, 3, 4)

    assert rd.get_obs_shape() == 3
    assert rd.get_var_shapes() == {"measname": 5, "raw": 7}


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

    assert rd.get_obs_shape() == 3
    assert rd.get_var_shapes() == {"measname": 5, "raw": 7}


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

    assert rd.get_obs_shape() == 3
    assert rd.get_var_shapes() == {"measname": 5, "raw": 7}


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

        nobs = rd.get_obs_shape()
        nvars = rd.get_var_shapes()
        tiledbsoma.io.resize_experiment(experiment_uri, nobs=nobs, nvars=nvars)

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

    assert rd.get_obs_shape() == 12
    assert rd.get_var_shapes() == {"measname": 7, "raw": 10}

    # Now do the ingestion per se.  Note that once registration is done sequentially, ingest order
    # mustn't matter, and in fact, can be done in parallel. This is why we test various permutations
    # of the ordering of the h5ad file names.
    for h5ad_file_name in [
        h5ad_file_names[permutation[0]],
        h5ad_file_names[permutation[1]],
        h5ad_file_names[permutation[2]],
        h5ad_file_names[permutation[3]],
    ]:

        if tiledbsoma.Experiment.exists(experiment_uri):
            tiledbsoma.io.resize_experiment(
                experiment_uri,
                nobs=rd.get_obs_shape(),
                nvars=rd.get_var_shapes(),
            )

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

    assert rd.get_obs_shape() == 12
    assert rd.get_var_shapes() == {"measname": 7, "raw": 10}


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

    assert rd.get_obs_shape() == 6
    assert rd.get_var_shapes() == {"measname": 5, "raw": 7}

    adata2 = ad.read_h5ad(h5ad2)

    original = adata2.copy()

    tiledbsoma.io.resize_experiment(
        soma1,
        nobs=rd.get_obs_shape(),
        nvars=rd.get_var_shapes(),
    )

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

    assert_adata_equal(original, adata2)

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

    tiledbsoma.io.resize_experiment(
        soma_uri,
        nobs=rd.get_obs_shape(),
        nvars=rd.get_var_shapes(),
    )

    tiledbsoma.io.from_anndata(
        soma_uri,
        anndata2,
        measurement_name="two",
        registration_mapping=rd,
    )

    assert_adata_equal(original, anndata2)

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

    assert rd.get_obs_shape() == 1000
    assert rd.get_var_shapes() == {"measname": 6}


def test_ealm_expose():
    """Checks that this is exported from tiledbsoma.io._registration"""
    # All we want to check is that the import doesn't throw. Job done. Period.
    # However, the pre-commit hook will strip out this import statement as "unused".
    # So, assert something.
    assert tiledbsoma.io.ExperimentAmbientLabelMapping is not None


def test_append_registration_with_nonexistent_storage(tmp_path):
    anndata1 = create_anndata_canned(1, "obs_id", "var_id")
    anndata2 = create_anndata_canned(2, "obs_id", "var_id")
    soma_uri = tmp_path.as_posix()

    tiledbsoma.io.from_anndata(soma_uri, anndata1, measurement_name="RNA")

    with pytest.raises(ValueError):
        tiledbsoma.io.register_anndatas(
            soma_uri + "-nonesuch",
            [anndata2],
            measurement_name="RNA",
            obs_field_name="obs_id",
            var_field_name="var_id",
        )


@pytest.mark.parametrize("obs_field_name", ["obs_id", "cell_id"])
@pytest.mark.parametrize("var_field_name", ["var_id", "gene_id"])
@pytest.mark.parametrize(
    "dataset_ids_and_exc",
    [
        [None, 2],
        [ValueError, 8],
        [ValueError, 9],
    ],
)
def test_append_with_nonunique_field_values(
    tmp_path,
    obs_field_name,
    var_field_name,
    dataset_ids_and_exc,
):
    """Verifies that we do a proactive check for uniqueness of obs/var registration-field values"""
    ida = 1
    exc, idb = dataset_ids_and_exc
    measurement_name = "test"

    anndataa = create_anndata_canned(ida, obs_field_name, var_field_name)
    anndatab = create_anndata_canned(idb, obs_field_name, var_field_name)
    soma_uri = tmp_path.as_posix()

    tiledbsoma.io.from_anndata(soma_uri, anndataa, measurement_name=measurement_name)

    ctx = pytest.raises(exc) if exc else nullcontext()
    with ctx:
        tiledbsoma.io.register_anndatas(
            soma_uri,
            [anndatab],
            measurement_name=measurement_name,
            obs_field_name=obs_field_name,
            var_field_name=var_field_name,
        )


@pytest.mark.parametrize("all_at_once", [False, True])
@pytest.mark.parametrize("nobs_a", [50, 300])
@pytest.mark.parametrize("nobs_b", [60, 400])
def test_enum_bit_width_append(tmp_path, all_at_once, nobs_a, nobs_b):
    """Creates an obs column whose bit width might naively be inferred to be int8
    by tiledbsoma.io, and another which could be inferred to int16.  Then
    ensures the dataframes are appendable regardless of which one was written
    first."""
    obs_ids_a = [("a_%08d" % e) for e in range(nobs_a)]
    obs_ids_b = [("b_%08d" % e) for e in range(nobs_b)]
    var_ids = ["W", "X", "Y", "Z"]
    obs_field_name = "cell_id"
    var_field_name = "gene_id"
    measurement_name = "meas"

    adata = _create_anndata(
        obs_ids=obs_ids_a,
        var_ids=var_ids,
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
        X_value_base=0,
    )

    bdata = _create_anndata(
        obs_ids=obs_ids_b,
        var_ids=var_ids,
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
        X_value_base=100,
    )

    adata.obs["enum"] = pd.Categorical(obs_ids_a, categories=obs_ids_a)
    bdata.obs["enum"] = pd.Categorical(obs_ids_b, categories=obs_ids_b)

    soma_uri = tmp_path.as_posix()

    if all_at_once:
        rd = tiledbsoma.io.register_anndatas(
            None,
            [adata, bdata],
            measurement_name=measurement_name,
            obs_field_name=obs_field_name,
            var_field_name=var_field_name,
        )

        assert rd.get_obs_shape() == nobs_a + nobs_b
        assert rd.get_var_shapes() == {"meas": 4, "raw": 0}

        tiledbsoma.io.from_anndata(
            soma_uri, adata, measurement_name=measurement_name, registration_mapping=rd
        )

        tiledbsoma.io.resize_experiment(
            soma_uri,
            nobs=rd.get_obs_shape(),
            nvars=rd.get_var_shapes(),
        )

        tiledbsoma.io.from_anndata(
            soma_uri, bdata, measurement_name=measurement_name, registration_mapping=rd
        )

    else:
        tiledbsoma.io.from_anndata(soma_uri, adata, measurement_name=measurement_name)

        rd = tiledbsoma.io.register_anndatas(
            soma_uri,
            [bdata],
            measurement_name=measurement_name,
            obs_field_name=obs_field_name,
            var_field_name=var_field_name,
        )

        assert rd.get_obs_shape() == nobs_a + nobs_b
        assert rd.get_var_shapes() == {"meas": 4}

        tiledbsoma.io.resize_experiment(
            soma_uri,
            nobs=rd.get_obs_shape(),
            nvars=rd.get_var_shapes(),
        )

        tiledbsoma.io.from_anndata(
            soma_uri, bdata, measurement_name=measurement_name, registration_mapping=rd
        )

    with tiledbsoma.Experiment.open(soma_uri) as exp:
        obs = exp.obs.read().concat()

        cell_ids = obs[obs_field_name].to_pylist()

        readback_a = cell_ids[:nobs_a]
        readback_b = cell_ids[nobs_a:]

        assert readback_a == obs_ids_a
        assert readback_b == obs_ids_b


def test_multimodal_names(tmp_path, conftest_pbmc3k_adata):
    uri = tmp_path.as_posix()

    # Data for "RNA" measurement of SOMA experiment
    adata_rna = conftest_pbmc3k_adata.copy()
    adata_rna.obs.index.name = "cell_id"
    adata_rna.var.index.name = "first_adata_var_index"

    # Non-appendable
    del adata_rna.obsm
    del adata_rna.varm
    del adata_rna.obsp
    del adata_rna.varp

    # Simulate for "protein" measurement of SOMA experiment:
    # * Different var values
    # * Different number of var rows
    # * Different var field name for registration

    adata_protein = adata_rna.copy()
    adata_protein.var = pd.DataFrame(
        {"assay_type": ["protein"] * adata_protein.n_vars},
        index=[f"p{i}" for i in range(adata_protein.n_vars)],
    )
    adata_protein = adata_protein[:, :500]

    adata_protein.obs["batch_id"] = np.nan
    adata_protein.obs["batch_id"] = adata_protein.obs["batch_id"].astype("string")
    assert adata_protein.obs["batch_id"].dtype == pd.StringDtype()

    adata_protein.var.index.name = "second_adata_var_index"

    tiledbsoma.io.from_anndata(
        experiment_uri=uri,
        anndata=adata_rna,
        measurement_name="RNA",
        uns_keys=[],
    )

    with tiledbsoma.Experiment.open(uri) as exp:
        assert "RNA" in exp.ms
        assert "protein" not in exp.ms

    with tiledbsoma.open(uri, "w") as exp:
        exp.ms.add_new_collection(
            "protein",
            kind=tiledbsoma.Measurement,
            uri="protein",  # relative path for local-disk operations
        )

    # Register the second anndata object in the protein measurement
    rd = tiledbsoma.io.register_anndatas(
        experiment_uri=uri,
        adatas=[adata_protein],
        measurement_name="protein",
        obs_field_name=adata_protein.obs.index.name,
        var_field_name=adata_protein.var.index.name,
    )

    assert rd.get_obs_shape() == 2638
    assert rd.get_var_shapes() == {"protein": 500, "raw": 13714}

    # Ingest the second anndata object into the protein measurement
    tiledbsoma.io.from_anndata(
        experiment_uri=uri,
        anndata=adata_protein,
        measurement_name="protein",
        registration_mapping=rd,
        uns_keys=[],
    )

    with tiledbsoma.Experiment.open(uri) as exp:
        assert "RNA" in exp.ms
        assert "protein" in exp.ms

        assert exp.obs.count == len(adata_rna.obs)
        assert exp.obs.count == len(adata_protein.obs)
        assert exp.ms["RNA"].var.count == len(adata_rna.var)
        assert exp.ms["protein"].var.count == len(adata_protein.var)


def test_registration_lists_and_tuples(tmp_path):
    obs_field_name = "cell_id"
    var_field_name = "gene_id"

    exp_uri = create_soma_canned(1, obs_field_name, var_field_name)
    adata = create_anndata_canned(2, obs_field_name, var_field_name)
    h5ad_file_name = create_h5ad_canned(2, obs_field_name, var_field_name)

    rd1 = tiledbsoma.io.register_anndatas(
        experiment_uri=exp_uri,
        adatas=[adata],
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )

    rd2 = tiledbsoma.io.register_anndatas(
        experiment_uri=exp_uri,
        adatas=(adata,),
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )

    rd3 = tiledbsoma.io.register_anndatas(
        experiment_uri=exp_uri,
        adatas=adata,
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )
    assert rd1 == rd2
    assert rd2 == rd3

    rd4 = tiledbsoma.io.register_h5ads(
        experiment_uri=exp_uri,
        h5ad_file_names=[h5ad_file_name],
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )

    rd5 = tiledbsoma.io.register_h5ads(
        experiment_uri=exp_uri,
        h5ad_file_names=(h5ad_file_name,),
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )

    rd6 = tiledbsoma.io.register_h5ads(
        experiment_uri=exp_uri,
        h5ad_file_names=h5ad_file_name,
        measurement_name="measname",
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
    )

    assert rd4 == rd5
    assert rd5 == rd6


@pytest.mark.parametrize(
    "version_and_shaped",
    [
        ["1.7.3", False],
        ["1.12.3", False],
        ["1.14.5", False],
        ["1.15.0", True],
        ["1.15.7", True],
    ],
)
def test_extend_enmr_to_older_experiments_64521(tmp_path, version_and_shaped):
    version, shaped = version_and_shaped

    import os
    import shutil

    from ._util import ROOT_DATA_DIR

    original_data_uri = str(
        ROOT_DATA_DIR / "soma-experiment-versions" / version / "pbmc3k_unprocessed"
    )

    if not os.path.isdir(original_data_uri):
        raise RuntimeError(
            f"Missing '{original_data_uri}' directory. Try running `make data` "
            "from the TileDB-SOMA project root directory."
        )

    with tiledbsoma.Experiment.open(original_data_uri) as exp:
        assert exp.obs.count == 2700

    # Make a copy of the Experiment as to not write over the data in ROOT_DATA_DIR
    uri = (tmp_path / version).as_posix()
    shutil.copytree(original_data_uri, uri)

    with tiledbsoma.Experiment.open(uri) as exp:
        assert exp.obs.count == 2700
        obs = exp.obs.read().concat().to_pandas()
        assert "new_ident" not in obs["orig.ident"].cat.categories

    with tiledbsoma.Experiment.open(uri) as exp:
        adata = tiledbsoma.io.to_anndata(exp, "RNA")

    # Make obs_id accessible via adata["obs_id]]
    adata.obs.reset_index(inplace=True)

    adata.obs["orig.ident"] = pd.Series(
        ["new_ident"] * len(adata.obs), dtype="category"
    )
    adata.obs["obs_id"] = adata.obs["obs_id"] + "_2"

    rd = tiledbsoma.io.register_anndatas(
        experiment_uri=uri,
        adatas=[adata],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id",
    )

    assert rd.get_obs_shape() == 5400

    if shaped:
        tiledbsoma.io.resize_experiment(
            uri,
            nobs=rd.get_obs_shape(),
            nvars=rd.get_var_shapes(),
        )

    tiledbsoma.io.from_anndata(
        experiment_uri=uri,
        anndata=adata,
        measurement_name="RNA",
        registration_mapping=rd,
    )

    with tiledbsoma.Experiment.open(uri) as exp:
        assert "RNA" in exp.ms

        assert exp.obs.count == 5400
        obs = exp.obs.read().concat().to_pandas()
        assert "pbmc3k" in obs["orig.ident"].cat.categories
        assert "new_ident" in obs["orig.ident"].cat.categories
