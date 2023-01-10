import sys
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from scipy import sparse

import tiledbsoma as soma
from tiledbsoma.experiment_query import AxisQuery, ExperimentAxisQuery, X_as_series

"""
WIP tracker - delete when complete.  Missing tests:
* query by coords when there is >1 dimension
* error check handling of Dense X/obsm/varm
"""


@pytest.fixture(scope="function")
def obs(tmp_path, n_obs) -> soma.DataFrame:
    return make_dataframe((tmp_path / "obs").as_posix(), n_obs)


@pytest.fixture(scope="function")
def var(tmp_path, n_vars) -> soma.DataFrame:
    return make_dataframe((tmp_path / "var").as_posix(), n_vars)


@pytest.fixture
def X_layer_names():
    return ["raw"]


@pytest.fixture
def obsp_layer_names():
    return None


@pytest.fixture
def varp_layer_names():
    return None


@pytest.fixture(scope="function")
def soma_experiment(
    tmp_path,
    n_obs,
    n_vars,
    obs,
    var,
    X_layer_names,
    obsp_layer_names,
    varp_layer_names,
):
    return make_experiment(
        tmp_path,
        n_obs,
        n_vars,
        obs,
        var,
        X_layer_names,
        obsp_layer_names,
        varp_layer_names,
    )


@pytest.mark.xfail(
    # This test fails on Python 3.10+ due to a bug in typeguard. The bug
    # is tripped due to use of a TypedDict. Remove
    # work-around when the typeguard issue is fixed AND released.
    # Underlying issue:
    #   https://github.com/agronholm/typeguard/issues/242
    # Related to _when_ we might see a release:
    #   https://github.com/agronholm/typeguard/issues/257
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars,X_layer_names", [(101, 11, ("raw", "extra"))])
def test_experiment_query_all(soma_experiment):
    """Test a query with default obs_query / var_query - i.e., query all."""
    assert soma_experiment.exists()

    with ExperimentAxisQuery(soma_experiment, "RNA") as query:
        assert query.n_obs == 101
        assert query.n_vars == 11

        assert np.array_equal(query.obs_joinids().to_numpy(), np.arange(101))
        assert np.array_equal(query.var_joinids().to_numpy(), np.arange(11))

        assert len(query.obs().concat()) == 101
        assert len(query.var().concat()) == 11

        assert query.obs().concat().to_pydict() == {
            "soma_joinid": list(range(101)),
            "label": [str(i) for i in range(101)],
        }
        assert query.var().concat().to_pydict() == {
            "soma_joinid": list(range(11)),
            "label": [str(i) for i in range(11)],
        }
        assert pa.concat_tables(query.X("raw").tables()) == pa.concat_tables(
            soma_experiment.ms["RNA"].X["raw"].read((slice(None), slice(None))).tables()
        )
        assert query.X("raw").tables().concat() == pa.concat_tables(
            query.X("raw").tables()
        )

        # read as table
        arrow_reads = query.read("raw")
        assert "X_layers" not in arrow_reads
        assert isinstance(arrow_reads["obs"], pa.Table)
        assert isinstance(arrow_reads["var"], pa.Table)
        assert isinstance(arrow_reads["X"], pa.Table)
        assert arrow_reads["obs"] == query.obs().concat()
        assert arrow_reads["var"] == query.var().concat()
        assert (
            arrow_reads["X"]
            == soma_experiment.ms["RNA"]
            .X["raw"]
            .read((slice(None), slice(None)))
            .tables()
            .concat()
        )

        # read as anndata
        ad = query.to_anndata("raw")
        assert ad.n_obs == query.n_obs and ad.n_vars == query.n_vars
        assert set(ad.obs.keys().to_list()) == set(["soma_joinid", "label"])
        assert set(ad.var.keys().to_list()) == set(["soma_joinid", "label"])

        raw_X = (
            soma_experiment.ms["RNA"]
            .X["raw"]
            .read((slice(None), slice(None)))
            .tables()
            .concat()
        )
        ad_X_coo = ad.X.tocoo()
        assert np.array_equal(raw_X["soma_dim_0"], ad_X_coo.row)
        assert np.array_equal(raw_X["soma_dim_1"], ad_X_coo.col)
        assert np.array_equal(raw_X["soma_data"], ad_X_coo.data)


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_coords(soma_experiment):
    """Test query by dimension coordinates"""
    obs_slice = slice(3, 72)
    var_slice = slice(7, 21)
    with soma_experiment.axis_query(
        "RNA",
        obs_query=AxisQuery(coords=(obs_slice,)),
        var_query=AxisQuery(coords=(var_slice,)),
    ) as query:
        assert query.n_obs == obs_slice.stop - obs_slice.start + 1
        assert query.n_vars == var_slice.stop - var_slice.start + 1
        assert np.array_equal(
            query.obs_joinids().to_numpy(),
            np.arange(obs_slice.start, obs_slice.stop + 1),
        )
        assert np.array_equal(
            query.var_joinids().to_numpy(),
            np.arange(var_slice.start, var_slice.stop + 1),
        )

        assert np.array_equal(
            query.obs(column_names=["soma_joinid"]).concat()["soma_joinid"].to_numpy(),
            np.arange(obs_slice.start, obs_slice.stop + 1),
        )
        assert np.array_equal(
            query.var(column_names=["soma_joinid"]).concat()["soma_joinid"].to_numpy(),
            np.arange(var_slice.start, var_slice.stop + 1),
        )

        raw_X = (
            soma_experiment.ms["RNA"]
            .X["raw"]
            .read((obs_slice, var_slice))
            .tables()
            .concat()
        )
        assert query.X("raw").tables().concat() == raw_X
        assert query.X("raw").coos().concat() == pa.SparseCOOTensor.from_numpy(
            raw_X["soma_data"].to_numpy(),
            np.array(
                [
                    raw_X["soma_dim_0"].to_numpy(),
                    raw_X["soma_dim_1"].to_numpy(),
                ]
            ).T,
            shape=soma_experiment.ms["RNA"].X["raw"].shape,
        )


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_value_filter(soma_experiment):
    """Test query by value filter"""
    obs_label_values = ["3", "7", "38", "99"]
    var_label_values = ["18", "34", "67"]
    with soma_experiment.axis_query(
        "RNA",
        obs_query=AxisQuery(value_filter=f"label in {obs_label_values}"),
        var_query=AxisQuery(value_filter=f"label in {var_label_values}"),
    ) as query:
        assert query.n_obs == len(obs_label_values)
        assert query.n_vars == len(var_label_values)
        assert query.obs().concat()["label"].to_pylist() == obs_label_values
        assert query.var().concat()["label"].to_pylist() == var_label_values


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_combo(soma_experiment):
    """Test query by combinations of coords and value_filter"""
    obs_label_values = ["3", "7", "38", "99"]
    var_label_values = ["18", "34", "67"]
    obs_slice = slice(3, 101)
    var_slice = slice(7, 80)

    with soma_experiment.axis_query(
        "RNA",
        obs_query=AxisQuery(coords=(obs_slice,)),
        var_query=AxisQuery(value_filter=f"label in {var_label_values}"),
    ) as query:
        assert query.n_obs == obs_slice.stop - obs_slice.start + 1
        assert query.var().concat()["label"].to_pylist() == var_label_values

    with soma_experiment.axis_query(
        "RNA",
        obs_query=AxisQuery(value_filter=f"label in {obs_label_values}"),
        var_query=AxisQuery(coords=(var_slice,)),
    ) as query:
        assert query.obs().concat()["label"].to_pylist() == obs_label_values
        assert np.array_equal(
            query.var_joinids().to_numpy(),
            np.arange(var_slice.start, var_slice.stop + 1),
        )

    with soma_experiment.axis_query(
        "RNA",
        obs_query=AxisQuery(
            coords=(obs_slice,), value_filter=f"label in {obs_label_values}"
        ),
        var_query=AxisQuery(
            coords=(var_slice,), value_filter=f"label in {var_label_values}"
        ),
    ) as query:
        assert query.obs().concat()["label"].to_pylist() == obs_label_values
        assert query.var().concat()["label"].to_pylist() == var_label_values


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars,X_layer_names", [(1001, 99, ["A", "B", "C"])])
def test_X_layers(soma_experiment):
    """Verify multi-layer-X handling"""
    assert soma_experiment.exists()
    A = pa.concat_tables(
        soma_experiment.ms["RNA"].X["A"].read((slice(None), slice(None))).tables()
    )
    B = pa.concat_tables(
        soma_experiment.ms["RNA"].X["B"].read((slice(None), slice(None))).tables()
    )

    with soma_experiment.axis_query("RNA") as query:
        arrow_reads = query.read("A", X_layers=["B"])
        assert arrow_reads["X"] == A
        assert arrow_reads["X_layers"]["B"] == B

        ad = query.to_anndata("B", X_layers=["A"])
        ad_X_coo = ad.X.tocoo()
        assert np.array_equal(B["soma_dim_0"], ad_X_coo.row)
        assert np.array_equal(B["soma_dim_1"], ad_X_coo.col)
        assert np.array_equal(B["soma_data"], ad_X_coo.data)

        ad_lyr_coo = ad.layers["A"].tocoo()
        assert np.array_equal(A["soma_dim_0"], ad_lyr_coo.row)
        assert np.array_equal(A["soma_dim_1"], ad_lyr_coo.col)
        assert np.array_equal(A["soma_data"], ad_lyr_coo.data)


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_indexer(soma_experiment):
    """Test result indexer"""
    assert soma_experiment.exists()

    with ExperimentAxisQuery(
        soma_experiment,
        "RNA",
        obs_query=AxisQuery(coords=(slice(1, 10),)),
        var_query=AxisQuery(coords=(slice(1, 10),)),
    ) as query:
        indexer = query.get_indexer()

        # coords outside of our query should return -1
        assert np.array_equal(
            indexer.obs_index(np.array([-1, 0, 11, 1003])),
            np.array([-1, -1, -1, -1]),
        )
        assert np.array_equal(
            indexer.var_index(np.array([-1, 0, 11, 1003])),
            np.array([-1, -1, -1, -1]),
        )

        # inside results, indexed
        assert np.array_equal(
            indexer.obs_index(np.array([1, 4, 2])), np.array([0, 3, 1])
        )
        assert np.array_equal(indexer.var_index(np.array([10, 1])), np.array([9, 0]))


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(2833, 107)])
def test_error_corners(soma_experiment: soma.Experiment):
    """Verify a couple of error conditions / corner cases."""
    assert soma_experiment.exists()

    with pytest.raises(ValueError):
        soma_experiment.axis_query("no-such-measurement")

    with pytest.raises(ValueError):
        soma.Experiment(uri="foobar").axis_query("foobar")

    with pytest.raises(ValueError):
        with soma_experiment.axis_query("RNA") as query:
            next(query.X("no-such-layer"))

    with pytest.raises(ValueError):
        with soma_experiment.axis_query("RNA") as query:
            next(query.obsp("no-such-layer"))

    with pytest.raises(ValueError):
        with soma_experiment.axis_query("RNA") as query:
            next(query.varp("no-such-layer"))


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_query_cleanup(soma_experiment: soma.Experiment):
    """
    Verify soma.Experiment.query works as context manager and stand-alone,
    and that it cleans up correct.
    """
    from contextlib import closing

    with soma_experiment.axis_query("RNA") as query:
        assert query.n_obs == 1001
        assert query.n_vars == 99
        assert query.to_anndata("raw") is not None
        assert query._threadpool_ is not None

    assert query._threadpool_ is None

    with closing(soma_experiment.axis_query("RNA")) as query:
        assert query.to_anndata("raw") is not None
        assert query._threadpool_ is not None

    assert query._threadpool_ is None


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize(
    "n_obs,n_vars,obsp_layer_names,varp_layer_names", [(1001, 99, ["foo"], ["bar"])]
)
def test_experiment_query_obsp_varp(soma_experiment):
    obs_slice = slice(3, 72)
    var_slice = slice(7, 21)
    with ExperimentAxisQuery(
        soma_experiment,
        "RNA",
        obs_query=AxisQuery(coords=(obs_slice,)),
        var_query=AxisQuery(coords=(var_slice,)),
    ) as query:
        assert query.n_obs == obs_slice.stop - obs_slice.start + 1
        assert query.n_vars == var_slice.stop - var_slice.start + 1

        with pytest.raises(ValueError):
            next(query.obsp("no-such-layer"))

        with pytest.raises(ValueError):
            next(query.varp("no-such-layer"))

        assert (
            query.obsp("foo").tables().concat()
            == soma_experiment.ms["RNA"]
            .obsp["foo"]
            .read((obs_slice, obs_slice))
            .tables()
            .concat()
        )

        assert (
            query.varp("bar").tables().concat()
            == soma_experiment.ms["RNA"]
            .varp["bar"]
            .read((var_slice, var_slice))
            .tables()
            .concat()
        )


def test_axis_query():
    """Basic test of the AxisQuery class"""
    assert AxisQuery().coords == (slice(None),)
    assert AxisQuery().value_filter is None
    assert AxisQuery() == AxisQuery(coords=(slice(None),))

    assert AxisQuery(coords=(1,)).coords == (1,)
    assert AxisQuery(coords=(slice(1, 2),)).coords == (slice(1, 2),)
    assert AxisQuery(coords=((1, 88),)).coords == ((1, 88),)

    assert AxisQuery(coords=(1, 2)).coords == (1, 2)
    assert AxisQuery(coords=(slice(1, 2), slice(None))).coords == (
        slice(1, 2),
        slice(None),
    )
    assert AxisQuery(coords=(slice(1, 2),)).value_filter is None

    assert AxisQuery(value_filter="foo == 'bar'").value_filter == "foo == 'bar'"
    assert AxisQuery(value_filter="foo == 'bar'").coords == (slice(None),)

    assert AxisQuery(coords=(slice(1, 100),), value_filter="foo == 'bar'").coords == (
        slice(1, 100),
    )
    assert (
        AxisQuery(coords=(slice(1, 100),), value_filter="foo == 'bar'").value_filter
        == "foo == 'bar'"
    )

    with pytest.raises(TypeError):
        AxisQuery(coords=True)

    with pytest.raises(TypeError):
        AxisQuery(value_filter=[])

    with pytest.raises(TypeError):
        AxisQuery(coords=({},))


def test_X_as_series():
    soma_dim_0 = np.arange(0, 100, dtype=np.int64)
    soma_dim_1 = np.arange(200, 300, dtype=np.int64)
    soma_data = np.random.default_rng().standard_normal(100, dtype=np.float32)
    ser = X_as_series(
        pa.Table.from_arrays(
            [soma_dim_0, soma_dim_1, soma_data],
            names=["soma_dim_0", "soma_dim_1", "soma_data"],
        )
    )

    assert isinstance(ser, pd.Series)
    assert np.array_equal(ser.to_numpy(), soma_data)
    assert np.array_equal(
        ser.index.get_level_values("soma_dim_0").to_numpy(), soma_dim_0
    )
    assert np.array_equal(
        ser.index.get_level_values("soma_dim_1").to_numpy(), soma_dim_1
    )


"""
Fixture support & utility functions below.
"""


def make_dataframe(path: str, sz: int) -> soma.DataFrame:
    df = soma.DataFrame(uri=path)
    df.create(
        schema=pa.schema(
            [
                ("soma_joinid", pa.int64()),
                ("label", pa.large_string()),
            ]
        ),
        index_column_names=["soma_joinid"],
    )
    df.write(
        pa.Table.from_pydict(
            {
                "soma_joinid": [i for i in range(sz)],
                "label": [str(i) for i in range(sz)],
            }
        )
    )
    return df


def make_sparse_array(path: str, shape: Tuple[int, int]) -> soma.SparseNDArray:
    a = soma.SparseNDArray(path).create(pa.float32(), shape)
    tensor = pa.SparseCOOTensor.from_scipy(
        sparse.random(
            shape[0],
            shape[1],
            density=0.1,
            format="coo",
            dtype=np.float32,
            random_state=np.random.default_rng(),
        )
    )
    a.write(tensor)
    return a


def make_experiment(
    root: str,
    n_obs: int,
    n_vars: int,
    obs: soma.DataFrame,
    var: soma.DataFrame,
    X_layer_names: List[str],  # will create a random matrix per layer name
    obsp_layer_names: Optional[List[str]] = None,
    varp_layer_names: Optional[List[str]] = None,
) -> soma.Experiment:

    assert obs.exists()
    assert var.exists()
    assert len(obs) == n_obs
    assert len(var) == n_vars

    experiment = soma.Experiment((root / "exp").as_posix()).create()
    experiment["ms"] = soma.Collection((root / "ms").as_posix()).create()
    experiment.ms["RNA"] = soma.Measurement(uri=f"{experiment.ms.uri}/RNA").create()
    experiment["obs"] = obs

    measurement = experiment.ms["RNA"]
    measurement["var"] = var
    measurement["X"] = soma.Collection(uri=f"{measurement.uri}/X").create()

    (root / "X").mkdir()
    for X_layer_name in X_layer_names:
        path = root / "X" / X_layer_name
        path.mkdir()
        measurement.X[X_layer_name] = make_sparse_array(
            path.as_posix(), (n_obs, n_vars)
        )

    if obsp_layer_names:
        (root / "obsp").mkdir()
        measurement["obsp"] = soma.Collection(uri=f"{measurement.uri}/obsp").create()
        for obsp_layer_name in obsp_layer_names:
            path = root / "obsp" / obsp_layer_name
            path.mkdir()
            measurement.obsp[obsp_layer_name] = make_sparse_array(
                path.as_posix(), (n_obs, n_obs)
            )

    if varp_layer_names:
        (root / "varp").mkdir()
        measurement["varp"] = soma.Collection(uri=f"{measurement.uri}/varp").create()
        for varp_layer_name in varp_layer_names:
            path = root / "varp" / varp_layer_name
            path.mkdir()
            measurement.varp[varp_layer_name] = make_sparse_array(
                path.as_posix(), (n_vars, n_vars)
            )

    return experiment
