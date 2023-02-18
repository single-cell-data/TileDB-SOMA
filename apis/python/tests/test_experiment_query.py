import sys
from typing import Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from scipy import sparse

import tiledbsoma as soma
from tiledbsoma import _factory
from tiledbsoma._collection import CollectionBase
from tiledbsoma.experiment_query import X_as_series


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
    X_layer_names,
    obsp_layer_names,
    varp_layer_names,
):
    with soma.Experiment.create((tmp_path / "exp").as_posix()) as exp:
        add_dataframe(exp, "obs", n_obs)
        ms = exp.add_new_collection("ms")
        rna = ms.add_new_collection("RNA", soma.Measurement)
        add_dataframe(rna, "var", n_vars)
        rna_x = rna.add_new_collection("X", soma.Collection)
        for X_layer_name in X_layer_names:
            add_sparse_array(rna_x, X_layer_name, (n_obs, n_vars))

        if obsp_layer_names:
            obsp = rna.add_new_collection("obsp")
            for obsp_layer_name in obsp_layer_names:
                add_sparse_array(obsp, obsp_layer_name, (n_obs, n_obs))

        if varp_layer_names:
            varp = rna.add_new_collection("varp")
            for varp_layer_name in varp_layer_names:
                add_sparse_array(varp, varp_layer_name, (n_vars, n_vars))
    return _factory.open((tmp_path / "exp").as_posix())


@pytest.mark.xfail(
    # This test fails on Python 3.10+ due to a bug in typeguard. The bug
    # is tripped due to use of a TypedDict. Remove
    # work-around when the typeguard issue is fixed AND released.
    # Underlying issue:
    #   https://github.com/agronholm/typeguard/issues/242
    # Related to _when_ we might see a release:
    #   https://github.com/agronholm/typeguard/issues/257
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars,X_layer_names", [(101, 11, ("raw", "extra"))])
def test_experiment_query_all(soma_experiment):
    """Test a query with default obs_query / var_query - i.e., query all."""
    with soma.ExperimentAxisQuery(soma_experiment, "RNA") as query:
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

        # read as anndata
        ad = query.to_anndata("raw")
        assert ad.n_obs == query.n_obs and ad.n_vars == query.n_vars

        assert set(ad.obs.keys().to_list()) == set(["soma_joinid", "label"])
        assert set(ad.var.keys().to_list()) == set(["soma_joinid", "label"])

        obs = soma_experiment.obs.read().concat().to_pandas()
        obs.index = obs.index.map(str)
        assert (obs == ad.obs).all().all()

        var = soma_experiment.ms["RNA"].var.read().concat().to_pandas()
        var.index = var.index.map(str)
        assert (var == ad.var).all().all()

        assert len(ad.layers) == 0

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
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_coords(soma_experiment):
    """Test query by dimension coordinates"""
    obs_slice = slice(3, 72)
    var_slice = slice(7, 21)
    with soma_experiment.axis_query(
        "RNA",
        obs_query=soma.AxisQuery(coords=(obs_slice,)),
        var_query=soma.AxisQuery(coords=(var_slice,)),
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
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_value_filter(soma_experiment):
    """Test query by value filter"""
    obs_label_values = ["3", "7", "38", "99"]
    var_label_values = ["18", "34", "67"]
    with soma_experiment.axis_query(
        "RNA",
        obs_query=soma.AxisQuery(value_filter=f"label in {obs_label_values}"),
        var_query=soma.AxisQuery(value_filter=f"label in {var_label_values}"),
    ) as query:
        assert query.n_obs == len(obs_label_values)
        assert query.n_vars == len(var_label_values)
        assert query.obs().concat()["label"].to_pylist() == obs_label_values
        assert query.var().concat()["label"].to_pylist() == var_label_values


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info >= (3, 10),
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
        obs_query=soma.AxisQuery(coords=(obs_slice,)),
        var_query=soma.AxisQuery(value_filter=f"label in {var_label_values}"),
    ) as query:
        assert query.n_obs == obs_slice.stop - obs_slice.start + 1
        assert query.var().concat()["label"].to_pylist() == var_label_values

    with soma_experiment.axis_query(
        "RNA",
        obs_query=soma.AxisQuery(value_filter=f"label in {obs_label_values}"),
        var_query=soma.AxisQuery(coords=(var_slice,)),
    ) as query:
        assert query.obs().concat()["label"].to_pylist() == obs_label_values
        assert np.array_equal(
            query.var_joinids().to_numpy(),
            np.arange(var_slice.start, var_slice.stop + 1),
        )

    with soma_experiment.axis_query(
        "RNA",
        obs_query=soma.AxisQuery(
            coords=(obs_slice,), value_filter=f"label in {obs_label_values}"
        ),
        var_query=soma.AxisQuery(
            coords=(var_slice,), value_filter=f"label in {var_label_values}"
        ),
    ) as query:
        assert query.obs().concat()["label"].to_pylist() == obs_label_values
        assert query.var().concat()["label"].to_pylist() == var_label_values


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_none(soma_experiment):
    """Test query resulting in empty result"""

    with soma_experiment.axis_query(
        "RNA",
        obs_query=soma.AxisQuery(value_filter="label == 'no-such-label'"),
        var_query=soma.AxisQuery(value_filter="label == 'no-such-label'"),
    ) as query:
        obs = query.obs().concat()
        var = query.var().concat()
        assert len(obs) == 0
        assert len(var) == 0
        assert set(next(query.obs()).column_names) == {"soma_joinid", "label"}
        assert set(next(query.var()).column_names) == {"soma_joinid", "label"}

        assert len(query.obs_joinids()) == 0
        assert len(query.var_joinids()) == 0

        ad = query.to_anndata("raw")
        assert ad.n_obs == 0 and ad.n_vars == 0

        assert len(query.X("raw").tables().concat()) == 0


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars,X_layer_names", [(1001, 99, ["A"])])
def test_joinid_caching(soma_experiment):
    """
    Verify that results are the same regardless of invocation order, which
    influences caching
    """

    obs_query = soma.AxisQuery(value_filter="label in ['17', '19', '21']")
    var_query = soma.AxisQuery(coords=(slice(0, 100),))

    with soma_experiment.axis_query(
        "RNA", obs_query=obs_query, var_query=var_query
    ) as query1:
        obs = query1.obs().concat()
        var = query1.var().concat()

    with soma_experiment.axis_query(
        "RNA", obs_query=obs_query, var_query=var_query
    ) as query2:
        query2.X("A").coos().concat().to_scipy()

    with soma_experiment.axis_query(
        "RNA", obs_query=obs_query, var_query=var_query
    ) as query3:
        ad = query3.to_anndata("A", column_names={"obs": ["label"], "var": ["label"]})

    assert query1 != query2 and query2 != query3 and query1 != query3
    assert np.array_equal(obs.to_pandas().label, ad.obs.label)
    assert np.array_equal(var.to_pandas().label, ad.var.label)
    assert ad.n_obs == len(obs)
    assert ad.n_vars == len(var)


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars,X_layer_names", [(1001, 99, ["A", "B", "C"])])
def test_X_layers(soma_experiment):
    """Verify multi-layer-X handling"""
    A = pa.concat_tables(
        soma_experiment.ms["RNA"].X["A"].read((slice(None), slice(None))).tables()
    )
    B = pa.concat_tables(
        soma_experiment.ms["RNA"].X["B"].read((slice(None), slice(None))).tables()
    )

    with soma_experiment.axis_query("RNA") as query:
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
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_indexer(soma_experiment):
    """Test result indexer"""
    with soma.ExperimentAxisQuery(
        soma_experiment,
        "RNA",
        obs_query=soma.AxisQuery(coords=(slice(1, 10),)),
        var_query=soma.AxisQuery(coords=(slice(1, 10),)),
    ) as query:
        # TODO: remove this work-around once a new `somacore` is released.
        # workaround:
        indexer = getattr(query, "indexer", query._indexer)
        # future version:
        # indexer = query.indexer

        # coords outside of our query should return -1
        assert np.array_equal(
            indexer.by_obs(np.array([-1, 0, 11, 1003])),
            np.array([-1, -1, -1, -1]),
        )
        assert np.array_equal(
            indexer.by_var(np.array([-1, 0, 11, 1003])),
            np.array([-1, -1, -1, -1]),
        )

        # inside results, indexed
        assert np.array_equal(indexer.by_obs(np.array([1, 4, 2])), np.array([0, 3, 1]))
        assert np.array_equal(indexer.by_var(np.array([10, 1])), np.array([9, 0]))

        # should be able to consume multiple types
        base_arg = np.array([1, 4, 2])
        expected_result = np.array([0, 3, 1])
        for arg in (
            base_arg,
            pa.array(base_arg),
            pa.chunked_array([base_arg]),
        ):
            assert np.array_equal(indexer.by_obs(arg), expected_result)
            assert np.array_equal(indexer.by_var(arg), expected_result)


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(2833, 107)])
def test_error_corners(soma_experiment: soma.Experiment):
    """Verify a couple of error conditions / corner cases."""
    # Unknown Measurement name
    with pytest.raises(ValueError):
        soma_experiment.axis_query("no-such-measurement")

    # Unknown X layer name
    with pytest.raises(KeyError):
        with soma_experiment.axis_query("RNA") as query:
            next(query.X("no-such-layer"))

    # Unknown X layer name
    with pytest.raises(ValueError):
        with soma_experiment.axis_query("RNA") as query:
            query.to_anndata("no-such-layer")

    # Unknown obsp layer name
    with pytest.raises(ValueError):
        with soma_experiment.axis_query("RNA") as query:
            next(query.obsp("no-such-layer"))

    # Unknown varp layer name
    with pytest.raises(ValueError):
        with soma_experiment.axis_query("RNA") as query:
            next(query.varp("no-such-layer"))

    # Illegal layer name type
    for lyr_name in [True, 3, 99.3]:
        with soma_experiment.axis_query("RNA") as query:
            with pytest.raises(KeyError):
                next(query.X(lyr_name))
            with pytest.raises(ValueError):
                next(query.obsp(lyr_name))
            with pytest.raises(ValueError):
                next(query.varp(lyr_name))


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info >= (3, 10),
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
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize(
    "n_obs,n_vars,obsp_layer_names,varp_layer_names", [(1001, 99, ["foo"], ["bar"])]
)
def test_experiment_query_obsp_varp(soma_experiment):
    obs_slice = slice(3, 72)
    var_slice = slice(7, 21)
    with soma.ExperimentAxisQuery(
        soma_experiment,
        "RNA",
        obs_query=soma.AxisQuery(coords=(obs_slice,)),
        var_query=soma.AxisQuery(coords=(var_slice,)),
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
    assert soma.AxisQuery().coords == (slice(None),)
    assert soma.AxisQuery().value_filter is None
    assert soma.AxisQuery() == soma.AxisQuery(coords=(slice(None),))

    assert soma.AxisQuery(coords=(1,)).coords == (1,)
    assert soma.AxisQuery(coords=(slice(1, 2),)).coords == (slice(1, 2),)
    assert soma.AxisQuery(coords=((1, 88),)).coords == ((1, 88),)

    assert soma.AxisQuery(coords=(1, 2)).coords == (1, 2)
    assert soma.AxisQuery(coords=(slice(1, 2), slice(None))).coords == (
        slice(1, 2),
        slice(None),
    )
    assert soma.AxisQuery(coords=(slice(1, 2),)).value_filter is None

    assert soma.AxisQuery(value_filter="foo == 'bar'").value_filter == "foo == 'bar'"
    assert soma.AxisQuery(value_filter="foo == 'bar'").coords == (slice(None),)

    assert soma.AxisQuery(
        coords=(slice(1, 100),), value_filter="foo == 'bar'"
    ).coords == (
        slice(1, 100),
    )
    assert (
        soma.AxisQuery(
            coords=(slice(1, 100),), value_filter="foo == 'bar'"
        ).value_filter
        == "foo == 'bar'"
    )

    with pytest.raises(TypeError):
        soma.AxisQuery(coords=True)

    with pytest.raises(TypeError):
        soma.AxisQuery(value_filter=[])

    with pytest.raises(TypeError):
        soma.AxisQuery(coords=({},))


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


@pytest.mark.xfail(
    # see comment on test_experiment_query_all
    sys.version_info >= (3, 10),
    reason="typeguard bug #242",
)
@pytest.mark.parametrize(
    "n_obs,n_vars,obsp_layer_names,varp_layer_names", [(101, 99, ["foo"], ["bar"])]
)
def test_experiment_query_column_names(soma_experiment):
    """
    Verify that column_names is correctly handled in the various obs/var accessors.

    Returned columns should be the union of the columns specifically requested via
    column_names, and the columns implicitly requested via value_filter.
    """

    # default
    with soma_experiment.axis_query("RNA") as query:
        assert set(next(query.obs()).column_names) == {"soma_joinid", "label"}
        assert set(next(query.var()).column_names) == {"soma_joinid", "label"}
        ad = query.to_anndata("raw")
        assert set(ad.obs.keys()) == {"soma_joinid", "label"}
        assert set(ad.var.keys()) == {"soma_joinid", "label"}

    # column_names only
    with soma_experiment.axis_query("RNA") as query:
        assert set(next(query.obs(column_names=["soma_joinid"])).column_names) == {
            "soma_joinid"
        }
        assert set(next(query.var(column_names=["soma_joinid"])).column_names) == {
            "soma_joinid"
        }
        ad = query.to_anndata(
            "raw", column_names={"obs": ["soma_joinid"], "var": ["soma_joinid"]}
        )
        assert set(ad.obs.keys()) == {"soma_joinid"}
        assert set(ad.var.keys()) == {"soma_joinid"}

        assert set(next(query.obs(column_names=["label"])).column_names) == {"label"}
        assert set(next(query.var(column_names=["label"])).column_names) == {"label"}
        ad = query.to_anndata("raw", column_names={"obs": ["label"], "var": ["label"]})
        assert set(ad.obs.keys()) == {"label"}
        assert set(ad.var.keys()) == {"label"}

        assert set(
            next(query.obs(column_names=["soma_joinid", "label"])).column_names
        ) == {"soma_joinid", "label"}
        assert set(
            next(query.var(column_names=["soma_joinid", "label"])).column_names
        ) == {"soma_joinid", "label"}
        ad = query.to_anndata(
            "raw",
            column_names={
                "obs": ["soma_joinid", "label"],
                "var": ["soma_joinid", "label"],
            },
        )
        assert set(ad.obs.keys()) == {"soma_joinid", "label"}
        assert set(ad.var.keys()) == {"soma_joinid", "label"}

    # column_names and value_filter
    with soma_experiment.axis_query(
        "RNA",
        obs_query=soma.AxisQuery(
            value_filter="label in [" + ",".join(f"'{i}'" for i in range(101)) + "]"
        ),
        var_query=soma.AxisQuery(
            value_filter="label in [" + ",".join(f"'{i}'" for i in range(99)) + "]"
        ),
    ) as query:
        assert set(next(query.obs(column_names=["soma_joinid"])).column_names) == {
            "soma_joinid",
            "label",
        }
        assert set(next(query.var(column_names=["soma_joinid"])).column_names) == {
            "soma_joinid",
            "label",
        }
        ad = query.to_anndata(
            "raw", column_names={"obs": ["soma_joinid"], "var": ["soma_joinid"]}
        )
        assert set(ad.obs.keys()) == {"soma_joinid", "label"}
        assert set(ad.var.keys()) == {"soma_joinid", "label"}

        assert set(next(query.obs(column_names=["label"])).column_names) == {"label"}
        assert set(next(query.var(column_names=["label"])).column_names) == {"label"}
        ad = query.to_anndata("raw", column_names={"obs": ["label"], "var": ["label"]})
        assert set(ad.obs.keys()) == {"label"}
        assert set(ad.var.keys()) == {"label"}

        assert set(
            next(query.obs(column_names=["soma_joinid", "label"])).column_names
        ) == {"soma_joinid", "label"}
        assert set(
            next(query.var(column_names=["soma_joinid", "label"])).column_names
        ) == {"soma_joinid", "label"}
        ad = query.to_anndata(
            "raw",
            column_names={
                "obs": ["soma_joinid", "label"],
                "var": ["soma_joinid", "label"],
            },
        )
        assert set(ad.obs.keys()) == {"soma_joinid", "label"}
        assert set(ad.var.keys()) == {"soma_joinid", "label"}


"""
Fixture support & utility functions below.
"""


def add_dataframe(coll: CollectionBase, key: str, sz: int) -> None:
    df = coll.add_new_dataframe(
        key,
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


def add_sparse_array(coll: CollectionBase, key: str, shape: Tuple[int, int]) -> None:
    a = coll.add_new_sparse_ndarray(key, type=pa.float32(), shape=shape)
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
