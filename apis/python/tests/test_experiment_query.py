import re
from concurrent import futures
from contextlib import nullcontext
from unittest import mock

import attrs
import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from pyarrow import ArrowInvalid
from scipy import sparse
from somacore import AxisQuery, options

import tiledbsoma as soma
from tiledbsoma import (
    Experiment,
    ExperimentAxisQuery,
    SOMATileDBContext,
    pytiledbsoma,
)
from tiledbsoma._collection import CollectionBase
from tiledbsoma.experiment_query import X_as_series

from tests._util import raises_no_typeguard

# Number of features for the embeddings layer
N_FEATURES = 50


@pytest.fixture
def X_layer_names():
    return ["raw"]


@pytest.fixture
def obsp_layer_names():
    return None


@pytest.fixture
def varp_layer_names():
    return None


@pytest.fixture
def obsm_layer_names():
    return None


@pytest.fixture
def varm_layer_names():
    return None


@pytest.fixture(scope="function")
def soma_experiment(
    tmp_path,
    n_obs,
    n_vars,
    X_layer_names,
    obsp_layer_names,
    varp_layer_names,
    obsm_layer_names,
    varm_layer_names,
) -> Experiment:
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

        if obsm_layer_names:
            obsm = rna.add_new_collection("obsm")
            for obsm_layer_name in obsm_layer_names:
                add_sparse_array(obsm, obsm_layer_name, (n_obs, N_FEATURES))

        if varm_layer_names:
            varm = rna.add_new_collection("varm")
            for varm_layer_name in varm_layer_names:
                add_sparse_array(varm, varm_layer_name, (n_vars, N_FEATURES))

    return Experiment.open((tmp_path / "exp").as_posix())


def get_soma_experiment_with_context(soma_experiment, context) -> Experiment:
    soma_experiment.close()
    return Experiment.open(soma_experiment.uri, context=context)


@pytest.mark.parametrize("n_obs,n_vars,X_layer_names", [(101, 11, ("raw", "extra"))])
def test_experiment_query_all(soma_experiment):
    """Test a query with default obs_query / var_query -- i.e., query all."""
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
        raw = query.X("raw")
        blockwise = raw.blockwise(axis=0, reindex_disable_on_axis=[1])
        assert sparse.vstack([sp for sp, _ in blockwise.scipy()]).shape == (
            query.n_obs,
            query.n_vars,
        )

        # read as anndata
        ad = query.to_anndata("raw")
        assert ad.n_obs == query.n_obs and ad.n_vars == query.n_vars

        assert set(ad.obs.keys()) == {"soma_joinid", "label"}
        assert set(ad.var.keys()) == {"soma_joinid", "label"}

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


@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_value_filter2(soma_experiment):
    """Test query by value filter"""
    obs_label_values = ["3", "7", "38", "99"]
    var_label_values = ["18", "34", "67"]
    with soma_experiment.axis_query(
        "RNA",
        obs_query=soma.AxisQuery(value_filter=f"label not in {obs_label_values}"),
        var_query=soma.AxisQuery(value_filter=f"label not in {var_label_values}"),
    ) as query:
        assert query.n_obs == soma_experiment.obs.count - len(obs_label_values)
        assert query.n_vars == soma_experiment.ms["RNA"].var.count - len(
            var_label_values
        )
        all_obs_values = set(
            soma_experiment.obs.read(column_names=["label"])
            .concat()
            .to_pandas()["label"]
        )
        all_var_values = set(
            soma_experiment.ms["RNA"]
            .var.read(column_names=["label"])
            .concat()
            .to_pandas()["label"]
        )
        qry_obs_values = set(query.obs().concat()["label"].to_pylist())
        qry_var_values = set(query.var().concat()["label"].to_pylist())
        assert qry_obs_values == all_obs_values.difference(set(obs_label_values))
        assert qry_var_values == all_var_values.difference(set(var_label_values))


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


@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_batch_size(soma_experiment):
    """
    batch_size is currently not supported by this implementation of SOMA.
    This test merely verifies that the batch_size parameter is accepted
    but as a no-op.
    """
    with ExperimentAxisQuery(soma_experiment, "RNA") as query:
        tbls = query.obs(batch_size=options.BatchSize(count=100))
        assert len(list(tbls)) == 1  # batch_size currently not implemented


@pytest.mark.parametrize("n_obs,n_vars", [(10, 10)])
def test_experiment_query_partitions(soma_experiment):
    """
    partitions is currently not supported by this implementation of SOMA.
    This test checks if a ValueError is raised if a partitioning is requested.
    """
    with ExperimentAxisQuery(soma_experiment, "RNA") as query:
        with pytest.raises(ValueError):
            query.obs(partitions=options.IOfN(i=0, n=3)).concat()

        with pytest.raises(ValueError):
            query.var(partitions=options.IOfN(i=0, n=3)).concat()

        with pytest.raises(ValueError):
            query.X("raw", partitions=options.IOfN(i=0, n=3)).concat()


@pytest.mark.parametrize("n_obs,n_vars", [(10, 10)])
def test_experiment_query_result_order(soma_experiment):
    with ExperimentAxisQuery(soma_experiment, "RNA") as query:
        # Since obs is 1-dimensional, row-major and column-major should be the same
        obs_data_row_major = (
            query.obs(result_order="row-major").concat()["label"].to_numpy()
        )
        obs_data_col_major = (
            query.obs(result_order="column-major").concat()["label"].to_numpy()
        )
        assert np.array_equal(obs_data_row_major, obs_data_col_major)
        assert np.array_equal(np.sort(obs_data_row_major), obs_data_row_major)
        assert np.array_equal(np.sort(obs_data_col_major), obs_data_col_major)

        # The same for var
        var_data_row_major = (
            query.var(result_order="row-major").concat()["label"].to_numpy()
        )
        var_data_col_major = (
            query.var(result_order="column-major").concat()["label"].to_numpy()
        )
        assert np.array_equal(var_data_row_major, var_data_col_major)
        assert np.array_equal(np.sort(var_data_row_major), var_data_row_major)
        assert np.array_equal(np.sort(var_data_col_major), var_data_col_major)

        X_tbl = query.X("raw", result_order="row-major").tables().concat()
        row = X_tbl["soma_dim_0"].to_numpy()
        data_row_major = X_tbl["soma_data"].to_numpy()
        assert np.array_equal(np.sort(row), row)

        X_tbl = query.X("raw", result_order="column-major").tables().concat()
        col = X_tbl["soma_dim_1"].to_numpy()
        data_col_major = X_tbl["soma_data"].to_numpy()
        assert np.array_equal(np.sort(col), col)
        assert not np.array_equal(data_row_major, data_col_major)


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


@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_axis_query_with_none(soma_experiment):
    """Test query by value filter"""
    obs_label_values = ["3", "7", "38", "99"]

    with ExperimentAxisQuery(
        experiment=soma_experiment,
        measurement_name="RNA",
        obs_query=soma.AxisQuery(value_filter=f"label in {obs_label_values}"),
    ) as query:
        assert query.n_obs == len(obs_label_values)
        assert query.obs().concat()["label"].to_pylist() == obs_label_values


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


@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_indexer(soma_experiment):
    """Test result indexer"""
    with ExperimentAxisQuery(
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


@pytest.mark.parametrize("n_obs,n_vars", [(2833, 107)])
def test_error_corners(soma_experiment: Experiment):
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
            with raises_no_typeguard(KeyError):
                next(query.X(lyr_name))
            with raises_no_typeguard(ValueError):
                next(query.obsp(lyr_name))
            with raises_no_typeguard(ValueError):
                next(query.varp(lyr_name))


@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_query_cleanup(soma_experiment: soma.Experiment):
    """
    Verify soma.Experiment.query works as context manager and stand-alone,
    and that it cleans up correctly.
    """
    from contextlib import closing

    context = SOMATileDBContext()
    soma_experiment = get_soma_experiment_with_context(soma_experiment, context)

    with soma_experiment.axis_query("RNA") as query:
        assert query.n_obs == 1001
        assert query.n_vars == 99
        assert query.to_anndata("raw") is not None

    with closing(soma_experiment.axis_query("RNA")) as query:
        assert query.to_anndata("raw") is not None


@pytest.mark.parametrize(
    "n_obs,n_vars,obsp_layer_names,varp_layer_names,obsm_layer_names,varm_layer_names",
    [(1001, 99, ["foo"], ["bar"], ["baz"], ["quux"])],
)
def test_experiment_query_obsp_varp_obsm_varm(soma_experiment):
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

        with raises_no_typeguard(ValueError):
            next(query.obsp("no-such-layer"))

        with pytest.raises(ValueError):
            next(query.varp("no-such-layer"))

        with pytest.raises(ValueError):
            next(query.obsm("no-such-layer"))

        with pytest.raises(ValueError):
            next(query.varm("no-such-layer"))

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

        assert (
            query.obsm("baz").tables().concat()
            == soma_experiment.ms["RNA"]
            .obsm["baz"]
            .read((obs_slice, range(N_FEATURES)))
            .tables()
            .concat()
        )

        assert (
            query.varm("quux").tables().concat()
            == soma_experiment.ms["RNA"]
            .varm["quux"]
            .read((var_slice, range(N_FEATURES)))
            .tables()
            .concat()
        )


@pytest.mark.parametrize(
    "n_obs,n_vars,obsm_layer_names,varm_layer_names", [(1001, 99, ["foo"], ["bar"])]
)
def test_experiment_query_to_anndata_obsm_varm(soma_experiment):
    with soma_experiment.axis_query("RNA") as query:
        ad = query.to_anndata("raw", obsm_layers=["foo"], varm_layers=["bar"])
        assert set(ad.obsm.keys()) == {"foo"}
        obsm = ad.obsm["foo"]
        assert isinstance(obsm, np.ndarray)
        assert obsm.shape == (query.n_obs, N_FEATURES)

        assert np.array_equal(
            query.obsm("foo").coos().concat().to_scipy().todense(), obsm
        )

        assert set(ad.varm.keys()) == {"bar"}
        varm = ad.varm["bar"]
        assert isinstance(varm, np.ndarray)
        assert varm.shape == (query.n_vars, N_FEATURES)
        assert np.array_equal(
            query.varm("bar").coos().concat().to_scipy().todense(), varm
        )


@pytest.mark.parametrize(
    "n_obs,n_vars,obsp_layer_names,varp_layer_names", [(1001, 99, ["foo"], ["bar"])]
)
def test_experiment_query_to_anndata_obsp_varp(soma_experiment):
    with soma_experiment.axis_query("RNA") as query:
        ad = query.to_anndata("raw", obsp_layers=["foo"], varp_layers=["bar"])
        assert set(ad.obsp.keys()) == {"foo"}
        obsp = ad.obsp["foo"]
        assert isinstance(obsp, sparse.spmatrix)
        assert sparse.isspmatrix_csr(obsp)
        assert obsp.shape == (query.n_obs, query.n_obs)

        assert (query.obsp("foo").coos().concat().to_scipy() != obsp).nnz == 0
        assert np.array_equal(
            query.obsp("foo").coos().concat().to_scipy().todense(), obsp.todense()
        )

        assert set(ad.varp.keys()) == {"bar"}
        varp = ad.varp["bar"]
        assert isinstance(varp, sparse.spmatrix)
        assert sparse.isspmatrix_csr(varp)
        assert varp.shape == (query.n_vars, query.n_vars)
        assert (query.varp("bar").coos().concat().to_scipy() != varp).nnz == 0
        assert np.array_equal(
            query.varp("bar").coos().concat().to_scipy().todense(), varp.todense()
        )


def test_axis_query():
    """Basic test of the AxisQuery class"""
    assert AxisQuery().coords == ()
    assert AxisQuery().value_filter is None
    assert AxisQuery() == AxisQuery(coords=())

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
    assert AxisQuery(value_filter="foo == 'bar'").coords == ()

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
        obs_query=AxisQuery(
            value_filter="label in [" + ",".join(f"'{i}'" for i in range(101)) + "]"
        ),
        var_query=AxisQuery(
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


@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_mp_disjoint_arrow_coords(soma_experiment):
    """
    Verify Pyarrow join ids that are offset are correctly handled.
    """
    pa.array(range(30))
    slices = [
        pa.array(range(10)),
        pa.array(range(10, 20)),
        pa.array(range(20, 30)),
    ]

    for ids in slices:
        with soma_experiment.axis_query(
            "RNA",
            obs_query=AxisQuery(coords=(ids,)),
        ) as query:
            assert query.obs_joinids() == ids


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
        domain=[[0, sz - 1]],
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


def add_sparse_array(coll: CollectionBase, key: str, shape: tuple[int, int]) -> None:
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


@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_uses_threadpool_from_context(soma_experiment):
    """
    Verify that ExperimentAxisQuery uses the threadpool from the context
    """

    pool = mock.Mock(wraps=futures.ThreadPoolExecutor(max_workers=4))
    pool.submit.assert_not_called()

    context = SOMATileDBContext(threadpool=pool)
    soma_experiment = get_soma_experiment_with_context(soma_experiment, context=context)

    with soma_experiment.axis_query("RNA") as query:
        # to_anndata uses the threadpool
        adata = query.to_anndata(X_name="raw")
        assert adata is not None

        pool.submit.assert_called()


def test_empty_categorical_query(conftest_pbmc_small_exp):
    q = conftest_pbmc_small_exp.axis_query(
        measurement_name="RNA", obs_query=AxisQuery(value_filter='groups == "g1"')
    )
    obs = q.obs().concat()
    assert len(obs) == 44

    adata = q.to_anndata(column_names={"obs": ["groups"]}, X_name="data")
    cat = adata.obs["groups"].cat.categories
    assert "g1" in cat
    assert "g2" in cat

    adata = q.to_anndata(
        column_names={"obs": ["groups"]}, X_name="data", drop_levels=True
    )
    cat = adata.obs["groups"].cat.categories
    assert "g1" in cat
    # Unused categories should not appear
    assert "g2" not in cat

    q = conftest_pbmc_small_exp.axis_query(
        measurement_name="RNA", obs_query=AxisQuery(value_filter='groups == "foo"')
    )
    # Empty query on a categorical column raised ArrowInvalid before TileDB 2.21; see https://github.com/single-cell-data/TileDB-SOMA/pull/2299
    m = re.fullmatch(r"libtiledb=(\d+\.\d+\.\d+)", pytiledbsoma.version())
    version = m.group(1).split(".")
    major, minor = int(version[0]), int(version[1])

    ctx = nullcontext() if (major, minor) >= (2, 21) else pytest.raises(ArrowInvalid)
    with ctx:
        obs = q.obs().concat()
        assert len(obs) == 0


@attrs.define(frozen=True)
class IHaveObsVarStuff:
    obs: int
    var: int
    the_obs_suf: str
    the_var_suf: str
