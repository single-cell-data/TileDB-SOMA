import sys

import numpy as np
import pyarrow as pa
import pytest
from scipy import sparse

import tiledbsoma as soma
from tiledbsoma.experiment_query import experiment_query, AxisQuery


def make_dataframe(path, sz) -> soma.DataFrame:
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


def make_sparse_array(path, shape) -> soma.SparseNDArray:
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
    a.write_sparse_tensor(tensor)
    return a


def make_experiment(path, n_obs, n_vars, obs, var):
    experiment = soma.Experiment((path / "exp").as_posix()).create()
    experiment["ms"] = soma.Collection((path / "ms").as_posix()).create()
    experiment.ms["RNA"] = soma.Measurement(uri=f"{experiment.ms.uri}/RNA").create()
    experiment["obs"] = obs

    measurement = experiment.ms["RNA"]
    measurement["var"] = var
    measurement["X"] = soma.Collection(uri=f"{measurement.uri}/X").create()

    measurement.X["raw"] = make_sparse_array((path / "raw").as_posix(), (n_obs, n_vars))

    return experiment


@pytest.fixture
def obs(tmp_path, n_obs) -> soma.DataFrame:
    return make_dataframe((tmp_path / "obs").as_posix(), n_obs)


@pytest.fixture
def var(tmp_path, n_vars) -> soma.DataFrame:
    return make_dataframe((tmp_path / "var").as_posix(), n_vars)


@pytest.fixture(scope="function")
def soma_experiment(tmp_path, n_obs, n_vars, obs, var):
    return make_experiment(tmp_path, n_obs, n_vars, obs, var)


# This test fails on Python 3.10+ due to a bug in typeguard. Remove this work-around
# when the typeguard issue is fixed AND released. Underlying issue:
#   https://github.com/agronholm/typeguard/issues/242
# Related to _when_ we might see a release:
#   https://github.com/agronholm/typeguard/issues/257
#
@pytest.mark.xfail(
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(101, 11)])
def test_experiment_query_all(soma_experiment):
    assert soma_experiment.exists()

    with experiment_query(soma_experiment, "RNA") as query:
        ad = query.read_as_anndata("raw")
        assert ad.n_obs == 101 and ad.n_vars == 11

    with experiment_query(soma_experiment, "RNA") as query:
        assert query.n_obs == 101
        assert query.n_vars == 11
        assert np.array_equal(query.obs_joinids().to_numpy(), np.arange(101))
        assert np.array_equal(query.var_joinids().to_numpy(), np.arange(11))
        assert len(query.obs()) == 101
        assert len(query.var()) == 11

        ad = query.read_as_anndata("raw")
        assert ad.n_obs == query.n_obs and ad.n_vars == query.n_vars

    with experiment_query(soma_experiment, "RNA") as query:
        ad = query.read_as_anndata("raw")
        assert ad.n_obs == 101 and ad.n_vars == 11


@pytest.mark.xfail(
    sys.version_info.major == 3 and sys.version_info.minor >= 10,
    reason="typeguard bug #242",
)
@pytest.mark.parametrize("n_obs,n_vars", [(1001, 99)])
def test_experiment_query_indexer(soma_experiment):
    """Test result indexer"""
    assert soma_experiment.exists()

    with experiment_query(
        soma_experiment,
        "RNA",
        obs_query=AxisQuery(coords=(slice(1, 10),)),
        var_query=AxisQuery(coords=(slice(1, 10),)),
    ) as query:
        indexer = query.get_indexer()

        # coords outside of our query should return -1
        assert np.array_equal(
            indexer.obs_index(np.array([-1, 0, 11, 1003])), np.array([-1, -1, -1, -1])
        )
        assert np.array_equal(
            indexer.var_index(np.array([-1, 0, 11, 1003])), np.array([-1, -1, -1, -1])
        )

        # inside results, indexed
        assert np.array_equal(
            indexer.obs_index(np.array([1, 4, 2])), np.array([0, 3, 1])
        )
        assert np.array_equal(indexer.var_index(np.array([10, 1])), np.array([9, 0]))


"""
Missing tests:
* X
* read
"""
