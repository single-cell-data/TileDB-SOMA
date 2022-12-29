import numpy as np
import pyarrow as pa
import pytest
from scipy import sparse

import tiledbsoma as soma
from tiledbsoma.experiment_query import experiment_query


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


@pytest.fixture
def obs(request, tmp_path) -> soma.DataFrame:
    return make_dataframe((tmp_path / "obs").as_posix(), 100)


@pytest.fixture
def var(request, tmp_path) -> soma.DataFrame:
    return make_dataframe((tmp_path / "var").as_posix(), 10)


@pytest.fixture(scope="function")
def soma_experiment(request, tmp_path, obs, var):
    experiment = soma.Experiment((tmp_path / "exp").as_posix()).create()
    experiment["ms"] = soma.Collection((tmp_path / "ms").as_posix()).create()
    experiment.ms["RNA"] = soma.Measurement(uri=f"{experiment.ms.uri}/RNA").create()
    experiment["obs"] = obs

    measurement = experiment.ms["RNA"]
    measurement["var"] = var
    measurement["X"] = soma.Collection(uri=f"{measurement.uri}/X").create()

    measurement.X["raw"] = make_sparse_array((tmp_path / "raw").as_posix(), (100, 10))

    return experiment


def test_experiment_query_all(soma_experiment):
    assert soma_experiment.exists()

    with experiment_query(soma_experiment, "RNA") as query:
        assert query.n_obs == 100
        assert query.n_vars == 10
        assert np.array_equal(query.obs_joinids().to_numpy(), np.arange(100))
        assert np.array_equal(query.var_joinids().to_numpy(), np.arange(10))
        assert len(query.obs()) == 100
        assert len(query.var()) == 10

        ad = query.read_as_anndata("raw")
        assert ad.n_obs == 100 and ad.n_vars == 10
