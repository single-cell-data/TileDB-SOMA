from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma


# ----------------------------------------------------------------
def create_and_populate_obs(obs: soma.DataFrame) -> soma.DataFrame:

    obs_arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )

    # TODO: indexing option ...
    obs.create(schema=obs_arrow_schema)

    pydict = {}
    pydict["soma_joinid"] = [0, 1, 2, 3, 4]
    pydict["foo"] = [10, 20, 30, 40, 50]
    pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
    pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
    rb = pa.Table.from_pydict(pydict)
    obs.write(rb)

    return obs


# ----------------------------------------------------------------
def create_and_populate_var(var: soma.DataFrame) -> soma.DataFrame:

    var_arrow_schema = pa.schema(
        [
            ("quux", pa.large_string()),
            ("xyzzy", pa.float64()),
        ]
    )

    var.create(schema=var_arrow_schema)

    pydict = {}
    pydict["soma_joinid"] = [0, 1, 2, 3]
    pydict["quux"] = ["zebra", "yak", "xylophone", "wapiti"]
    pydict["xyzzy"] = [12.3, 23.4, 34.5, 45.6]
    rb = pa.Table.from_pydict(pydict)
    var.write(rb)

    return var


# ----------------------------------------------------------------
def create_and_populate_sparse_nd_array(
    sparse_nd_array: soma.SparseNdArray,
) -> soma.SparseNdArray:
    nr = 5
    nc = 3

    #   0 1 2
    # 0 . . 7
    # 1 . . .
    # 2 . . .
    # 3 . 8 .
    # 4 . . 9

    sparse_nd_array.create(pa.int64(), [nr, nc])

    tensor = pa.SparseCOOTensor.from_numpy(
        data=np.asarray([7, 8, 9]),
        coords=[[0, 2], [3, 1], [4, 2]],
        shape=(nr, nc),
    )
    sparse_nd_array.write_sparse_tensor(tensor)

    return sparse_nd_array


# ----------------------------------------------------------------
def test_experiment_basic(tmp_path):
    basedir = tmp_path.as_uri()

    # ----------------------------------------------------------------
    experiment = soma.Experiment(basedir)
    experiment.create()

    experiment["obs"] = create_and_populate_obs(
        soma.DataFrame(uri=urljoin(basedir, "obs"))
    )
    experiment["ms"] = soma.Collection(uri=urljoin(basedir, "ms")).create()

    measurement = soma.Measurement(uri=f"{experiment.ms.uri}/mRNA")
    measurement.create()
    experiment.ms.set("mRNA", measurement)

    measurement["var"] = create_and_populate_var(
        soma.DataFrame(uri=urljoin(measurement.uri, "var"))
    )
    measurement["X"] = soma.Collection(uri=urljoin(measurement.uri, "X")).create()

    nda = create_and_populate_sparse_nd_array(
        soma.SparseNdArray(uri=urljoin(measurement.X.uri, "data"))
    )
    measurement.X.set("data", nda)

    # ----------------------------------------------------------------
    assert len(experiment) == 2
    assert isinstance(experiment.obs, soma.DataFrame)
    assert isinstance(experiment.ms, soma.Collection)
    assert "obs" in experiment
    assert "ms" in experiment
    assert "nonesuch" not in experiment

    assert experiment.obs == experiment["obs"]
    assert experiment.ms == experiment["ms"]

    assert len(experiment.ms) == 1
    assert isinstance(experiment.ms["mRNA"], soma.Measurement)

    assert len(experiment.ms["mRNA"]) == 2
    assert "mRNA" in experiment.ms
    assert "meas2" not in experiment.ms
    assert isinstance(experiment.ms["mRNA"].var, soma.DataFrame)
    assert isinstance(experiment.ms["mRNA"].X, soma.Collection)

    assert experiment.ms["mRNA"].var == experiment["ms"]["mRNA"]["var"]
    assert experiment.ms["mRNA"].X == experiment["ms"]["mRNA"]["X"]

    assert len(experiment.ms["mRNA"].X) == 1
    assert "data" in experiment.ms["mRNA"].X
    assert "nonesuch" not in experiment.ms["mRNA"].X
    assert isinstance(experiment.ms["mRNA"].X["data"], soma.SparseNdArray)

    # >>> experiment.ms.mRNA.X.data._tiledb_open().df[:]
    #    __dim_0  __dim_1  data
    # 0        0        2     7
    # 1        3        1     8
    # 2        4        2     9

    # ----------------------------------------------------------------
    # Paths exist and are of the right type
    assert experiment.exists()
    assert experiment.obs.exists()
    assert experiment.ms.exists()
    assert experiment.ms["mRNA"].exists()
    assert experiment.ms["mRNA"].X.exists()
    assert experiment.ms["mRNA"].X["data"].exists()

    # Paths exist but are not of the right type
    assert not soma.DataFrame(experiment.uri).exists()
    assert not soma.Collection(experiment.obs.uri).exists()

    # Paths do not exist
    assert not soma.Experiment("/nonesuch/no/nope/nope/never").exists()
    assert not soma.DataFrame("/nonesuch/no/nope/nope/never").exists()

    # ----------------------------------------------------------------
    # TODO: check more things


def test_experiment_obs_type_constraint(tmp_path):
    """
    The obs and ms keys are special props, and should
    only allow a constrained set of types to be set in their slots.
    """

    se = soma.Experiment(uri=tmp_path.as_uri()).create()

    with pytest.raises(TypeError):
        se["obs"] = soma.Collection(uri=(tmp_path / "A").as_uri()).create()
    with pytest.raises(TypeError):
        se["obs"] = soma.SparseNdArray(uri=(tmp_path / "B").as_uri()).create(
            type=pa.float32(), shape=(10,)
        )
    with pytest.raises(TypeError):
        se["obs"] = soma.DenseNdArray(uri=(tmp_path / "C").as_uri()).create(
            type=pa.float32(), shape=(10,)
        )
    with pytest.raises(TypeError):
        se["obs"] = soma.Measurement(uri=(tmp_path / "D").as_uri()).create()
    se["obs"] = soma.DataFrame(uri=(tmp_path / "E").as_uri()).create(
        schema=pa.schema([("A", pa.int32())])
    )
    se["obs"] = soma.DataFrame(uri=(tmp_path / "F").as_uri()).create(
        schema=pa.schema([("A", pa.int32())]), index_column_names=["A"]
    )


def test_experiment_ms_type_constraint(tmp_path):
    se = soma.Experiment(uri=tmp_path.as_uri()).create()

    se["ms"] = soma.Collection(uri=(tmp_path / "A").as_uri()).create()
    with pytest.raises(TypeError):
        se["ms"] = soma.SparseNdArray(uri=(tmp_path / "B").as_uri()).create(
            type=pa.float32(), shape=(10,)
        )
    with pytest.raises(TypeError):
        se["ms"] = soma.DenseNdArray(uri=(tmp_path / "C").as_uri()).create(
            type=pa.float32(), shape=(10,)
        )
    with pytest.raises(TypeError):
        se["ms"] = soma.Measurement(uri=(tmp_path / "D").as_uri()).create()
    with pytest.raises(TypeError):
        se["ms"] = soma.DataFrame(uri=(tmp_path / "E").as_uri()).create(
            schema=pa.schema([("A", pa.int32())])
        )
    with pytest.raises(TypeError):
        se["ms"] = soma.DataFrame(uri=(tmp_path / "F").as_uri()).create(
            schema=pa.schema([("A", pa.int32())]), index_column_names=["A"]
        )
