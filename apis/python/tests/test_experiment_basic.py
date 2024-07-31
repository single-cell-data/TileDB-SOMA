import datetime
from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma import _factory

from tests._util import raises_no_typeguard


# ----------------------------------------------------------------
def create_and_populate_obs(uri: str) -> soma.DataFrame:
    obs_arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )

    # TODO: indexing option ...
    with soma.DataFrame.create(uri, schema=obs_arrow_schema) as obs:
        pydict = {}
        pydict["soma_joinid"] = [0, 1, 2, 3, 4]
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.Table.from_pydict(pydict)
        obs.write(rb)

    return _factory.open(uri)


# ----------------------------------------------------------------
def create_and_populate_var(uri: str) -> soma.DataFrame:
    var_arrow_schema = pa.schema(
        [
            ("quux", pa.large_string()),
            ("xyzzy", pa.float64()),
        ]
    )

    with soma.DataFrame.create(uri, schema=var_arrow_schema) as var:
        pydict = {}
        pydict["soma_joinid"] = [0, 1, 2, 3]
        pydict["quux"] = ["zebra", "yak", "xylophone", "wapiti"]
        pydict["xyzzy"] = [12.3, 23.4, 34.5, 45.6]
        rb = pa.Table.from_pydict(pydict)
        var.write(rb)

    return _factory.open(uri)


# ----------------------------------------------------------------
def create_and_populate_sparse_nd_array(uri: str) -> soma.SparseNDArray:
    nr = 5
    nc = 3

    #   0 1 2
    # 0 . . 7
    # 1 . . .
    # 2 . . .
    # 3 . 8 .
    # 4 . . 9

    with soma.SparseNDArray.create(
        uri, type=pa.int64(), shape=[nr, nc]
    ) as sparse_nd_array:
        tensor = pa.SparseCOOTensor.from_numpy(
            data=np.asarray([7, 8, 9]),
            coords=[[0, 2], [3, 1], [4, 2]],
            shape=(nr, nc),
        )
        sparse_nd_array.write(tensor)

    return _factory.open(uri)


# ----------------------------------------------------------------
def test_experiment_basic(tmp_path):
    basedir = tmp_path.as_uri()

    # ----------------------------------------------------------------
    experiment = soma.Experiment.create(basedir)
    assert soma.Experiment.exists(basedir)
    assert not soma.Collection.exists(basedir)

    experiment["obs"] = create_and_populate_obs(urljoin(basedir, "obs"))
    ms = experiment.add_new_collection("ms", soma.Collection)
    measurement = ms.add_new_collection("RNA", soma.Measurement)
    assert soma.Measurement.exists(measurement.uri)
    assert not soma.Collection.exists(measurement.uri)

    measurement["var"] = create_and_populate_var(urljoin(measurement.uri, "var"))

    x = measurement.add_new_collection("X")

    nda = create_and_populate_sparse_nd_array(urljoin(x.uri, "data"))
    x.set("data", nda, use_relative_uri=False)

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
    assert isinstance(experiment.ms["RNA"], soma.Measurement)

    assert len(experiment.ms["RNA"]) == 2
    assert "RNA" in experiment.ms
    assert "meas2" not in experiment.ms
    assert isinstance(experiment.ms["RNA"].var, soma.DataFrame)
    assert isinstance(experiment.ms["RNA"].X, soma.Collection)

    assert experiment.ms["RNA"].var == experiment["ms"]["RNA"]["var"]
    assert experiment.ms["RNA"].X == experiment["ms"]["RNA"]["X"]

    assert len(experiment.ms["RNA"].X) == 1
    assert "data" in experiment.ms["RNA"].X
    assert "nonesuch" not in experiment.ms["RNA"].X
    assert isinstance(experiment.ms["RNA"].X["data"], soma.SparseNDArray)

    # >>> experiment.ms.RNA.X.data._tiledb_open().df[:]
    #    __dim_0  __dim_1  data
    # 0        0        2     7
    # 1        3        1     8
    # 2        4        2     9

    # ----------------------------------------------------------------
    # Paths exist and are of the right type
    assert experiment is not None
    assert experiment.obs is not None
    assert experiment.ms is not None
    assert experiment.ms["RNA"] is not None
    assert experiment.ms["RNA"].X is not None
    assert experiment.ms["RNA"].X["data"] is not None

    # Paths exist but are not of the right type
    # TODO: Restore once type verification is back
    # assert not soma.DataFrame(experiment.uri).exists()
    # assert not soma.Collection(experiment.obs.uri).exists()

    # Paths do not exist
    with pytest.raises(soma.DoesNotExistError):
        soma.Experiment.open("/nonesuch/no/nope/nope/never")
    with pytest.raises(soma.DoesNotExistError):
        soma.DataFrame.open("/nonesuch/no/nope/nope/never")

    # ----------------------------------------------------------------
    # TODO: check more things


def test_experiment_obs_type_constraint(tmp_path):
    """
    The obs and ms keys are special props, and should
    only allow a constrained set of types to be set in their slots.
    """

    se = soma.Experiment.create(tmp_path.as_uri())

    with pytest.raises(TypeError):
        se["obs"] = soma.Collection.create((tmp_path / "A").as_uri())
    with pytest.raises(TypeError):
        se["obs"] = soma.SparseNDArray.create(
            (tmp_path / "B").as_uri(), type=pa.float32(), shape=(10,)
        )
    with pytest.raises(TypeError):
        se["obs"] = soma.DenseNDArray.create(
            (tmp_path / "C").as_uri(), type=pa.float32(), shape=(10,)
        )
    with pytest.raises(TypeError):
        se["obs"] = soma.Measurement.create((tmp_path / "D").as_uri())
    se["obs"] = soma.DataFrame.create(
        (tmp_path / "E").as_uri(), schema=pa.schema([("A", pa.int32())])
    )


def test_experiment_ms_type_constraint(tmp_path):
    se = soma.Experiment.create(tmp_path.as_uri())

    se["ms"] = soma.Collection.create((tmp_path / "A").as_uri())
    with pytest.raises(TypeError):
        se["ms"] = soma.SparseNDArray.create(
            (tmp_path / "B").as_uri(), type=pa.float32(), shape=(10,)
        )
    with pytest.raises(TypeError):
        se["ms"] = soma.DenseNDArray.create(
            (tmp_path / "C").as_uri(), type=pa.float32(), shape=(10,)
        )
    with pytest.raises(TypeError):
        se["ms"] = soma.Measurement.create((tmp_path / "D").as_uri())
    with pytest.raises(TypeError):
        se["ms"] = soma.DataFrame.create(
            (tmp_path / "E").as_uri(), schema=pa.schema([("A", pa.int32())])
        )
    with pytest.raises(TypeError):
        se["ms"] = soma.DataFrame.create(
            (tmp_path / "F").as_uri(),
            schema=pa.schema([("A", pa.int32())]),
            index_column_names=["A"],
        )


def test_experiment_reopen(tmp_path):
    # Ensure that reopen uses the correct mode
    soma.Experiment.create(tmp_path.as_uri(), tiledb_timestamp=1)

    with soma.Experiment.open(tmp_path.as_posix(), "r", tiledb_timestamp=1) as exp1:
        with raises_no_typeguard(ValueError):
            exp1.reopen("invalid")

        with exp1.reopen("w", tiledb_timestamp=2) as exp2:
            with exp2.reopen("r", tiledb_timestamp=3) as exp3:
                assert exp1.mode == "r"
                assert exp2.mode == "w"
                assert exp3.mode == "r"
                assert exp1.tiledb_timestamp_ms == 1
                assert exp2.tiledb_timestamp_ms == 2
                assert exp3.tiledb_timestamp_ms == 3

    ts1 = datetime.datetime(2023, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    ts2 = datetime.datetime(2024, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    with soma.Experiment.open(tmp_path.as_posix(), "r", tiledb_timestamp=ts1) as exp1:
        with exp1.reopen("r", tiledb_timestamp=ts2) as exp2:
            assert exp1.mode == "r"
            assert exp2.mode == "r"
            assert exp1.tiledb_timestamp == ts1
            assert exp2.tiledb_timestamp == ts2

    with soma.Experiment.open(tmp_path.as_posix(), "w") as exp1:
        with exp1.reopen("w", tiledb_timestamp=None) as exp2:
            with exp2.reopen("w") as exp3:
                assert exp1.mode == "w"
                assert exp2.mode == "w"
                assert exp3.mode == "w"
                now = datetime.datetime.now(datetime.timezone.utc)
                assert exp1.tiledb_timestamp <= now
                assert exp2.tiledb_timestamp <= now
                assert exp2.tiledb_timestamp <= now
