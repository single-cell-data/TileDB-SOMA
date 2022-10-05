import numpy as np
import pyarrow as pa

import tiledbsoma as t


# ----------------------------------------------------------------
def create_and_populate_obs(obs: t.SOMADataFrame) -> None:

    obs_arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.string()),
        ]
    )

    # TODO: indexing option ...
    obs.create(schema=obs_arrow_schema)

    pydict = {}
    pydict["soma_rowid"] = [0, 1, 2, 3, 4]
    pydict["foo"] = [10, 20, 30, 40, 50]
    pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
    pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
    rb = pa.Table.from_pydict(pydict)
    obs.write(rb)


# ----------------------------------------------------------------
def create_and_populate_var(var: t.SOMADataFrame) -> None:

    var_arrow_schema = pa.schema(
        [
            ("quux", pa.string()),
            ("xyzzy", pa.float64()),
        ]
    )

    var.create(schema=var_arrow_schema)

    pydict = {}
    pydict["soma_rowid"] = [0, 1, 2, 3]
    pydict["quux"] = ["zebra", "yak", "xylophone", "wapiti"]
    pydict["xyzzy"] = [12.3, 23.4, 34.5, 45.6]
    rb = pa.Table.from_pydict(pydict)
    var.write(rb)


# ----------------------------------------------------------------
def create_and_populate_sparse_nd_array(sparse_nd_array: t.SOMASparseNdArray) -> None:
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


# ----------------------------------------------------------------
def test_soma_experiment_basic(tmp_path):
    basedir = tmp_path.as_posix()

    # ----------------------------------------------------------------
    experiment = t.SOMAExperiment(basedir)
    experiment.create()

    create_and_populate_obs(experiment.obs)
    experiment.set(experiment.obs)

    experiment.ms.create()
    experiment.set(experiment.ms)

    measurement = t.SOMAMeasurement(uri=f"{experiment.ms.uri}/mRNA")
    measurement.create()
    experiment.ms.set(measurement)

    create_and_populate_var(measurement.var)
    measurement.set(measurement.var)

    measurement.X.create()
    measurement.set(measurement.X)

    nda = t.SOMASparseNdArray(uri=f"{measurement.X.uri}/data")
    create_and_populate_sparse_nd_array(nda)
    measurement.X.set(nda)

    # ----------------------------------------------------------------
    assert len(experiment) == 2
    assert isinstance(experiment.obs, t.SOMADataFrame)
    assert isinstance(experiment.ms, t.SOMACollectionBase)
    assert "obs" in experiment
    assert "ms" in experiment
    assert "nonesuch" not in experiment

    assert experiment.obs == experiment["obs"]
    assert experiment.ms == experiment["ms"]

    assert len(experiment.ms) == 1
    assert isinstance(experiment.ms["mRNA"], t.SOMAMeasurement)

    assert len(experiment.ms["mRNA"]) == 2
    assert "mRNA" in experiment.ms
    assert "meas2" not in experiment.ms
    assert isinstance(experiment.ms["mRNA"].var, t.SOMADataFrame)
    assert isinstance(experiment.ms["mRNA"].X, t.SOMACollectionBase)

    assert experiment.ms["mRNA"].var == experiment["ms"]["mRNA"]["var"]
    assert experiment.ms["mRNA"].X == experiment["ms"]["mRNA"]["X"]

    assert len(experiment.ms["mRNA"].X) == 1
    assert "data" in experiment.ms["mRNA"].X
    assert "nonesuch" not in experiment.ms["mRNA"].X
    assert isinstance(experiment.ms["mRNA"].X["data"], t.SOMASparseNdArray)

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
    assert not t.SOMADataFrame(experiment.uri).exists()
    assert not t.SOMACollection(experiment.obs.uri).exists()

    # Paths do not exist
    assert not t.SOMAExperiment("/nonesuch/no/nope/nope/never").exists()
    assert not t.SOMADataFrame("/nonesuch/no/nope/nope/never").exists()

    # ----------------------------------------------------------------
    # TODO: check more things
