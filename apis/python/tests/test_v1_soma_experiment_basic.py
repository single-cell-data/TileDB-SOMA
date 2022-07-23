import numpy as np
import pyarrow as pa

import tiledbsc.v1 as t


# ----------------------------------------------------------------
def create_and_populate_obs(obs: t.SOMADataFrame) -> None:

    obs_arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.string()),
        ]
    )

    obs.create(schema=obs_arrow_schema, indexed=False)

    pydict = {}
    if not obs.get_is_indexed():
        pydict["soma_rowid"] = [0, 1, 2, 3, 4]
    pydict["foo"] = [10, 20, 30, 40, 50]
    pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
    pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
    rb = pa.RecordBatch.from_pydict(pydict)
    obs.write(rb)


# ----------------------------------------------------------------
def create_and_populate_var(var: t.SOMADataFrame) -> None:

    var_arrow_schema = pa.schema(
        [
            ("quux", pa.string()),
            ("xyzzy", pa.float64()),
        ]
    )

    var.create(schema=var_arrow_schema, indexed=False)

    pydict = {}
    if not var.get_is_indexed():
        pydict["soma_rowid"] = [0, 1, 2, 3]
    pydict["quux"] = ["zebra", "yak", "xylophone", "wapiti"]
    pydict["xyzzy"] = [12.3, 23.4, 34.5, 45.6]
    rb = pa.RecordBatch.from_pydict(pydict)
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
    sparse_nd_array.write(tensor)


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

    measurement = t.SOMAMeasurement(uri=f"{experiment.ms.get_uri()}/meas1")
    measurement.create()
    experiment.ms.set(measurement)

    create_and_populate_var(measurement.var)
    measurement.set(measurement.var)

    measurement.X.create()
    measurement.set(measurement.X)

    nda = t.SOMASparseNdArray(uri=f"{measurement.X.get_uri()}/data")
    create_and_populate_sparse_nd_array(nda)
    measurement.X.set(nda)

    # ----------------------------------------------------------------
    assert len(experiment) == 2
    assert isinstance(experiment.obs, t.SOMADataFrame)
    assert isinstance(experiment.ms, t.SOMACollection)

    assert len(experiment.ms) == 1
    assert isinstance(experiment.ms["meas1"], t.SOMAMeasurement)

    assert len(experiment.ms["meas1"]) == 2
    assert isinstance(experiment.ms["meas1"].var, t.SOMADataFrame)
    assert isinstance(experiment.ms["meas1"].X, t.SOMACollection)

    assert len(experiment.ms["meas1"].X) == 1
    assert isinstance(experiment.ms["meas1"].X["data"], t.SOMASparseNdArray)

    # >>> experiment.ms.meas1.X.data._tiledb_open().df[:]
    #    __dim_0  __dim_1  data
    # 0        0        2     7
    # 1        3        1     8
    # 2        4        2     9

    # TODO: check more things
