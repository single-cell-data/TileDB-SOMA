import os

import numpy as np
import pyarrow as pa

import tiledbsc.v1 as t

# ----------------------------------------------------------------
def create_and_populate_dataframe(dataframe: t.SOMADataFrame) -> None:

    arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.string()),
        ]
    )

    dataframe.create(schema=arrow_schema, indexed=False)

    pydict = {}
    if not dataframe.get_is_indexed():
        pydict["soma_rowid"] = [0, 1, 2, 3, 4]
    pydict["foo"] = [10, 20, 30, 40, 50]
    pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
    pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
    rb = pa.RecordBatch.from_pydict(pydict)
    dataframe.write(rb)

# ----------------------------------------------------------------
def create_and_populate_sparse_nd_array(sparse_nd_array: t.SOMASparseNdArray) -> None:
    nr = 10
    nc = 20
    sparse_nd_array.create(pa.int64(), [nr, nc])

    tensor = pa.SparseCOOTensor.from_numpy(
        data=np.asarray([7, 8, 9]),
        coords=[[0, 1], [2, 3], [3, 4]],
        shape=(nr, nc),
    )
    sparse_nd_array.write(tensor)

# ----------------------------------------------------------------
def test_soma_collection_basic(tmp_path):
    basedir = tmp_path.as_posix()
    collection = t.SOMACollection(basedir)

    # ----------------------------------------------------------------
    collection.create()

    dataframe = t.SOMADataFrame(os.path.join(basedir, "sdf"), parent=collection)
    create_and_populate_dataframe(dataframe)

    sparse_nd_array = t.SOMASparseNdArray(
        os.path.join(basedir, "snda"), parent=collection
    )
    create_and_populate_sparse_nd_array(sparse_nd_array)

    collection.set(dataframe)
    collection.set(sparse_nd_array)

    # ----------------------------------------------------------------
    readback_collection = t.SOMACollection(collection.get_uri())
    assert len(readback_collection) == 2

    readback_dataframe = readback_collection.get("sdf")
    with readback_dataframe._tiledb_open() as A:
        assert len(A.df[:]) == 5

    readback_sparse_nd_array = readback_collection.get("snda")
    with readback_sparse_nd_array._tiledb_open() as A:
        assert len(A.df[:]) == 3
