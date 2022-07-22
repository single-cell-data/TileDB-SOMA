import numpy as np
import pyarrow as pa

import tiledbsc.v1 as t


def test_soma_dataframe_row_indexed(tmp_path):
    # ----------------------------------------------------------------
    collection = t.SOMACollection(tmp_path.as_posix())
    collection.create()

    # ----------------------------------------------------------------
    dataframe = t.SOMADataFrame("collection/obsri", parent=collection)
    arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.string()),
        ]
    )
    dataframe.create(schema=arrow_schema, user_indexed=False)

    pydict = {}
    pydict["soma_rowid"] = [0, 1, 2, 3, 4]
    pydict["foo"] = [10, 20, 30, 40, 50]
    pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
    pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
    record_batch = pa.RecordBatch.from_pydict(pydict)
    dataframe.write(record_batch)

    # ----------------------------------------------------------------
    sparse_nd_array = t.SOMASparseNdArray("collection/snda", parent=collection)
    nr = 10
    nc = 20
    sparse_nd_array.create(pa.int64(), [nr, nc])

    tensor = pa.SparseCOOTensor.from_numpy(
        data=np.asarray([7, 8, 9]),
        coords=[[1, 2], [3, 4], [5, 6]],
        shape=(nr, nc),
    )
    sparse_nd_array.write(tensor)

    # ----------------------------------------------------------------
    collection.set(dataframe)
    collection.set(sparse_nd_array)

    # ----------------------------------------------------------------
    readback = t.SOMACollection(collection.uri)
    assert readback.get_name != "nonesuch"  # TEMP TEMP TEMP
