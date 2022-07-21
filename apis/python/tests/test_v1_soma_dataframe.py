import pyarrow as pa

import tiledbsc.v1 as t


def test_soma_dataframe_row_indexed(tmp_path):
    sdf = t.SOMADataFrame(uri=tmp_path.as_posix())

    asch = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.string()),
        ]
    )

    # Create
    sdf.create(schema=asch, user_indexed=False)

    # Write
    for i in range(3):
        pydict = {}
        pydict["rowid__"] = [0, 1, 2, 3, 4]
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.RecordBatch.from_pydict(pydict)
        sdf.write(rb)

    # Read all
    batches = []
    for batch in sdf.read(ids="all"):
        batches.append(batch)
    # Weird thing about pyarrow RecordBatch:
    # * We should have 5 "rows" with 3 "columns"
    # * Indeed batch.num_rows is 5 and batch.num_columns is 3
    # * But len(batch) is 3
    # * If you thought `for record in record_batch` would print records ... you would be wrong -- it loops over columns
    assert len(batches) == 1
    batch = batches[0]
    assert batch.num_rows == 5
    # We should be getting back the rowid__ column as well
    # TODO assert batch.num_columns == 4
    assert batch.num_columns == 3
    # TODO assert [e.as_py() for e in list(batch['rowid__'])] == [0,1,2,3,4]
    assert [e.as_py() for e in list(batch["foo"])] == pydict["foo"]
    assert [e.as_py() for e in list(batch["bar"])] == pydict["bar"]
    assert [e.as_py() for e in list(batch["baz"])] == pydict["baz"]

    # Read ids
    batches = []
    for batch in sdf.read(ids=[1, 2]):
        batches.append(batch)
    # Weird thing about pyarrow RecordBatch:
    # * We should have 5 "rows" with 3 "columns"
    # * Indeed batch.num_rows is 5 and batch.num_columns is 3
    # * But len(batch) is 3
    # * If you thought `for record in record_batch` would print records ... you would be wrong -- it loops over columns
    assert len(batches) == 1
    batch = batches[0]
    assert batch.num_rows == 2
    # We should be getting back the rowid__ column as well
    # TODO assert batch.num_columns == 4
    assert batch.num_columns == 3
    # TODO assert [e.as_py() for e in list(batch['rowid__'])] == [0,1,2,3,4]
    assert [e.as_py() for e in list(batch["foo"])] == [20, 30]
    assert [e.as_py() for e in list(batch["bar"])] == [5.2, 6.3]
    assert [e.as_py() for e in list(batch["baz"])] == ["ball", "cat"]

    # Read ids
    batches = []
    for batch in sdf.read(ids=slice(1, 2)):
        batches.append(batch)
    # Weird thing about pyarrow RecordBatch:
    # * We should have 5 "rows" with 3 "columns"
    # * Indeed batch.num_rows is 5 and batch.num_columns is 3
    # * But len(batch) is 3
    # * If you thought `for record in record_batch` would print records ... you would be wrong -- it loops over columns
    assert len(batches) == 1
    batch = batches[0]
    assert batch.num_rows == 2
    # We should be getting back the rowid__ column as well
    # TODO assert batch.num_columns == 4
    assert batch.num_columns == 3
    # TODO assert [e.as_py() for e in list(batch['rowid__'])] == [0,1,2,3,4]
    assert [e.as_py() for e in list(batch["foo"])] == [20, 30]
    assert [e.as_py() for e in list(batch["bar"])] == [5.2, 6.3]
    assert [e.as_py() for e in list(batch["baz"])] == ["ball", "cat"]


def test_soma_dataframe_user_indexed(tmp_path):
    sdf = t.SOMADataFrame(uri=tmp_path.as_posix())

    asch = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.string()),
        ]
    )

    # Create
    sdf.create(schema=asch, user_indexed=True, index_column_names=["foo"])

    # Write
    for i in range(3):
        pydict = {}
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.RecordBatch.from_pydict(pydict)
        sdf.write(rb)

    # Read all
    batches = []
    for batch in sdf.read(ids="all"):
        batches.append(batch)
    # Weird thing about pyarrow RecordBatch:
    # * We should have 5 "rows" with 3 "columns"
    # * Indeed batch.num_rows is 5 and batch.num_columns is 3
    # * But len(batch) is 3
    # * If you thought `for record in record_batch` would print records ... you would be wrong -- it loops over columns
    assert len(batches) == 1
    batch = batches[0]
    assert batch.num_rows == 5
    # We should be getting back the rowid__ column as well
    assert batch.num_columns == 4
    assert [e.as_py() for e in list(batch["rowid__"])] == [0, 1, 2, 3, 4]
    assert [e.as_py() for e in list(batch["foo"])] == pydict["foo"]
    assert [e.as_py() for e in list(batch["bar"])] == pydict["bar"]
    assert [e.as_py() for e in list(batch["baz"])] == pydict["baz"]

    # Read ids
    batches = []
    for batch in sdf.read(ids=[30, 10]):
        batches.append(batch)
    # Weird thing about pyarrow RecordBatch:
    # * We should have 5 "rows" with 3 "columns"
    # * Indeed batch.num_rows is 5 and batch.num_columns is 3
    # * But len(batch) is 3
    # * If you thought `for record in record_batch` would print records ... you would be wrong -- it loops over columns
    assert len(batches) == 1
    batch = batches[0]
    assert batch.num_rows == 2
    # We should be getting back the rowid__ column as well
    assert batch.num_columns == 4
    assert [e.as_py() for e in list(batch["rowid__"])] == [2, 0]
    assert [e.as_py() for e in list(batch["foo"])] == [30, 10]
    assert [e.as_py() for e in list(batch["bar"])] == [6.3, 4.1]
    assert [e.as_py() for e in list(batch["baz"])] == ["cat", "apple"]
