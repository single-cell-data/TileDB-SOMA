import pyarrow as pa
import pytest

import tiledbsoma as t


def test_soma_indexed_dataframe(tmp_path):
    sdf = t.SOMAIndexedDataFrame(uri=tmp_path.as_posix())

    asch = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.string()),
        ]
    )

    # Create
    sdf.create(schema=asch, index_column_names=["foo"])

    # Write
    for _ in range(3):
        pydict = {}
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.RecordBatch.from_pydict(pydict)
        sdf.write(rb)

    # Read all
    batch = sdf.read_all()
    # Weird thing about pyarrow RecordBatch:
    # * We should have 5 "rows" with 3 "columns"
    # * Indeed batch.num_rows is 5 and batch.num_columns is 3
    # * But len(batch) is 3
    # * If you thought `for record in record_batch` would print records ... you would be wrong -- it
    #   loops over columns
    assert batch.num_rows == 5
    assert batch.num_columns == 3
    assert [e.as_py() for e in list(batch["foo"])] == pydict["foo"]
    assert [e.as_py() for e in list(batch["bar"])] == pydict["bar"]
    assert [e.as_py() for e in list(batch["baz"])] == pydict["baz"]

    # Read ids
    batch = sdf.read_all(ids=[30, 10])
    # Weird thing about pyarrow RecordBatch:
    # * We should have 5 "rows" with 3 "columns"
    # * Indeed batch.num_rows is 5 and batch.num_columns is 3
    # * But len(batch) is 3
    # * If you thought `for record in record_batch` would print records ... you would be wrong -- it
    #   loops over columns
    assert batch.num_rows == 2
    assert batch.num_columns == 3
    assert sorted([e.as_py() for e in list(batch["foo"])]) == [10, 30]
    assert sorted([e.as_py() for e in list(batch["bar"])]) == [4.1, 6.3]
    assert sorted([e.as_py() for e in list(batch["baz"])]) == ["apple", "cat"]


@pytest.fixture
def simple_soma_indexed_data_frame(tmp_path):
    """
    A pytest fixture which creates a simple SOMAIndexedDataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            ("index", pa.uint64()),
            ("A", pa.int64()),
            ("B", pa.float64()),
            ("C", pa.string()),
        ]
    )
    index_column_names = ["index"]
    sdf = t.SOMAIndexedDataFrame(uri=tmp_path.as_posix())
    sdf.create(schema=schema, index_column_names=index_column_names)

    data = {
        "index": [0, 1, 2, 3],
        "A": [10, 11, 12, 13],
        "B": [100.1, 200.2, 300.3, 400.4],
        "C": ["this", "is", "a", "test"],
    }
    n_data = len(data["index"])
    rb = pa.RecordBatch.from_pydict(data)
    sdf.write(rb)
    yield (schema, sdf, n_data, index_column_names)
    sdf.delete()


@pytest.mark.parametrize(
    "ids",
    [
        None,
        [
            0,
        ],
        [1, 3],
    ],
)
@pytest.mark.parametrize(
    "col_names",
    [
        ["A"],
        ["B"],
        ["A", "B"],
        ["index"],
        ["index", "A", "B", "C"],
        None,
    ],
)
def test_SOMAIndexedDataFrame_read_column_names(
    simple_soma_indexed_data_frame, ids, col_names
):
    schema, sdf, n_data, index_column_names = simple_soma_indexed_data_frame
    assert sdf.exists()

    print(schema)
    print(sdf)
    print(n_data)
    print(index_column_names)
    print(col_names)
    print(ids)

    def _check_tbl(tbl, col_names, ids):
        print(tbl)
        assert tbl.num_columns == (
            len(schema.names) if col_names is None else len(col_names)
        )
        assert tbl.num_rows == (n_data if ids is None else len(ids))
        assert tbl.schema == pa.schema(
            [
                schema.field(f)
                for f in (col_names if col_names is not None else schema.names)
            ]
        )

    _check_tbl(
        pa.Table.from_batches(sdf.read(ids=ids, column_names=col_names)),
        col_names,
        ids,
    )
    _check_tbl(
        pa.Table.from_batches([sdf.read_all(column_names=col_names)]),
        col_names,
        None,
    )

    # TODO: currently unimplemented. Enable tests when issue #329 is resolved.
    #
    # _check_tbl(
    #     pa.Table.from_pandas(
    #         pd.concat(sdf.read_as_pandas(ids=ids, column_names=col_names))
    #     ),
    #     col_names,
    #     ids,
    # )
    # _check_tbl(
    #     pa.Table.from_pandas(sdf.read_as_pandas_all(column_names=col_names)),
    #     col_names,
    #     None,
    # )
