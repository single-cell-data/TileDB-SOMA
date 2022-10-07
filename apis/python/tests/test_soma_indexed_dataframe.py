import pyarrow as pa
import pytest

import tiledbsoma as t

@pytest.fixture
def arrow_schema():
    def _schema():
        return pa.schema([
            pa.field("foo", pa.int64()),
            pa.field("bar", pa.float64()),
            pa.field("baz", pa.string())
        ])
    return _schema

def test_soma_indexed_dataframe(tmp_path, arrow_schema):
    sdf = t.SOMAIndexedDataFrame(uri=tmp_path.as_posix())

    # Create
    asch = arrow_schema()
    sdf.create(schema=asch, index_column_names=["foo"])

    # Write
    for _ in range(3):
        pydict = {}
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.Table.from_pydict(pydict)
        sdf.write(rb)

    # Read all
    table = sdf.read_all()
    # Weird thing about pyarrow Table:
    # * We have table.num_rows is 5 and table.num_columns is 3
    # * But len(table) is 3
    # * `for column in table` loops over columns
    assert table.num_rows == 5
    assert table.num_columns == 3
    assert [e.as_py() for e in list(table["foo"])] == pydict["foo"]
    assert [e.as_py() for e in list(table["bar"])] == pydict["bar"]
    assert [e.as_py() for e in list(table["baz"])] == pydict["baz"]

    # Read ids
    table = sdf.read_all(ids=[30, 10])
    assert table.num_rows == 2
    assert table.num_columns == 3
    assert sorted([e.as_py() for e in list(table["foo"])]) == [10, 30]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [4.1, 6.3]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["apple", "cat"]

def test_soma_indexed_dataframe_with_float_dim(tmp_path, arrow_schema):
    sdf = t.SOMAIndexedDataFrame(uri=tmp_path.as_posix())
    asch = arrow_schema()
    sdf.create(schema=asch, index_column_names=["bar"])
    assert sdf.get_index_column_names() == ["bar"]

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
    rb = pa.Table.from_pydict(data)
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
        sdf.read_all(ids=ids, column_names=col_names),
        col_names,
        ids,
    )
    _check_tbl(
        sdf.read_all(column_names=col_names),
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
