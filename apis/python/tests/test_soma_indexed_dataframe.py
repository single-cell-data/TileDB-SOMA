import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

import tiledbsoma as soma


@pytest.fixture
def arrow_schema():
    def _schema():
        return pa.schema(
            [
                pa.field("foo", pa.int64()),
                pa.field("bar", pa.float64()),
                pa.field("baz", pa.string()),
            ]
        )

    return _schema


def test_soma_indexed_dataframe(tmp_path, arrow_schema):
    sdf = soma.SOMAIndexedDataFrame(uri=tmp_path.as_posix())

    asch = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )

    # Create
    asch = arrow_schema()
    sdf.create(schema=asch, index_column_names=["foo"])

    assert sorted(sdf.schema.names) == sorted(["foo", "bar", "baz", "soma_joinid"])
    assert sorted(sdf.keys()) == sorted(sdf.schema.names)

    # Write
    for _ in range(3):
        pydict = {}
        pydict["soma_joinid"] = [0, 1, 2, 3, 4]
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
    assert table.num_columns == 4
    assert [e.as_py() for e in list(table["soma_joinid"])] == pydict["soma_joinid"]
    assert [e.as_py() for e in list(table["foo"])] == pydict["foo"]
    assert [e.as_py() for e in list(table["bar"])] == pydict["bar"]
    assert [e.as_py() for e in list(table["baz"])] == pydict["baz"]

    # Read ids
    table = sdf.read_all(ids=[30, 10])
    assert table.num_rows == 2
    assert table.num_columns == 4
    assert sorted([e.as_py() for e in list(table["soma_joinid"])]) == [0, 2]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [10, 30]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [4.1, 6.3]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["apple", "cat"]


def test_soma_indexed_dataframe_with_float_dim(tmp_path, arrow_schema):
    sdf = soma.SOMAIndexedDataFrame(uri=tmp_path.as_posix())
    asch = arrow_schema()
    sdf.create(schema=asch, index_column_names=("bar",))
    assert sdf.get_index_column_names() == ("bar",)


@pytest.fixture
def simple_soma_indexed_data_frame(tmp_path):
    """
    A pytest fixture which creates a simple SOMAIndexedDataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            ("index", pa.int64()),
            ("soma_joinid", pa.int64()),
            ("A", pa.int64()),
            ("B", pa.float64()),
            ("C", pa.large_string()),
        ]
    )
    index_column_names = ["index"]
    sdf = soma.SOMAIndexedDataFrame(uri=tmp_path.as_posix())
    sdf.create(schema=schema, index_column_names=index_column_names)

    data = {
        "index": [0, 1, 2, 3],
        "soma_joinid": [10, 11, 12, 13],
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
        [0],
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
        ["soma_joinid"],
        ["soma_joinid", "A"],
        None,
    ],
)
def test_SOMAIndexedDataFrame_read_column_names(
    simple_soma_indexed_data_frame, ids, col_names
):
    schema, sdf, n_data, index_column_names = simple_soma_indexed_data_frame
    assert sdf.exists()

    def _check_tbl(tbl, col_names, ids, *, demote):
        assert tbl.num_columns == (
            len(schema.names) if col_names is None else len(col_names)
        )
        assert tbl.num_rows == (n_data if ids is None else len(ids))

        if demote:
            assert tbl.schema == pa.schema(
                [
                    pa.field(schema.field(f).name, pa.string())
                    if schema.field(f).type == pa.large_string()
                    else schema.field(f)
                    for f in (col_names if col_names is not None else schema.names)
                ]
            )
        else:
            assert tbl.schema == pa.schema(
                [
                    schema.field(f)
                    for f in (col_names if col_names is not None else schema.names)
                ]
            )

    # TileDB ASCII -> Arrow large_string
    _check_tbl(
        sdf.read_all(ids=ids, column_names=col_names), col_names, ids, demote=False
    )
    _check_tbl(sdf.read_all(column_names=col_names), col_names, None, demote=False)

    # TileDB ASCII -> Pandas string -> Arrow string (not large_string)
    _check_tbl(
        pa.Table.from_pandas(
            pd.concat(sdf.read_as_pandas(ids=ids, column_names=col_names))
        ),
        col_names,
        ids,
        demote=True,
    )
    _check_tbl(
        pa.Table.from_pandas(sdf.read_as_pandas_all(column_names=col_names)),
        col_names,
        None,
        demote=True,
    )


def test_soma_columns(tmp_path):
    """
    1. soma_joinid is int64
    2. soma_joinid will be added by default, if missing in call to create
    3. soma_joinid is explicit in keys/schema
    4. No other soma_ ids allowed
    """

    A = soma.SOMAIndexedDataFrame((tmp_path / "A").as_posix())
    A.create(pa.schema([("a", pa.int32())]), index_column_names=["a"])
    assert sorted(A.keys()) == sorted(["a", "soma_joinid"])
    assert A.schema.field("soma_joinid").type == pa.int64()
    A.delete()

    B = soma.SOMAIndexedDataFrame((tmp_path / "B").as_posix())
    with pytest.raises(ValueError):
        B.create(
            pa.schema([("a", pa.int32()), ("soma_joinid", pa.float32())]),
            index_column_names=["a"],
        )

    D = soma.SOMAIndexedDataFrame((tmp_path / "D").as_posix())
    D.create(
        pa.schema([("a", pa.int32()), ("soma_joinid", pa.int64())]),
        index_column_names=["a"],
    )
    assert sorted(D.keys()) == sorted(["a", "soma_joinid"])
    assert D.schema.field("soma_joinid").type == pa.int64()
    D.delete()

    E = soma.SOMAIndexedDataFrame((tmp_path / "E").as_posix())
    with pytest.raises(ValueError):
        E.create(
            pa.schema([("a", pa.int32()), ("soma_rowid", pa.bool_())]),
            index_column_names=["a"],
        )


@pytest.fixture
def make_dataframe(request):
    index_type = request.param

    index = {
        pa.string(): ["A", "B", "C"],
        pa.large_string(): ["A", "B", "C"],
        pa.binary(): [b"A", b"B", b"C"],
        pa.large_binary(): [b"A", b"B", b"C"],
        **{
            t: np.arange(3, dtype=t.to_pandas_dtype())
            for t in (
                pa.int32(),
                pa.uint32(),
                pa.int64(),
                pa.uint64(),
                pa.float32(),
                pa.float64(),
            )
        },
    }[index_type]

    df = pd.DataFrame(
        data={
            "index": index,
            "soma_joinid": np.arange(3, dtype=np.int64),
            "ascii": ["aa", "bbb", "cccccc"],
            "float32": np.array([0.0, 1.1, 2.2], np.float32),
        }
    )
    return pa.Table.from_pandas(df)


@pytest.mark.parametrize(
    "make_dataframe",
    [
        pa.float32(),
        pa.float64(),
        pa.int32(),
        pa.uint32(),
        pa.int64(),
        pa.uint64(),
        pytest.param(
            pa.string(), marks=pytest.mark.xfail
        ),  # TODO: remove xfail when #418 is fixed
        pytest.param(
            pa.large_string(), marks=pytest.mark.xfail
        ),  # TODO: remove xfail when #418 is fixed
        pytest.param(
            pa.binary(), marks=pytest.mark.xfail
        ),  # TODO: remove xfail when #419 is fixed
        pytest.param(
            pa.large_binary(), marks=pytest.mark.xfail
        ),  # TODO: remove xfail when #419 is fixed
    ],
    indirect=True,
)
def test_soma_index_types(tmp_path, make_dataframe):
    """Verify that the index columns can be of various types"""
    sdf = soma.SOMAIndexedDataFrame(tmp_path.as_posix())
    sdf.create(make_dataframe.schema, index_column_names=["index"])
    sdf.write(make_dataframe)
