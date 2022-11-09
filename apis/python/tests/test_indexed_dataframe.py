from typing import List

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


def test_indexed_dataframe(tmp_path, arrow_schema):
    sidf = soma.IndexedDataFrame(uri=tmp_path.as_posix())

    asch = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )

    # Create
    asch = arrow_schema()
    sidf.create(schema=asch, index_column_names=["foo"])

    assert sorted(sidf.schema.names) == sorted(["foo", "bar", "baz", "soma_joinid"])
    assert sorted(sidf.keys()) == sorted(sidf.schema.names)

    # Write
    for _ in range(3):
        pydict = {}
        pydict["soma_joinid"] = [0, 1, 2, 3, 4]
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.Table.from_pydict(pydict)
        sidf.write(rb)

    # Read all
    table = sidf.read_all()
    assert table.num_rows == 5
    assert table.num_columns == 4
    assert [e.as_py() for e in list(table["soma_joinid"])] == pydict["soma_joinid"]
    assert [e.as_py() for e in list(table["foo"])] == pydict["foo"]
    assert [e.as_py() for e in list(table["bar"])] == pydict["bar"]
    assert [e.as_py() for e in list(table["baz"])] == pydict["baz"]

    # Read ids
    table = sidf.read_all(ids=[[30, 10]])
    assert table.num_rows == 2
    assert table.num_columns == 4
    assert sorted([e.as_py() for e in list(table["soma_joinid"])]) == [0, 2]
    assert sorted([e.as_py() for e in list(table["foo"])]) == [10, 30]
    assert sorted([e.as_py() for e in list(table["bar"])]) == [4.1, 6.3]
    assert sorted([e.as_py() for e in list(table["baz"])]) == ["apple", "cat"]


def test_indexed_dataframe_with_float_dim(tmp_path, arrow_schema):
    sidf = soma.IndexedDataFrame(uri=tmp_path.as_posix())
    asch = arrow_schema()
    sidf.create(schema=asch, index_column_names=("bar",))
    assert sidf.get_index_column_names() == ("bar",)


@pytest.fixture
def simple_indexed_data_frame(tmp_path):
    """
    A pytest fixture which creates a simple IndexedDataFrame for use in tests below.
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
    sidf = soma.IndexedDataFrame(uri=tmp_path.as_posix())
    sidf.create(schema=schema, index_column_names=index_column_names)

    data = {
        "index": [0, 1, 2, 3],
        "soma_joinid": [10, 11, 12, 13],
        "A": [10, 11, 12, 13],
        "B": [100.1, 200.2, 300.3, 400.4],
        "C": ["this", "is", "a", "test"],
    }
    n_data = len(data["index"])
    rb = pa.Table.from_pydict(data)
    sidf.write(rb)
    return (schema, sidf, n_data, index_column_names)


@pytest.mark.parametrize(
    "ids",
    [
        [None],
        [0],
        [[1, 3]],
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
def test_IndexedDataFrame_read_column_names(simple_indexed_data_frame, ids, col_names):
    schema, sidf, n_data, index_column_names = simple_indexed_data_frame
    assert sidf.exists()

    def _check_tbl(tbl, col_names, ids, *, demote):
        assert tbl.num_columns == (
            len(schema.names) if col_names is None else len(col_names)
        )

        if ids is None:
            assert tbl.num_rows == n_data
        elif ids[0] is None:
            assert tbl.num_rows == n_data
        elif isinstance(ids[0], int):
            assert tbl.num_rows == 1
        else:
            assert tbl.num_rows == len(ids[0])

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
        sidf.read_all(ids=ids, column_names=col_names), col_names, ids, demote=False
    )
    _check_tbl(sidf.read_all(column_names=col_names), col_names, None, demote=False)

    # TileDB ASCII -> Pandas string -> Arrow string (not large_string)
    _check_tbl(
        pa.Table.from_pandas(
            pd.concat(sidf.read_as_pandas(ids=ids, column_names=col_names))
        ),
        col_names,
        ids,
        demote=True,
    )
    _check_tbl(
        pa.Table.from_pandas(sidf.read_as_pandas_all(column_names=col_names)),
        col_names,
        None,
        demote=True,
    )


def test_empty_indexed_dataframe(tmp_path):
    a = soma.IndexedDataFrame((tmp_path / "A").as_posix())
    a.create(pa.schema([("a", pa.int32())]), index_column_names=["a"])
    # Must not throw
    assert len(next(a.read())) == 0
    assert len(a.read_all()) == 0
    assert len(next(a.read_as_pandas())) == 0
    assert len(a.read_as_pandas_all()) == 0
    assert isinstance(a.read_as_pandas_all(), pd.DataFrame)


def test_columns(tmp_path):
    """
    1. soma_joinid is int64
    2. soma_joinid will be added by default, if missing in call to create
    3. soma_joinid is explicit in keys/schema
    4. No other soma_ ids allowed
    """

    A = soma.IndexedDataFrame((tmp_path / "A").as_posix())
    A.create(pa.schema([("a", pa.int32())]), index_column_names=["a"])
    assert sorted(A.keys()) == sorted(["a", "soma_joinid"])
    assert A.schema.field("soma_joinid").type == pa.int64()
    A.delete()

    B = soma.IndexedDataFrame((tmp_path / "B").as_posix())
    with pytest.raises(ValueError):
        B.create(
            pa.schema([("a", pa.int32()), ("soma_joinid", pa.float32())]),
            index_column_names=["a"],
        )

    D = soma.IndexedDataFrame((tmp_path / "D").as_posix())
    D.create(
        pa.schema([("a", pa.int32()), ("soma_joinid", pa.int64())]),
        index_column_names=["a"],
    )
    assert sorted(D.keys()) == sorted(["a", "soma_joinid"])
    assert D.schema.field("soma_joinid").type == pa.int64()
    D.delete()

    E = soma.IndexedDataFrame((tmp_path / "E").as_posix())
    with pytest.raises(ValueError):
        E.create(
            pa.schema([("a", pa.int32()), ("soma_rowid", pa.bool_())]),
            index_column_names=["a"],
        )


@pytest.fixture
def make_dataframe(request):
    index_type = request.param
    print()
    print("================================================================")
    print("INDEX_TYPE", index_type)

    foo = {
        pa.string(): ["A", "B", "C"],
        pa.large_string(): ["A", "B", "C"],
        pa.binary(): [b"A", b"B", b"C"],
        pa.large_binary(): [b"A", b"B", b"C"],
        **{
            t: np.arange(3, dtype=t.to_pandas_dtype())
            for t in (
                pa.int8(),
                pa.uint8(),
                pa.int16(),
                pa.uint16(),
                pa.int32(),
                pa.uint32(),
                pa.int64(),
                pa.uint64(),
                pa.float32(),
                pa.float64(),
            )
        },
    }
    print("FOO")
    print(foo)

    # {
    #    DataType(string): ['A', 'B', 'C'],
    #    DataType(large_string): ['A', 'B', 'C'],
    #    DataType(binary): [b'A', b'B', b'C'],
    #    DataType(large_binary): [b'A', b'B', b'C'],
    #    DataType(int8): array([0, 1, 2], dtype=int8),
    #    DataType(uint8): array([0, 1, 2], dtype=uint8),
    #    DataType(int16): array([0, 1, 2], dtype=int16),
    #    DataType(uint16): array([0, 1, 2], dtype=uint16),
    #    DataType(int32): array([0, 1, 2], dtype=int32),
    #    DataType(uint32): array([0, 1, 2], dtype=uint32),
    #    DataType(int64): array([0, 1, 2]),
    #    DataType(uint64): array([0, 1, 2], dtype=uint64),
    #    DataType(float): array([0., 1., 2.], dtype=float32),
    #    DataType(double): array([0., 1., 2.])
    # }

    print("BAR")
    print(foo[index_type])

    print("================================================================")

    index = {
        pa.string(): ["A", "B", "C"],
        pa.large_string(): ["A", "B", "C"],
        pa.binary(): [b"A", b"B", b"C"],
        pa.large_binary(): [b"A", b"B", b"C"],
        **{
            t: np.arange(3, dtype=t.to_pandas_dtype())
            for t in (
                pa.int8(),
                pa.uint8(),
                pa.int16(),
                pa.uint16(),
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
        pytest.param(pa.float32(), marks=pytest.mark.xfail),
        pytest.param(pa.float64(), marks=pytest.mark.xfail),
        pytest.param(
            pa.int8(), marks=pytest.mark.xfail
        ),  # TODO: remove xfail when #518 is fixed
        pytest.param(
            pa.uint8(), marks=pytest.mark.xfail
        ),  # TODO: remove xfail when #518 is fixed
        pa.int16(),
        pa.uint16(),
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
def test_index_types(tmp_path, make_dataframe):
    """Verify that the index columns can be of various types"""
    sidf = soma.IndexedDataFrame(tmp_path.as_posix())
    sidf.create(make_dataframe.schema, index_column_names=["index"])
    sidf.write(make_dataframe)


def make_multiply_indexed_dataframe(tmp_path, index_column_names: List[str]):
    """
    Creates a variably-indexed IndexedDataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            # TO DO: Support non-int index types when we have non-int index support
            # in libtiledbsoma's SOMAReader. See also
            # https://github.com/single-cell-data/TileDB-SOMA/issues/418
            # https://github.com/single-cell-data/TileDB-SOMA/issues/419
            ("index1", pa.int64()),
            ("index2", pa.int64()),
            ("index3", pa.int64()),
            ("index4", pa.int64()),
            ("soma_joinid", pa.int64()),
            ("A", pa.int64()),
        ]
    )

    sidf = soma.IndexedDataFrame(uri=tmp_path.as_posix())
    sidf.create(schema=schema, index_column_names=index_column_names)

    data = {
        "index1": [0, 1, 2, 3, 4, 5],
        "index2": [400, 400, 500, 500, 600, 600],
        "index3": [0, 1, 0, 1, 0, 1],
        "index4": [1000, 2000, 1000, 1000, 1000, 1000],
        "soma_joinid": [10, 11, 12, 13, 14, 15],
        "A": [10, 11, 12, 13, 14, 15],
    }

    n_data = len(data["index1"])
    rb = pa.Table.from_pydict(data)
    sidf.write(rb)
    return (schema, sidf, n_data)


@pytest.mark.parametrize(
    "io",
    [
        # 1D: indexing list is None
        {
            "index_column_names": ["index1"],
            "ids": None,
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        # 1D: indexing slot is None
        {
            "index_column_names": ["index1"],
            "ids": [None],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        # 1D: indexing slot is int
        {
            "index_column_names": ["index1"],
            "ids": [0],
            "A": [10],
            "throws": None,
        },
        {
            "index_column_names": ["index1"],
            "ids": [100],
            "A": [],
            "throws": None,
        },
        {
            "index_column_names": ["index1"],
            "ids": [-100],
            "A": [],
            "throws": None,
        },
        # 1D: indexing slot is list
        {
            "index_column_names": ["index1"],
            "ids": [[1, 3]],
            "A": [11, 13],
            "throws": None,
        },
        {
            "index_column_names": ["index1"],
            "ids": [[-100, 100]],
            "A": [],
            "throws": None,
        },
        pytest.param(
            # TODO: use after https://github.com/single-cell-data/TileDB-SOMA/issues/484 is resolved
            {
                "index_column_names": ["index1"],
                "ids": [[]],
                "A": [],
                "throws": None,
            },
            marks=pytest.mark.xfail,
        ),
        # 1D: indexing slot is tuple
        {
            "index_column_names": ["index1"],
            "ids": [(1, 3)],
            "A": [11, 13],
            "throws": None,
        },
        # 1D: indexing slot is range
        {
            "index_column_names": ["index1"],
            "ids": [range(1, 3)],
            "A": [11, 12],
            "throws": None,
        },
        # 1D: indexing slot is pa.ChunkedArray
        {
            "index_column_names": ["index1"],
            "ids": [pa.chunked_array(pa.array([1, 3]))],
            "A": [11, 13],
            "throws": None,
        },
        # 1D: indexing slot is pa.Array
        {
            "index_column_names": ["index1"],
            "ids": [pa.array([1, 3])],
            "A": [11, 13],
            "throws": None,
        },
        # 1D: indexing slot is pa.Array
        {
            "index_column_names": ["index1"],
            "ids": [pa.array([1, 3])],
            "A": [11, 13],
            "throws": None,
        },
        # 1D: indexing slot is np.ndarray
        {
            "index_column_names": ["index1"],
            "ids": [np.asarray([1, 3])],
            "A": [11, 13],
            "throws": None,
        },
        {
            "index_column_names": ["index1"],
            "ids": [np.asarray([[1, 3], [2, 4]])],  # Error since 2D array in the slot
            "A": [11, 13],
            "throws": ValueError,
        },
        # 1D: indexing slot is slice
        {
            "index_column_names": ["index1"],
            "ids": [
                slice(None)
            ],  # Indexing slot is none-slice i.e. `[:]` which is like None
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "index_column_names": ["index1"],
            "ids": [slice(1, 3)],  # Indexing slot is double-ended slice
            "A": [11, 12, 13],
            "throws": None,
        },
        {
            "index_column_names": ["index1"],
            "ids": [slice(None, None)],  # Indexing slot is slice-all
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "index_column_names": ["index1"],
            "ids": [slice(None, 3)],  # Half-slices are not supported yet
            "A": None,
            "throws": ValueError,
        },
        {
            "index_column_names": ["index1"],
            "ids": [slice(1, None)],  # Half-slices are not supported yet
            "A": None,
            "throws": ValueError,
        },
        {
            "index_column_names": ["index1"],
            "ids": [slice(1, 5, 2)],  # Slice step must be 1 or None
            "A": None,
            "throws": ValueError,
        },
        {
            "index_column_names": ["index1"],
            "ids": [slice(-2, -1)],  # Negative slices are not supported
            "A": None,
            "throws": ValueError,
        },
        # 1D: indexing slot is of invalid type
        # TODO: I want to test this but Typeguard fails the test since it already knows strings are not
        # valid until we implement
        # https://github.com/single-cell-data/TileDB-SOMA/issues/418
        # https://github.com/single-cell-data/TileDB-SOMA/issues/419
        pytest.param(
            {
                "index_column_names": ["index1"],
                "ids": ["nonesuch"],  # noqa
                "A": None,
                "throws": soma.SOMAError,
            },
            marks=pytest.mark.xfail,
        ),
        # 2D: indexing list is None
        {
            "index_column_names": ["index2", "index3"],
            "ids": None,
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        # 2D: indexing slot is None
        {
            "index_column_names": ["index2", "index3"],
            "ids": [None, None],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        # 2D: indexing slot is int
        {
            "index_column_names": ["index2", "index3"],
            "ids": [400, 0],
            "A": [10],
            "throws": None,
        },
        # 2D: indexing slot is list
        # TODO: at present SOMAReader only accepts int dims. See also:
        # https://github.com/single-cell-data/TileDB-SOMA/issues/418
        # https://github.com/single-cell-data/TileDB-SOMA/issues/419
        {
            "index_column_names": ["index2", "index3"],
            "ids": [[400, 600], None],
            "A": [10, 11, 14, 15],
            "throws": None,
        },
        # 3D: indexing slot is list
        {
            "index_column_names": ["index2", "index3", "index4"],
            "ids": [[400, 600], None, None],
            "A": [10, 11, 14, 15],
            "throws": None,
        },
        # 3D: indexing slot is mixed
        {
            "index_column_names": ["index2", "index3", "index4"],
            "ids": [range(400, 600), None, np.asarray([2000, 9999])],
            "A": [11],
            "throws": None,
        },
    ],
)
def test_read_indexing(tmp_path, io):
    """Test various ways of indexing on read"""

    schema, sidf, n_data = make_multiply_indexed_dataframe(
        tmp_path, io["index_column_names"]
    )
    assert sidf.exists()

    col_names = ["A"]

    if io["throws"] is not None:
        with pytest.raises(io["throws"]):
            next(sidf.read(ids=io["ids"], column_names=col_names))
    else:
        table = next(sidf.read(ids=io["ids"], column_names=col_names))
        assert table["A"].to_pylist() == io["A"]

    if io["throws"] is not None:
        with pytest.raises(io["throws"]):
            next(sidf.read_as_pandas(ids=io["ids"], column_names=col_names))
    else:
        table = next(sidf.read_as_pandas(ids=io["ids"], column_names=col_names))
        assert table["A"].to_list() == io["A"]

    sidf.delete()
