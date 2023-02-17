from typing import Dict, List

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
import somacore
import tiledb

import tiledbsoma as soma
from tiledbsoma.options import SOMATileDBContext


@pytest.fixture
def arrow_schema():
    def _schema():
        return pa.schema(
            [
                pa.field("foo", pa.int64()),
                pa.field("bar", pa.float64()),
                pa.field("baz", pa.string()),
                pa.field("quux", pa.bool_()),
            ]
        )

    return _schema


def test_dataframe(tmp_path, arrow_schema):
    uri = tmp_path.as_posix()

    asch = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
            ("quux", pa.bool_()),
        ]
    )

    # Create
    asch = arrow_schema()
    with pytest.raises(ValueError):
        # requires one or more index columns
        soma.DataFrame.create(uri, schema=asch, index_column_names=[])
    with pytest.raises(TypeError):
        # invalid schema type
        soma.DataFrame.create(uri, schema=asch.to_string(), index_column_names=[])
    with pytest.raises(ValueError):
        # nonexistent indexed column
        soma.DataFrame.create(uri, schema=asch, index_column_names=["bogus"])
    soma.DataFrame.create(uri, schema=asch, index_column_names=["foo"]).close()

    assert soma.DataFrame.exists(uri)
    assert not soma.Collection.exists(uri)
    assert not soma.SparseNDArray.exists(uri)

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 0
        assert len(sdf) == 0

        assert sorted(sdf.schema.names) == sorted(
            ["foo", "bar", "baz", "soma_joinid", "quux"]
        )
        assert sorted(sdf.keys()) == sorted(sdf.schema.names)

    with soma.DataFrame.open(uri, "w") as sdf:
        # Write
        for _ in range(3):
            pydict = {}
            pydict["soma_joinid"] = [0, 1, 2, 3, 4]
            pydict["foo"] = [10, 20, 30, 40, 50]
            pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
            pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
            pydict["quux"] = [True, False, False, True, False]
            rb = pa.Table.from_pydict(pydict)

            sdf.write(rb)

            with pytest.raises(TypeError):
                # non-arrow write
                sdf.write(rb.to_pandas)

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 5
        assert len(sdf) == 5

        # Read all
        table = sdf.read().concat()
        assert table.num_rows == 5
        assert table.num_columns == 5
        assert [e.as_py() for e in list(table["soma_joinid"])] == pydict["soma_joinid"]
        assert [e.as_py() for e in list(table["foo"])] == pydict["foo"]
        assert [e.as_py() for e in list(table["bar"])] == pydict["bar"]
        assert [e.as_py() for e in list(table["baz"])] == pydict["baz"]
        assert [e.as_py() for e in list(table["quux"])] == pydict["quux"]

        # Read ids
        table = sdf.read(coords=[[30, 10]]).concat()
        assert table.num_rows == 2
        assert table.num_columns == 5
        assert sorted([e.as_py() for e in list(table["soma_joinid"])]) == [0, 2]
        assert sorted([e.as_py() for e in list(table["foo"])]) == [10, 30]
        assert sorted([e.as_py() for e in list(table["bar"])]) == [4.1, 6.3]
        assert sorted([e.as_py() for e in list(table["baz"])]) == ["apple", "cat"]
        assert [e.as_py() for e in list(table["quux"])] == [True, False]

    # Validate TileDB array schema
    with tiledb.open(uri) as A:
        assert A.schema.sparse
        assert not A.schema.allows_duplicates


@pytest.mark.parametrize("appendable", [True, False])
def test_dataframe_append(tmp_path, appendable):
    uri = tmp_path.as_posix()

    soma.DataFrame.create(
        uri,
        schema=pa.schema([("foo", pa.large_string())]),
        shape=None if appendable else (5,),
    ).close()
    assert soma.DataFrame.exists(uri)

    batch1 = pa.Table.from_pydict(
        {
            "soma_joinid": [0, 1, 2, 3, 4],
            "foo": ["apple", "ball", "cat", "dog", "egg"],
        }
    )

    batch2 = pa.Table.from_pydict(
        {
            "soma_joinid": [5, 6, 7],
            "foo": ["fish", "gate", "house"],
        }
    )

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 0
        assert len(sdf) == 0

    with soma.DataFrame.open(uri, "w") as sdf:
        sdf.write(batch1)

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 5
        assert len(sdf) == 5

    if appendable:
        with soma.DataFrame.open(uri, "w") as sdf:
            sdf.write(batch2)
    else:
        # tiledb.cc.TileDBError: [TileDB::Dimension] Error: Coordinate 6
        # is out of domain bounds [0, 4] on dimension 'soma_joinid'
        with pytest.raises(tiledb.cc.TileDBError):
            with soma.DataFrame.open(uri, "w") as sdf:
                sdf.write(batch2)


def test_dataframe_with_float_dim(tmp_path, arrow_schema):
    sdf = soma.DataFrame.create(
        tmp_path.as_posix(), schema=arrow_schema(), index_column_names=("bar",)
    )
    assert sdf.index_column_names == ("bar",)


@pytest.fixture
def simple_data_frame(tmp_path):
    """
    A pytest fixture which creates a simple DataFrame for use in tests below.
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
    with soma.DataFrame.create(
        tmp_path.as_posix(), schema=schema, index_column_names=index_column_names
    ) as sdf:

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
    sdf = soma.DataFrame.open(tmp_path.as_posix())
    return (schema, sdf, n_data, index_column_names)


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
def test_DataFrame_read_column_names(simple_data_frame, ids, col_names):
    schema, sdf, n_data, index_column_names = simple_data_frame

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
        sdf.read(ids, column_names=col_names).concat(),
        col_names,
        ids,
        demote=False,
    )
    _check_tbl(sdf.read(column_names=col_names).concat(), col_names, None, demote=False)

    # TileDB ASCII -> Pandas string -> Arrow string (not large_string)
    _check_tbl(
        pa.Table.from_pandas(
            pd.concat(
                [tbl.to_pandas() for tbl in sdf.read(ids, column_names=col_names)]
            )
        ),
        col_names,
        ids,
        demote=True,
    )
    _check_tbl(
        pa.Table.from_pandas(sdf.read(column_names=col_names).concat().to_pandas()),
        col_names,
        None,
        demote=True,
    )


def test_empty_dataframe(tmp_path):
    soma.DataFrame.create(
        (tmp_path / "A").as_posix(),
        schema=pa.schema([("a", pa.int32())]),
        index_column_names=["a"],
    ).close()
    with soma.DataFrame.open((tmp_path / "A").as_posix()) as a:
        # Must not throw
        assert len(next(a.read())) == 0
        assert len(a.read().concat()) == 0
        assert len(next(a.read()).to_pandas()) == 0
        assert len(a.read().concat().to_pandas()) == 0
        assert isinstance(a.read().concat().to_pandas(), pd.DataFrame)

    with pytest.raises(ValueError):
        # illegal column name
        soma.DataFrame.create(
            (tmp_path / "B").as_posix(),
            schema=pa.schema([("a", pa.int32()), ("soma_bogus", pa.int32())]),
            index_column_names=["a"],
        )


def test_columns(tmp_path):
    """
    1. soma_joinid is int64
    2. soma_joinid will be added by default, if missing in call to create
    3. soma_joinid is explicit in keys/schema
    4. No other soma_ ids allowed
    """

    A = soma.DataFrame.create(
        (tmp_path / "A").as_posix(),
        schema=pa.schema([("a", pa.int32())]),
        index_column_names=["a"],
    )
    assert sorted(A.keys()) == sorted(["a", "soma_joinid"])
    assert A.schema.field("soma_joinid").type == pa.int64()

    with pytest.raises(ValueError):
        soma.DataFrame.create(
            (tmp_path / "B").as_posix(),
            schema=pa.schema([("a", pa.int32()), ("soma_joinid", pa.float32())]),
            index_column_names=["a"],
        )

    D = soma.DataFrame.create(
        (tmp_path / "D").as_posix(),
        schema=pa.schema([("a", pa.int32()), ("soma_joinid", pa.int64())]),
        index_column_names=["a"],
    )
    assert sorted(D.keys()) == sorted(["a", "soma_joinid"])
    assert D.schema.field("soma_joinid").type == pa.int64()

    with pytest.raises(ValueError):
        soma.DataFrame.create(
            (tmp_path / "E").as_posix(),
            schema=pa.schema(
                [("a", pa.int32()), ("soma_is_a_reserved_prefix", pa.bool_())]
            ),
            index_column_names=["a"],
        )


@pytest.fixture
def make_dataframe(request):
    index_type = request.param

    # TODO: https://github.com/single-cell-data/TileDB-SOMA/issues/518
    # Check against all `SUPPORTED_ARROW_TYPES` in tests/test_type_system.py`

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
        pa.float32(),
        pa.float64(),
        pa.int8(),
        pa.uint8(),
        pa.int16(),
        pa.uint16(),
        pa.int32(),
        pa.uint32(),
        pa.int64(),
        pa.uint64(),
        pa.string(),
        pa.large_string(),
        pa.binary(),
        pa.large_binary(),
    ],
    indirect=True,
)
def test_index_types(tmp_path, make_dataframe):
    """Verify that the index columns can be of various types"""
    sdf = soma.DataFrame.create(
        tmp_path.as_posix(), schema=make_dataframe.schema, index_column_names=["index"]
    )
    sdf.write(make_dataframe)


def make_multiply_indexed_dataframe(tmp_path, index_column_names: List[str]):
    """
    Creates a variably-indexed DataFrame for use in tests below.
    """
    schema = pa.schema(
        [
            # TO DO: Support other index types when we have support for more than int and string/bytes
            # index types in libtiledbsoma's SOMAReader.
            ("0_thru_5", pa.int64()),
            ("strings_aaa", pa.string()),
            ("zero_one", pa.int64()),
            ("thousands", pa.int64()),
            ("both_signs", pa.int64()),
            ("soma_joinid", pa.int64()),
            ("A", pa.int64()),
        ]
    )

    sdf = soma.DataFrame.create(
        uri=tmp_path.as_posix(), schema=schema, index_column_names=index_column_names
    )

    data: Dict[str, list] = {
        "0_thru_5": [0, 1, 2, 3, 4, 5],
        "strings_aaa": ["aaa", "aaa", "bbb", "bbb", "ccc", "ccc"],
        "zero_one": [0, 1, 0, 1, 0, 1],
        "thousands": [1000, 2000, 1000, 1000, 1000, 1000],
        "both_signs": [-1, -2, -3, 1, 2, 3],
        "soma_joinid": [10, 11, 12, 13, 14, 15],
        "A": [10, 11, 12, 13, 14, 15],
    }

    n_data = len(data["0_thru_5"])
    sdf.write(pa.Table.from_pandas(pd.DataFrame(data=data)))

    return (schema, sdf, n_data)


@pytest.mark.parametrize(
    "io",
    [
        {
            "name": "1D indexing slot is None",
            "index_column_names": ["0_thru_5"],
            "coords": [None],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing slot is int",
            "index_column_names": ["0_thru_5"],
            "coords": [0],
            "A": [10],
            "throws": None,
        },
        {
            "name": "1D no results for 100",
            "index_column_names": ["0_thru_5"],
            "coords": [100],
            "A": [],
            "throws": None,
        },
        {
            "name": "1D no results for -100",
            "index_column_names": ["0_thru_5"],
            "coords": [-100],
            "A": [],
            "throws": None,
        },
        {
            "name": "1D indexing slot is list",
            "index_column_names": ["0_thru_5"],
            "coords": [[1, 3]],
            "A": [11, 13],
            "throws": None,
        },
        {
            "name": "1D no results for -100, 100",
            "index_column_names": ["0_thru_5"],
            "coords": [[-100, 100]],
            "A": [],
            "throws": None,
        },
        {
            "name": "1D empty list returns empty results",
            "index_column_names": ["0_thru_5"],
            "coords": [[]],
            "A": [],
            "throws": None,
        },
        {
            "name": "1D indexing slot is tuple",
            "index_column_names": ["0_thru_5"],
            "coords": [(1, 3)],
            "A": [11, 13],
            "throws": None,
        },
        {
            "name": "1D indexing slot is range",
            "index_column_names": ["0_thru_5"],
            "coords": [range(1, 3)],
            "A": [11, 12],
            "throws": None,
        },
        {
            "name": "1D indexing slot is pa.ChunkedArray",
            "index_column_names": ["0_thru_5"],
            "coords": [pa.chunked_array(pa.array([1, 3]))],
            "A": [11, 13],
            "throws": None,
        },
        {
            "name": "1D indexing slot is pa.Array",
            "index_column_names": ["0_thru_5"],
            "coords": [pa.array([1, 3])],
            "A": [11, 13],
            "throws": None,
        },
        # 1D: indexing slot is np.ndarray
        {
            "name": "1D indexing slot is np.ndarray",
            "index_column_names": ["0_thru_5"],
            "coords": [np.asarray([1, 3])],
            "A": [11, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by 2D np.ndarray",
            "index_column_names": ["0_thru_5"],
            "coords": [
                np.asarray([[1, 3], [2, 4]])
            ],  # Error since 2D array in the slot
            "A": [11, 13],
            "throws": ValueError,
        },
        {
            "name": "1D indexing by slice(None)",
            "index_column_names": ["0_thru_5"],
            "coords": [
                slice(None)
            ],  # Indexing slot is none-slice i.e. `[:]` which is like None
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing by empty coords",
            "index_column_names": ["0_thru_5"],
            "coords": [],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing by 1:3",
            "index_column_names": ["0_thru_5"],
            "coords": [slice(1, 3)],  # Indexing slot is double-ended slice
            "A": [11, 12, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by [:3]",
            "index_column_names": ["0_thru_5"],
            "coords": [slice(None, 3)],  # Half-slice
            "A": [10, 11, 12, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by [2:]",
            "index_column_names": ["0_thru_5"],
            "coords": [slice(2, None)],  # Half-slice
            "A": [12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing with negatives",
            "index_column_names": ["both_signs"],
            "coords": [slice(-2, 1)],
            "A": [11, 10, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by ['bbb':'c']",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": [slice("bbb", "c")],
            "A": [12, 13],
            "throws": None,
        },
        {
            "name": "1D indexing by ['ccc':]",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": [slice("ccc", None)],
            "A": [14, 15],
            "throws": None,
        },
        {
            "name": "1D indexing by [:'bbd']",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": [slice("bbd")],
            "A": [10, 11, 12, 13],
            "throws": None,
        },
        {
            "name": "1D indexing with one partition",
            "index_column_names": ["0_thru_5"],
            "coords": [slice(2, None)],
            "partitions": somacore.IOfN(0, 1),
            "A": [12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "partitioned reads unimplemented",
            "index_column_names": ["0_thru_5"],
            "coords": [],
            "partitions": somacore.IOfN(1, 2),
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "steps forbidden",
            "index_column_names": ["0_thru_5"],
            "coords": [slice(1, 5, 2)],
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "slice must overlap domain (negative)",
            "index_column_names": ["soma_joinid"],
            "coords": [slice(-2, -1)],
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "backwards slice",
            "index_column_names": ["0_thru_5"],
            "coords": [slice(1, 0)],
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "too many columns",
            "index_column_names": ["0_thru_5"],
            "coords": [(1,), (2,)],
            "A": None,
            "throws": ValueError,
        },
        {
            "name": "wrong coords type",
            "index_column_names": ["0_thru_5"],
            "coords": "bogus",
            "A": None,
            "throws": TypeError,
        },
        {
            "name": "bad index type dict",
            "index_column_names": ["0_thru_5"],
            "coords": [{"bogus": True}],
            "A": None,
            "throws": TypeError,
        },
        {
            "name": "bad index type bool",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": [[True], slice(None)],
            "A": None,
            "throws": (RuntimeError, tiledb.cc.TileDBError),
        },
        {
            "name": "2D index empty",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": (),
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "2D index None",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": [None, None],
            "A": [10, 11, 12, 13, 14, 15],
            "throws": None,
        },
        {
            "name": "2D index 0, 0",
            "index_column_names": ["0_thru_5", "zero_one"],
            "coords": [0, 0],
            "A": [10],
            "throws": None,
        },
        {
            "name": "2D index str, int",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": [["aaa"], 0],
            "A": [10],
            "throws": None,
        },
        {
            "name": "2D index str, not sequence[str]",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": ["aaa", 0],
            "A": [10],
            "throws": None,
        },
        {
            "name": "2D index List[str]",
            "index_column_names": ["strings_aaa", "zero_one"],
            "coords": [["aaa", "ccc"], None],
            "A": [10, 11, 14, 15],
            "throws": None,
        },
        {
            "name": "3D index List[str]",
            "index_column_names": ["strings_aaa", "zero_one", "thousands"],
            "coords": [["aaa", "ccc"], None, None],
            "A": [10, 11, 14, 15],
            "throws": None,
        },
        {
            "name": "3D index mixed",
            "index_column_names": ["strings_aaa", "zero_one", "thousands"],
            "coords": [("aaa", "ccc"), None, np.asarray([2000, 9999])],
            "A": [11],
            "throws": None,
        },
        {
            "name": "value filter good",
            "index_column_names": ["0_thru_5", "strings_aaa"],
            "coords": [None, ("ccc", "zzz")],
            "value_filter": "soma_joinid > 13",
            "A": [14, 15],
        },
        {
            "name": "value filter bad",
            "index_column_names": ["0_thru_5", "strings_aaa"],
            "coords": [None, ("bbb", "zzz")],
            "value_filter": "quick brown fox",
            "A": None,
            "throws": soma.SOMAError,
        },
    ],
    ids=lambda d: d.get("name"),
)
def test_read_indexing(tmp_path, io):
    """Test various ways of indexing on read"""

    schema, sdf, n_data = make_multiply_indexed_dataframe(
        tmp_path, io["index_column_names"]
    )
    with soma.DataFrame.open(uri=sdf.uri) as sdf:
        assert list(sdf.index_column_names) == io["index_column_names"]

        read_kwargs = {"column_names": ["A"]}
        read_kwargs.update(
            {k: io[k] for k in ("coords", "partitions", "value_filter") if k in io}
        )
        if io.get("throws", None):
            with pytest.raises(io["throws"]):
                next(sdf.read(**read_kwargs))
        else:
            table = next(sdf.read(**read_kwargs))
            assert table["A"].to_pylist() == io["A"]

        if io.get("throws", None):
            with pytest.raises(io["throws"]):
                next(sdf.read(**read_kwargs)).to_pandas()
        else:
            table = next(sdf.read(**read_kwargs)).to_pandas()
            assert table["A"].to_list() == io["A"]


@pytest.mark.parametrize(
    "schema",
    [
        pa.schema(
            [
                (
                    "A",
                    pa.dictionary(
                        value_type=pa.string(), index_type=pa.int8(), ordered=True
                    ),
                ),
            ]
        ),
        pa.schema(
            [
                (
                    "A",
                    pa.dictionary(
                        value_type=pa.string(), index_type=pa.int8(), ordered=False
                    ),
                ),
            ]
        ),
        pa.Schema.from_pandas(
            pd.DataFrame(
                data={
                    "A": pd.Categorical(
                        ["a", "b", "a", "b"], ordered=True, categories=["b", "a"]
                    )
                }
            )
        ),
        pa.Schema.from_pandas(
            pd.DataFrame(
                data={
                    "A": pd.Categorical(
                        ["a", "b", "a", "b"], ordered=False, categories=["b", "a"]
                    )
                }
            )
        ),
    ],
)
def test_create_categorical_types(tmp_path, schema):
    """
    Verify that `create` throws expected error on (unsupported) dictionary/categorical types.
    """
    schema = schema.insert(0, pa.field("soma_joinid", pa.int64()))

    # Test exception as normal column
    with pytest.raises(TypeError):
        soma.DataFrame.create(
            tmp_path.as_posix(), schema=schema, index_column_names=["soma_joinid"]
        )

    # test as index column
    with pytest.raises(TypeError):
        soma.DataFrame.create(
            tmp_path.as_posix(), schema=schema, index_column_names=["A"]
        )


def test_write_categorical_types(tmp_path):
    """
    Verify that write path accepts categoricals
    """
    schema = pa.schema([("soma_joinid", pa.int64()), ("A", pa.large_string())])
    with soma.DataFrame.create(
        tmp_path.as_posix(), schema=schema, index_column_names=["soma_joinid"]
    ) as sdf:

        df = pd.DataFrame(
            data={
                "soma_joinid": [0, 1, 2, 3],
                "A": pd.Categorical(
                    ["a", "b", "a", "b"], ordered=True, categories=["b", "a"]
                ),
            }
        )
        sdf.write(pa.Table.from_pandas(df))

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        assert (df == sdf.read().concat().to_pandas()).all().all()


def test_result_order(tmp_path):
    # cf. https://docs.tiledb.com/main/background/key-concepts-and-data-format#data-layout
    schema = pa.schema(
        [
            ("row", pa.int64()),
            ("col", pa.int64()),
            ("soma_joinid", pa.int64()),
        ]
    )
    with soma.DataFrame.create(
        uri=tmp_path.as_posix(), schema=schema, index_column_names=["row", "col"]
    ) as sdf:
        data = {
            "row": [0] * 4 + [1] * 4 + [2] * 4 + [3] * 4,
            "col": [0, 1, 2, 3] * 4,
            "soma_joinid": list(range(16)),
        }
        sdf.write(pa.Table.from_pydict(data))

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        table = sdf.read(result_order="row-major").concat().to_pandas()
        assert table["soma_joinid"].to_list() == list(range(16))

        table = sdf.read(result_order="column-major").concat().to_pandas()
        assert table["soma_joinid"].to_list() == [
            0,
            4,
            8,
            12,
            1,
            5,
            9,
            13,
            2,
            6,
            10,
            14,
            3,
            7,
            11,
            15,
        ]

        with pytest.raises(ValueError):
            next(sdf.read(result_order="bogus"))


@pytest.mark.parametrize(
    "create_options,expected_schema_fields",
    (
        (
            {"allows_duplicates": True},
            {
                "validity_filters": tiledb.FilterList([tiledb.RleFilter()]),
                "allows_duplicates": True,
            },
        ),
        (
            {"allows_duplicates": False},
            {
                "validity_filters": tiledb.FilterList([tiledb.RleFilter()]),
                "allows_duplicates": False,
            },
        ),
        (
            {"validity_filters": ["NoOpFilter"], "allows_duplicates": False},
            {
                "validity_filters": tiledb.FilterList([tiledb.NoOpFilter()]),
                "allows_duplicates": False,
            },
        ),
    ),
)
def test_create_platform_config_overrides(
    tmp_path, create_options, expected_schema_fields
):
    uri = tmp_path.as_posix()
    soma.DataFrame.create(
        uri,
        schema=pa.schema([pa.field("colA", pa.string())]),
        platform_config={"tiledb": {"create": {**create_options}}},
    ).close()
    with tiledb.open(uri) as D:
        for k, v in expected_schema_fields.items():
            assert getattr(D.schema, k) == v


def test_timestamped_ops(tmp_path):
    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("B", pa.float64()),
            ("C", pa.large_string()),
        ]
    )
    with soma.DataFrame.create(
        tmp_path.as_posix(),
        schema=schema,
        index_column_names=["soma_joinid"],
        context=SOMATileDBContext(write_timestamp=10),
    ) as sidf:
        data = {
            "soma_joinid": [0],
            "B": [100.1],
            "C": ["foo"],
        }
        sidf.write(pa.Table.from_pydict(data))

    with soma.DataFrame.open(
        uri=tmp_path.as_posix(), mode="w", context=SOMATileDBContext(write_timestamp=20)
    ) as sidf:
        data = {
            "soma_joinid": [0, 1],
            "B": [200.2, 300.3],
            "C": ["bar", "bas"],
        }
        sidf.write(pa.Table.from_pydict(data))

    # read without timestamp & see final image
    with soma.DataFrame.open(tmp_path.as_posix()) as sidf:
        tab = sidf.read().concat()
        assert list(x.as_py() for x in tab["soma_joinid"]) == [0, 1]
        assert list(x.as_py() for x in tab["B"]) == [200.2, 300.3]
        assert list(x.as_py() for x in tab["C"]) == ["bar", "bas"]

    # read at t=15 & see only the first write
    with soma.DataFrame.open(
        tmp_path.as_posix(), context=SOMATileDBContext(read_timestamp=15)
    ) as sidf:
        tab = sidf.read().concat()
        assert list(x.as_py() for x in tab["soma_joinid"]) == [0]
        assert list(x.as_py() for x in tab["B"]) == [100.1]
        assert list(x.as_py() for x in tab["C"]) == ["foo"]
