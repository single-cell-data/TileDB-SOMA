import contextlib
import datetime
import os
from typing import Dict, List

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
import somacore
from pandas.api.types import union_categoricals

import tiledbsoma as soma
import tiledb

from tests._util import raises_no_typeguard


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
    with raises_no_typeguard(TypeError):
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

            with raises_no_typeguard(TypeError):
                # non-arrow write
                sdf.write(rb.to_pandas)

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 5
        assert len(sdf) == 5

        # More to come on https://github.com/single-cell-data/TileDB-SOMA/issues/2407
        assert (
            sdf.tiledbsoma_has_upgraded_domain
            == soma._flags.NEW_SHAPE_FEATURE_FLAG_ENABLED
        )

        with pytest.raises(AttributeError):
            assert sdf.shape is None

        # soma_joinid is not a dim here
        assert sdf._maybe_soma_joinid_shape is None
        assert sdf._maybe_soma_joinid_maxshape is None

        # Read all
        table = sdf.read().concat()
        assert table.num_rows == 5
        assert table.num_columns == 5
        assert [e.as_py() for e in table["soma_joinid"]] == pydict["soma_joinid"]
        assert [e.as_py() for e in table["foo"]] == pydict["foo"]
        assert [e.as_py() for e in table["bar"]] == pydict["bar"]
        assert [e.as_py() for e in table["baz"]] == pydict["baz"]
        assert [e.as_py() for e in table["quux"]] == pydict["quux"]

        # Read ids
        table = sdf.read(coords=[[30, 10]]).concat()
        assert table.num_rows == 2
        assert table.num_columns == 5
        assert sorted([e.as_py() for e in table["soma_joinid"]]) == [0, 2]
        assert sorted([e.as_py() for e in table["foo"]]) == [10, 30]
        assert sorted([e.as_py() for e in table["bar"]]) == [4.1, 6.3]
        assert sorted([e.as_py() for e in table["baz"]]) == ["apple", "cat"]
        assert [e.as_py() for e in table["quux"]] == [True, False]

    # Open and read with bindings
    with contextlib.closing(
        soma.pytiledbsoma.SOMADataFrame.open(
            uri, soma.pytiledbsoma.OpenMode.read, soma.pytiledbsoma.SOMAContext()
        )
    ) as sdf:
        table = sdf.read_next()
        assert table.num_rows == 5
        assert table.num_columns == 5
        assert [e.as_py() for e in table["soma_joinid"]] == pydict["soma_joinid"]
        assert [e.as_py() for e in table["foo"]] == pydict["foo"]
        assert [e.as_py() for e in table["bar"]] == pydict["bar"]
        assert [e.as_py() for e in table["baz"]] == pydict["baz"]
        assert [e.as_py() for e in table["quux"]] == pydict["quux"]

    # Validate TileDB array schema
    with tiledb.open(uri) as A:
        assert A.schema.sparse
        assert not A.schema.allows_duplicates
        assert A.dim("foo").filters == [tiledb.ZstdFilter(level=3)]
        assert A.attr("bar").filters == [tiledb.ZstdFilter()]

    with soma.DataFrame.open(uri) as sdf:
        assert sdf.count == 5
        assert len(sdf) == 5

    # Ensure read mode uses clib object
    with soma.DataFrame.open(tmp_path.as_posix(), "r") as A:
        assert isinstance(A._handle._handle, soma.pytiledbsoma.SOMADataFrame)

    # Ensure write mode uses clib object
    with soma.DataFrame.open(tmp_path.as_posix(), "w") as A:
        assert isinstance(A._handle._handle, soma.pytiledbsoma.SOMADataFrame)


def test_dataframe_reopen(tmp_path, arrow_schema):
    soma.DataFrame.create(
        tmp_path.as_posix(), schema=arrow_schema(), tiledb_timestamp=1
    )

    with soma.DataFrame.open(tmp_path.as_posix(), "r", tiledb_timestamp=1) as sdf1:
        with raises_no_typeguard(ValueError):
            sdf1.reopen("invalid")

        with sdf1.reopen("w", tiledb_timestamp=2) as sdf2:
            with sdf2.reopen("r", tiledb_timestamp=3) as sdf3:
                assert sdf1.mode == "r"
                assert sdf2.mode == "w"
                assert sdf3.mode == "r"
                assert sdf1.tiledb_timestamp_ms == 1
                assert sdf2.tiledb_timestamp_ms == 2
                assert sdf3.tiledb_timestamp_ms == 3

    ts1 = datetime.datetime(2023, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    ts2 = datetime.datetime(2024, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    with soma.DataFrame.open(tmp_path.as_posix(), "r", tiledb_timestamp=ts1) as sdf1:
        with sdf1.reopen("r", tiledb_timestamp=ts2) as sdf2:
            assert sdf1.mode == "r"
            assert sdf2.mode == "r"
            assert sdf1.tiledb_timestamp == ts1
            assert sdf2.tiledb_timestamp == ts2

    with soma.DataFrame.open(tmp_path.as_posix(), "w") as sdf1:
        with sdf1.reopen("w", tiledb_timestamp=None) as sdf2:
            with sdf1.reopen("w") as sdf3:
                assert sdf1.mode == "w"
                assert sdf2.mode == "w"
                assert sdf3.mode == "w"
                now = datetime.datetime.now(datetime.timezone.utc)
                assert sdf1.tiledb_timestamp <= now
                assert sdf2.tiledb_timestamp <= now
                assert sdf3.tiledb_timestamp <= now


def test_dataframe_with_float_dim(tmp_path, arrow_schema):
    sdf = soma.DataFrame.create(
        tmp_path.as_posix(), schema=arrow_schema(), index_column_names=("bar",)
    )
    assert sdf.index_column_names == ("bar",)


def test_dataframe_with_enumeration(tmp_path):
    schema = pa.schema(
        [
            pa.field("foo", pa.dictionary(pa.int64(), pa.large_string())),
            pa.field("bar", pa.dictionary(pa.int64(), pa.large_string())),
        ]
    )
    enums = {"enmr1": ("a", "bb", "ccc"), "enmr2": ("cat", "dog")}
    with soma.DataFrame.create(tmp_path.as_posix(), schema=schema) as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["foo"] = ["a", "bb", "ccc", "bb", "a"]
        data["bar"] = ["cat", "dog", "cat", "cat", "cat"]
        with pytest.raises(soma.SOMAError):
            sdf.write(pa.Table.from_pydict(data))

        data["foo"] = pd.Categorical(["a", "bb", "ccc", "bb", "a"])
        data["bar"] = pd.Categorical(["cat", "dog", "cat", "cat", "cat"])
        sdf.write(pa.Table.from_pydict(data))

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        df = sdf.read().concat()
        np.testing.assert_array_equal(df["foo"].chunk(0).dictionary, enums["enmr1"])
        np.testing.assert_array_equal(df["bar"].chunk(0).dictionary, enums["enmr2"])


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
                    (
                        pa.field(schema.field(f).name, pa.string())
                        if schema.field(f).type == pa.large_string()
                        else schema.field(f)
                    )
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
            # Note: exhaustive type-coverage is in a separate test case.
            # (As of this writing: test_dataframe_column_indexing.py)
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
            "throws": TypeError,
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


def test_write_categorical_types(tmp_path):
    """
    Verify that write path accepts categoricals
    """
    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            (
                "string-ordered",
                pa.dictionary(pa.int8(), pa.large_string(), ordered=True),
            ),
            ("string-unordered", pa.dictionary(pa.int8(), pa.large_string())),
            ("string-compat", pa.large_string()),
            ("int-ordered", pa.dictionary(pa.int8(), pa.int64(), ordered=True)),
            ("int-unordered", pa.dictionary(pa.int8(), pa.int64())),
            ("int-compat", pa.int64()),
            ("bool-ordered", pa.dictionary(pa.int8(), pa.bool_(), ordered=True)),
            ("bool-unordered", pa.dictionary(pa.int8(), pa.bool_())),
            ("bool-compat", pa.bool_()),
        ]
    )
    with soma.DataFrame.create(
        tmp_path.as_posix(), schema=schema, index_column_names=["soma_joinid"]
    ) as sdf:
        df = pd.DataFrame(
            data={
                "soma_joinid": [0, 1, 2, 3],
                "string-ordered": pd.Categorical(
                    ["a", "b", "a", "b"], ordered=True, categories=["b", "a"]
                ),
                "string-unordered": pd.Categorical(
                    ["a", "b", "a", "b"], ordered=False, categories=["b", "a"]
                ),
                "string-compat": pd.Categorical(
                    ["a", "b", "a", "b"], ordered=False, categories=["a", "b"]
                ),
                "int-ordered": pd.Categorical(
                    [777777777, 888888888, 777777777, 888888888],
                    ordered=True,
                    categories=[888888888, 777777777],
                ),
                "int-unordered": pd.Categorical(
                    [777777777, 888888888, 777777777, 888888888],
                    ordered=False,
                    categories=[888888888, 777777777],
                ),
                "int-compat": pd.Categorical(
                    [777777777, 888888888, 777777777, 888888888],
                    ordered=False,
                    categories=[777777777, 888888888],
                ),
                "bool-ordered": pd.Categorical(
                    [True, False, True, False],
                    ordered=True,
                    categories=[True, False],
                ),
                "bool-unordered": pd.Categorical(
                    [True, False, True, False],
                    ordered=False,
                    categories=[True, False],
                ),
                "bool-compat": pd.Categorical(
                    [True, False, True, False],
                    ordered=False,
                    categories=[True, False],
                ),
            }
        )
        sdf.write(pa.Table.from_pandas(df))

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        assert (df == sdf.read().concat().to_pandas()).all().all()


# def test_write_categorical_dims(tmp_path):
#     """
#     Categories are not supported as dims. Here we test our handling of what we
#     do when we are given them as input.
#     """
#     schema = pa.schema(
#         [
#             ("soma_joinid", pa.int64()),
#             ("string", pa.dictionary(pa.int8(), pa.large_string())),
#         ]
#     )
#     with soma.DataFrame.create(
#         tmp_path.as_posix(),
#         schema=schema,
#         index_column_names=["soma_joinid"],
#     ) as sdf:
#         df = pd.DataFrame(
#             data={
#                 "soma_joinid": pd.Categorical([0, 1, 2, 3], categories=[0, 1, 2, 3]),
#                 "string": pd.Categorical(["a", "b", "a", "b"], categories=["b", "a"]),
#             }
#         )
#         sdf.write(pa.Table.from_pandas(df))

#     with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
#         assert (df == sdf.read().concat().to_pandas()).all().all()


@pytest.mark.parametrize("index_type", [pa.int8(), pa.int16(), pa.int32(), pa.int64()])
def test_write_categorical_dim_extend(tmp_path, index_type):
    """
    Introduce new categorical values in each subsequent write.
    """
    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("string", pa.dictionary(pa.int8(), pa.large_string())),
        ]
    )

    df1 = pd.DataFrame(
        data={
            "soma_joinid": [0, 1, 2, 3],
            "string": pd.Categorical(["a", "b", "a", "b"], categories=["b", "a"]),
        }
    )

    df2 = pd.DataFrame(
        data={
            "soma_joinid": [4, 5],
            "string": pd.Categorical(["c", "b"], categories=["b", "c"]),
        }
    )

    with soma.DataFrame.create(
        tmp_path.as_posix(),
        schema=schema,
        index_column_names=["soma_joinid"],
    ) as sdf:
        table = pa.Table.from_pandas(df1)
        dtype = pa.dictionary(index_type, pa.string())
        set_index_type = table.set_column(
            1,
            pa.field("string", dtype),
            pa.array(["a", "b", "a", "b"], dtype),
        )
        sdf.write(set_index_type)

    with soma.DataFrame.open(tmp_path.as_posix(), "w") as sdf:
        sdf.write(pa.Table.from_pandas(df2))

    # https://stackoverflow.com/questions/45639350/retaining-categorical-dtype-upon-dataframe-concatenation
    uc = union_categoricals([df1.string, df2.string])
    df1.string = pd.Categorical(df1.string, categories=uc.categories)
    df2.string = pd.Categorical(df2.string, categories=uc.categories)
    expected_df = pd.concat((df1, df2), ignore_index=True)

    with soma.DataFrame.open(tmp_path.as_posix()) as sdf:
        data = sdf.read().concat()
        assert expected_df.compare(data.to_pandas()).empty


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

        with raises_no_typeguard(ValueError):
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


@pytest.mark.parametrize("allows_duplicates", [False, True])
@pytest.mark.parametrize("consolidate", [False, True])
def test_timestamped_ops(tmp_path, allows_duplicates, consolidate):
    uri = tmp_path.as_posix()

    platform_config = {"tiledb": {"create": {"allows_duplicates": allows_duplicates}}}
    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("float", pa.float64()),
            ("string", pa.large_string()),
        ]
    )

    start = datetime.datetime(2021, 3, 10, 19, 0, tzinfo=datetime.timezone.utc)
    with soma.DataFrame.create(
        uri,
        schema=schema,
        index_column_names=["soma_joinid"],
        tiledb_timestamp=start,
        platform_config=platform_config,
    ) as sidf:
        data = {
            "soma_joinid": [0],
            "float": [100.1],
            "string": ["apple"],
        }
        sidf.write(pa.Table.from_pydict(data))
        assert sidf.tiledb_timestamp_ms == 1615402800000
        assert sidf.tiledb_timestamp.isoformat() == "2021-03-10T19:00:00+00:00"

    end = start + datetime.timedelta(minutes=3, seconds=25)
    with soma.DataFrame.open(uri=uri, mode="w", tiledb_timestamp=end) as sidf:
        data = {
            "soma_joinid": [0, 1],
            "float": [200.2, 300.3],
            "string": ["ball", "cat"],
        }
        sidf.write(pa.Table.from_pydict(data))
        assert sidf.tiledb_timestamp_ms == 1615403005000
        assert sidf.tiledb_timestamp.isoformat() == "2021-03-10T19:03:25+00:00"

    # Without consolidate:
    # * There are two fragments:
    #   o One with tiledb.fragment.FragmentInfoList[i].timestamp_range = (10, 10)
    #   o One with tiledb.fragment.FragmentInfoList[i].timestamp_range = (20, 20)
    # With consolidate:
    # * There is one fragment:
    #   o One with tiledb.fragment.FragmentInfoList[i].timestamp_range = (10, 20)
    if consolidate:
        tiledb.consolidate(uri)
        tiledb.vacuum(uri)

    # read without timestamp (i.e., after final write) & see final image
    with soma.DataFrame.open(uri) as sidf:
        table = sidf.read().concat()
        if allows_duplicates:
            assert sorted(list(x.as_py() for x in table["soma_joinid"])) == [0, 0, 1]
            assert sorted(list(x.as_py() for x in table["float"])) == [
                100.1,
                200.2,
                300.3,
            ]
            assert sorted(list(x.as_py() for x in table["string"])) == [
                "apple",
                "ball",
                "cat",
            ]
            assert sidf.count == 3
        else:
            assert list(x.as_py() for x in table["soma_joinid"]) == [0, 1]
            assert list(x.as_py() for x in table["float"]) == [200.2, 300.3]
            assert list(x.as_py() for x in table["string"]) == ["ball", "cat"]
            assert sidf.count == 2

    middle = 1615402887987
    # read at t=15 & see only the first write
    with soma.DataFrame.open(tmp_path.as_posix(), tiledb_timestamp=middle) as sidf:
        tab = sidf.read().concat()
        assert list(x.as_py() for x in tab["soma_joinid"]) == [0]
        assert list(x.as_py() for x in tab["float"]) == [100.1]
        assert list(x.as_py() for x in tab["string"]) == ["apple"]
        assert sidf.tiledb_timestamp_ms == 1615402887987
        assert sidf.tiledb_timestamp.isoformat() == "2021-03-10T19:01:27.987000+00:00"


def test_extend_enumerations(tmp_path):
    written_df = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2, 3, 4, 5], dtype=np.int64),
            "str": pd.Series(["A", "B", "A", "B", "B", "B"], dtype="category"),
            "byte": pd.Series([b"A", b"B", b"A", b"B", b"B", b"B"], dtype="category"),
            "bool": pd.Series(
                [True, False, True, False, False, False], dtype="category"
            ),
            "int64": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.int64), dtype="category"
            ),
            "uint64": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.uint64), dtype="category"
            ),
            "int32": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.int32), dtype="category"
            ),
            "uint32": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.uint32), dtype="category"
            ),
            "int16": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.int16), dtype="category"
            ),
            "uint16": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.uint16), dtype="category"
            ),
            "int8": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.int8), dtype="category"
            ),
            "uint8": pd.Series(
                np.array([0, 1, 2, 0, 1, 2], dtype=np.uint8), dtype="category"
            ),
            "float32": pd.Series(
                np.array([0, 1.1, 2.1, 0, 1.1, 2.1], dtype=np.float32), dtype="category"
            ),
            "float64": pd.Series(
                np.array([0, 1.1, 2.1, 0, 1.1, 2.1], dtype=np.float64), dtype="category"
            ),
            "float64_w_non_finite": pd.Series(
                np.array([0, 1.1, 2.1, 0, np.Inf, np.NINF], dtype=np.float64),
                dtype="category",
            ),
            "str_ordered": pd.Series(
                pd.Categorical(
                    ["A", "B", "A", "B", "B", "B"],
                    categories=["B", "A", "C"],
                    ordered=True,
                ),
            ),
            "int64_ordered": pd.Series(
                pd.Categorical(
                    [1, 2, 3, 3, 2, 1],
                    categories=np.array([3, 2, 1], dtype=np.int64),
                    ordered=True,
                ),
            ),
            "uint64_ordered": pd.Series(
                pd.Categorical(
                    [1, 2, 3, 3, 2, 1],
                    categories=np.array([3, 2, 1], dtype=np.uint64),
                    ordered=True,
                ),
            ),
            "float64_ordered": pd.Series(
                pd.Categorical(
                    [0, 1.1, 2.1, 0, 1.1, 2.1],
                    categories=np.array([1.1, 0, 2.1], dtype=np.float64),
                    ordered=True,
                ),
            ),
        },
    )

    schema = pa.Schema.from_pandas(written_df, preserve_index=False)

    with soma.DataFrame.create(str(tmp_path), schema=schema) as soma_dataframe:
        tbl = pa.Table.from_pandas(written_df, preserve_index=False)
        soma_dataframe.write(tbl)

    with soma.open(str(tmp_path)) as soma_dataframe:
        readback_df = soma_dataframe.read().concat().to_pandas()
        for c in readback_df:
            assert readback_df[c].dtype == written_df[c].dtype
            if readback_df[c].dtype == "category":
                assert (
                    readback_df[c].cat.categories.dtype
                    == written_df[c].cat.categories.dtype
                )
            assert (readback_df[c] == written_df[c]).all()


def test_multiple_writes_with_str_enums(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            (
                "obs",
                pa.dictionary(
                    index_type=pa.int8(), value_type=pa.string(), ordered=False
                ),
            ),
        ]
    )
    soma.DataFrame.create(uri, schema=schema).close()

    df1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2], dtype=np.int64),
            "obs": pd.Series(["A", "B", "A"], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df1, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    df2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([3, 4, 5], dtype=np.int64),
            "obs": pd.Series(["B", "C", "B"], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df2, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    # https://stackoverflow.com/questions/45639350/retaining-categorical-dtype-upon-dataframe-concatenation
    uc = union_categoricals([df1.obs, df2.obs])
    df1.obs = pd.Categorical(df1.obs, categories=uc.categories)
    df2.obs = pd.Categorical(df2.obs, categories=uc.categories)
    expected_df = pd.concat((df1, df2), ignore_index=True)

    assert df.equals(expected_df)

    # No new enumerations introduced but we still need to remap
    # the indexes
    df3 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([6, 7], dtype=np.int64),
            "obs": pd.Series(["C", "C"], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df3, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    uc = union_categoricals([df1.obs, df2.obs, df3.obs])
    df1.obs = pd.Categorical(df1.obs, categories=uc.categories)
    df2.obs = pd.Categorical(df2.obs, categories=uc.categories)
    df3.obs = pd.Categorical(df3.obs, categories=uc.categories)
    expected_df = pd.concat((df1, df2, df3), ignore_index=True)

    assert df.equals(expected_df)


def test_multiple_writes_with_int_enums(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            (
                "obs",
                pa.dictionary(
                    index_type=pa.int8(), value_type=pa.int64(), ordered=False
                ),
            ),
        ]
    )
    soma.DataFrame.create(uri, schema=schema).close()

    df1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2], dtype=np.int64),
            "obs": pd.Series([1, 2, 1], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df1, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    df2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([3, 4, 5], dtype=np.int64),
            "obs": pd.Series([2, 3, 2], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df2, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    # https://stackoverflow.com/questions/45639350/retaining-categorical-dtype-upon-dataframe-concatenation
    uc = union_categoricals([df1.obs, df2.obs])
    df1.obs = pd.Categorical(df1.obs, categories=uc.categories)
    df2.obs = pd.Categorical(df2.obs, categories=uc.categories)
    expected_df = pd.concat((df1, df2), ignore_index=True)

    assert df.equals(expected_df)

    # No new enumerations introduced but we still need to remap
    # the indexes
    df3 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([6, 7], dtype=np.int64),
            "obs": pd.Series([3, 3], dtype="category"),
        }
    )
    tbl = pa.Table.from_pandas(df3, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    uc = union_categoricals([df1.obs, df2.obs, df3.obs])
    df1.obs = pd.Categorical(df1.obs, categories=uc.categories)
    df2.obs = pd.Categorical(df2.obs, categories=uc.categories)
    df3.obs = pd.Categorical(df3.obs, categories=uc.categories)
    expected_df = pd.concat((df1, df2, df3), ignore_index=True)

    assert df.equals(expected_df)


def test_multichunk(tmp_path):
    uri = tmp_path.as_posix()

    # --- three dataframes, all with identical schema
    df_0 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2, 3], dtype=np.int64),
            "obs": pd.Series(["A", "B", "A", "B"], dtype="str"),
        }
    )
    df_1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([4, 5, 6, 7], dtype=np.int64),
            "obs": pd.Series(["A", "A", "B", "B"], dtype="str"),
        }
    )
    df_2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([8, 9, 10, 11], dtype=np.int64),
            "obs": pd.Series(["B", "C", "B", "C"], dtype="str"),
        }
    )
    expected_df = pd.concat([df_0, df_1, df_2], ignore_index=True)

    soma.DataFrame.create(
        uri, schema=pa.Schema.from_pandas(df_0, preserve_index=False)
    ).close()

    with soma.open(uri, mode="w") as A:
        # one-chunk table
        A.write(pa.Table.from_pandas(df_0, preserve_index=False))

    with soma.open(uri, mode="w") as A:
        # two-chunk table
        A.write(
            pa.concat_tables(
                [
                    pa.Table.from_pandas(df_1, preserve_index=False),
                    pa.Table.from_pandas(df_2, preserve_index=False),
                ]
            )
        )

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    assert df.equals(expected_df)


def test_multichunk_with_enums(tmp_path):
    uri = tmp_path.as_posix()

    # --- three dataframes, all with identical schema
    df_0 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2, 3], dtype=np.int64),
            "obs": pd.Series(["A", "B", "A", "B"], dtype="category"),
        }
    )
    df_1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([4, 5, 6, 7], dtype=np.int64),
            "obs": pd.Series(["A", "A", "B", "B"], dtype="category"),
        }
    )
    df_2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series([8, 9, 10, 11], dtype=np.int64),
            "obs": pd.Series(["B", "C", "B", "C"], dtype="category"),
        }
    )
    expected_df = pd.concat([df_0, df_1, df_2], ignore_index=True)

    soma.DataFrame.create(
        uri, schema=pa.Schema.from_pandas(df_0, preserve_index=False)
    ).close()

    with soma.open(uri, mode="w") as A:
        # one-chunk table
        A.write(pa.Table.from_pandas(df_0, preserve_index=False))

    with soma.open(uri, mode="w") as A:
        # two-chunk table
        A.write(
            pa.concat_tables(
                [
                    pa.Table.from_pandas(df_1, preserve_index=False),
                    pa.Table.from_pandas(df_2, preserve_index=False),
                ]
            )
        )

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()

    # https://stackoverflow.com/questions/45639350/retaining-categorical-dtype-upon-dataframe-concatenation
    uc = union_categoricals([df_0.obs, df_1.obs, df_2.obs])
    df_0.obs = pd.Categorical(df_0.obs, categories=uc.categories)
    df_1.obs = pd.Categorical(df_1.obs, categories=uc.categories)
    df_2.obs = pd.Categorical(df_2.obs, categories=uc.categories)
    expected_df = pd.concat((df_0, df_1, df_2), ignore_index=True)

    assert df.equals(expected_df)


def test_enum_extend_past_numerical_limit(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            (
                "obs",
                pa.dictionary(
                    index_type=pa.int8(), value_type=pa.large_string(), ordered=False
                ),
            ),
        ]
    )
    soma.DataFrame.create(uri, schema=schema).close()

    n_elem = 132
    n_cats = 127
    df1 = pd.DataFrame(
        {
            "soma_joinid": pd.Series(np.arange(n_elem), dtype=np.int64),
            "obs": pd.Series(
                [f"enum_{i % n_cats}" for i in range(n_elem)], dtype="category"
            ),
        }
    )

    # use max number of possible categories
    tbl = pa.Table.from_pandas(df1, preserve_index=False)
    with soma.open(uri, mode="w") as A:
        A.write(tbl)

    more_elem = 4
    df2 = pd.DataFrame(
        {
            "soma_joinid": pd.Series(
                np.arange(n_elem, n_elem + more_elem), dtype=np.int64
            ),
            "obs": pd.Series(["TEST"] * more_elem, dtype="category"),
        }
    )

    # cannot add additional categories as already maxed out earlier
    tbl = pa.Table.from_pandas(df2, preserve_index=False)
    with pytest.raises(soma.SOMAError):
        with soma.open(uri, mode="w") as A:
            A.write(tbl)


def test_write_str_empty_ned(tmp_path):
    tmp_path.as_posix()


def test_enum_schema_report(tmp_path):
    uri = tmp_path.as_posix()

    pandas_df = pd.DataFrame(
        {
            "soma_joinid": pd.Series([0, 1, 2, 3, 4, 5], dtype=np.int64),
            "int_cat": pd.Series([10, 20, 10, 20, 20, 20], dtype="category"),
            "int": pd.Series([10, 20, 10, 20, 20, 20]),
            "str_cat": pd.Series(["A", "B", "A", "B", "B", "B"], dtype="category"),
            "str": pd.Series(["A", "B", "A", "B", "B", "B"]),
            "byte_cat": pd.Series(
                [b"A", b"B", b"A", b"B", b"B", b"B"], dtype="category"
            ),
            "byte": pd.Series([b"A", b"B", b"A", b"B", b"B", b"B"]),
        },
    )

    arrow_schema = pa.Schema.from_pandas(pandas_df, preserve_index=False)

    with soma.DataFrame.create(uri, schema=arrow_schema) as sdf:
        arrow_table = pa.Table.from_pandas(pandas_df, preserve_index=False)
        sdf.write(arrow_table)

    # Double-check against TileDB-Py reporting
    with tiledb.open(uri) as A:
        for i in range(A.schema.nattr):
            attr = A.schema.attr(i)
            try:
                index_type = attr.dtype
                value_type = A.enum(attr.name).dtype
            except tiledb.cc.TileDBError:
                pass  # not an enum attr
            if attr.name == "int_cat":
                assert index_type.name == "int8"
                assert value_type.name == "int64"
            elif attr.name == "str_cat":
                assert index_type.name == "int8"
                assert value_type.name == "str32"
            elif attr.name == "byte_cat":
                assert index_type.name == "int8"
                assert value_type.name == "bytes8"

    # Verify SOMA Arrow schema
    with soma.open(uri) as sdf:
        f = sdf.schema.field("int_cat")
        assert f.type.index_type == pa.int8()
        assert f.type.value_type == pa.int64()

        f = sdf.schema.field("str_cat")
        assert f.type.index_type == pa.int8()
        assert f.type.value_type == pa.string()

        f = sdf.schema.field("byte_cat")
        assert f.type.index_type == pa.int8()
        assert f.type.value_type == pa.binary()


def test_nullable(tmp_path):
    uri = tmp_path.as_posix()

    # Arrow fields are nullable by default.  They can be explicitly set nullable
    # or non-nullable via the nullable kwarg to pa.field.  Also, they can be
    # explicitly set nullable via metadata. The latter, if present, overrides
    # the former.
    asch = pa.schema(
        [
            pa.field("int", pa.int32()),
            pa.field("bool", pa.bool_()),
            pa.field("ord", pa.dictionary(pa.int64(), pa.string())),
            pa.field("no-meta-flag-unspecified", pa.int32()),
            pa.field("no-meta-flag-true", pa.int32(), nullable=True),
            pa.field("no-meta-flag-false", pa.int32(), nullable=False),
            pa.field("yes-meta-flag-unspecified", pa.int32()),
            pa.field("yes-meta-flag-true", pa.int32(), nullable=True),
            pa.field("yes-meta-flag-false", pa.int32(), nullable=False),
        ],
        metadata={
            "int": "nullable",
            "bool": "nullable",
            "ord": "nullable",
            "yes-meta-flag-unspecified": "nullable",
            "yes-meta-flag-true": "nullable",
            "yes-meta-flag-false": "nullable",
        },
    )

    pydict = {}
    pydict["soma_joinid"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    pydict["int"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["bool"] = [True, True, True, False, True, False, None, False, None, None]
    pydict["ord"] = pd.Categorical(
        ["g1", "g2", "g3", None, "g2", "g3", "g1", None, "g3", "g1"]
    )
    pydict["no-meta-flag-unspecified"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["no-meta-flag-true"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["no-meta-flag-false"] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    pydict["yes-meta-flag-unspecified"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["yes-meta-flag-true"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["yes-meta-flag-false"] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    data = pa.Table.from_pydict(pydict)

    with soma.DataFrame.create(uri, schema=asch) as sdf:
        sdf.write(data)

    with soma.DataFrame.open(uri, "r") as sdf:
        df = sdf.read().concat().to_pandas()
        assert df.compare(data.to_pandas()).empty


def test_only_evolve_schema_when_enmr_is_extended(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            pa.field("foo", pa.dictionary(pa.int64(), pa.large_string())),
            pa.field("bar", pa.large_string()),
        ]
    )

    # +1 creating the schema
    # +1 evolving the schema
    with soma.DataFrame.create(uri, schema=schema) as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["foo"] = pd.Categorical(["a", "bb", "ccc", "bb", "a"])
        data["bar"] = ["cat", "dog", "cat", "cat", "cat"]
        sdf.write(pa.Table.from_pydict(data))

    # +1 evolving the schema
    with soma.DataFrame.open(uri, "w") as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["foo"] = pd.Categorical(["a", "bb", "ccc", "d", "a"])
        data["bar"] = ["cat", "dog", "cat", "cat", "cat"]
        sdf.write(pa.Table.from_pydict(data))

    # +0 no changes to enumeration values
    with soma.DataFrame.open(uri, "w") as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["foo"] = pd.Categorical(["a", "bb", "ccc", "d", "a"])
        data["bar"] = ["cat", "dog", "cat", "cat", "cat"]
        sdf.write(pa.Table.from_pydict(data))

    # +0 no changes enumeration values
    with soma.DataFrame.open(uri, "w") as sdf:
        data = {}
        data["soma_joinid"] = [0, 1, 2, 3, 4]
        data["foo"] = pd.Categorical(["a", "bb", "ccc", "d", "d"])
        data["bar"] = ["cat", "dog", "cat", "cat", "cat"]
        sdf.write(pa.Table.from_pydict(data))

    # total 3 fragment files

    vfs = tiledb.VFS()
    # subtract 1 for the __schema/__enumerations directory;
    # only looking at fragment files
    assert len(vfs.ls(os.path.join(uri, "__schema"))) - 1 == 3


def test_fix_update_dataframe_with_var_strings(tmp_path):
    uri = tmp_path.as_posix()

    tbl = pa.table(
        {
            "soma_joinid": pa.array([0, 1, 2, 3], pa.int64()),
            "mystring": pa.array(["a", "bb", "ccc", "dddd"], pa.large_utf8()),
            "myint": pa.array([33, 44, 55, 66], pa.int32()),
            "myfloat": pa.array([4.5, 5.5, 6.5, 7.5], pa.float32()),
        }
    )

    with soma.DataFrame.create(uri, schema=tbl.schema) as sdf:
        sdf.write(tbl)

    with soma.DataFrame.open(uri, "r") as sdf:
        updated_sdf = sdf.read().concat().to_pandas()
    updated_sdf["newattr"] = np.array(["a", "b", "c", "d"])

    with soma.DataFrame.open(uri, "w") as sdf:
        soma.io.ingest._update_dataframe(
            sdf,
            updated_sdf,
            "testing",
            platform_config=None,
            context=None,
            default_index_name="mystring",
        )

    with soma.DataFrame.open(uri, "r") as sdf:
        results = sdf.read().concat().to_pandas()
        assert results.equals(updated_sdf)
