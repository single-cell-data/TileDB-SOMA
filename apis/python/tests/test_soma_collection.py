import os

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma.soma_exception import SOMADoesNotExistError


# ----------------------------------------------------------------
def create_and_populate_dataframe(dataframe: soma.SOMADataFrame) -> None:

    arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )

    dataframe.create(schema=arrow_schema)

    pydict = {}
    pydict["soma_rowid"] = [0, 1, 2, 3, 4]
    pydict["soma_joinid"] = [0, 1, 2, 3, 4]
    pydict["foo"] = [10, 20, 30, 40, 50]
    pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
    pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
    rb = pa.Table.from_pydict(pydict)
    dataframe.write(rb)


# ----------------------------------------------------------------
def create_and_populate_sparse_nd_array(
    sparse_nd_array: soma.SOMASparseNdArray,
) -> None:
    nr = 10
    nc = 20
    sparse_nd_array.create(pa.int64(), [nr, nc])

    tensor = pa.SparseCOOTensor.from_numpy(
        data=np.asarray([7, 8, 9]),
        coords=[[0, 1], [2, 3], [3, 4]],
        shape=(nr, nc),
    )
    sparse_nd_array.write_sparse_tensor(tensor)


# ----------------------------------------------------------------
def test_soma_collection_basic(tmp_path):
    basedir = tmp_path.as_uri()
    collection = soma.SOMACollection(basedir)

    # ----------------------------------------------------------------
    assert not collection.exists()
    with pytest.raises(soma.SOMADoesNotExistError):
        assert "foobar" not in collection

    collection.create()
    assert collection.exists()
    assert collection.uri == basedir
    assert "foobar" not in collection

    dataframe = soma.SOMADataFrame(os.path.join(basedir, "sdf"), parent=collection)
    create_and_populate_dataframe(dataframe)

    sparse_nd_array = soma.SOMASparseNdArray(
        os.path.join(basedir, "snda"), parent=collection
    )
    create_and_populate_sparse_nd_array(sparse_nd_array)

    collection.set("sdf", dataframe)
    collection.set("snda", sparse_nd_array)

    # ----------------------------------------------------------------
    readback_collection = soma.SOMACollection(collection.uri)
    assert len(readback_collection) == 2

    readback_dataframe = readback_collection.get("sdf")
    with readback_dataframe._tiledb_open() as A:
        assert len(A.df[:]) == 5

    readback_sparse_nd_array = readback_collection.get("snda")
    with readback_sparse_nd_array._tiledb_open() as A:
        assert len(A.df[:]) == 3


@pytest.fixture(
    scope="function",
    params=[
        "SOMACollection",
        "SOMADataFrame",
        "SOMAIndexedDataFrame",
        "SOMADenseNdArray",
        "SOMASparseNdArray",
    ],
)
def soma_object(request, tmp_path):
    """
    Make an empty test object of the given foundational class name.
    """
    uri = tmp_path.joinpath("object").as_uri()
    class_name = request.param

    if class_name == "SOMACollection":
        so = soma.SOMACollection(uri=uri)
        so.create()

    elif class_name == "SOMADataFrame":
        so = soma.SOMADataFrame(uri=uri)
        so.create(pa.schema([("A", pa.int32()), ("B", pa.large_string())]))

    elif class_name == "SOMAIndexedDataFrame":
        so = soma.SOMAIndexedDataFrame(uri=uri)
        so.create(
            schema=pa.schema([("C", pa.float32()), ("D", pa.uint32())]),
            index_column_names=["D"],
        )

    elif class_name == "SOMADenseNdArray":
        so = soma.SOMADenseNdArray(uri=uri)
        so.create(type=pa.float64(), shape=(100, 10, 1))

    elif class_name == "SOMASparseNdArray":
        so = soma.SOMASparseNdArray(uri=uri)
        so.create(type=pa.int8(), shape=(11,))

    assert so is not None, f"Unknown class name: {class_name}"
    yield so
    so.delete()


def test_soma_collection_mapping(soma_object, tmp_path):
    c = soma.SOMACollection(uri=(tmp_path / "collection").as_uri())
    assert not c.exists()
    with pytest.raises(soma.SOMADoesNotExistError):
        assert "foobar" not in c

    c.create()
    assert c.exists()
    assert "foobar" not in c

    assert soma_object.exists()
    c["mumble"] = soma_object
    assert "mumble" in c
    assert c["mumble"] == soma_object

    assert list(c.keys()) == ["mumble"]
    assert list(c.items()) == [("mumble", soma_object)]
    assert list(c.values()) == [soma_object]
    assert list(k for k in c) == list(c.keys())
    assert len(c) == 1
    assert [k for k in c] == ["mumble"]

    del c["mumble"]
    assert "mumble" not in c
    assert not c.get("mumble", False)

    c.delete()
    assert not c.exists()
    assert not (tmp_path / "collection").exists()


@pytest.mark.parametrize("relative", [False, True])
def test_collection_repr(tmp_path, relative):
    a = soma.SOMACollection(uri=(tmp_path / "A").as_uri())
    a.create()
    assert a.exists()
    assert a.uri == (tmp_path / "A").as_uri()

    b = soma.SOMACollection(uri=(tmp_path / "A" / "B").as_uri())
    b.create()
    assert b.exists()
    assert b.uri == (tmp_path / "A" / "B").as_uri()

    a.set("Another_Name", b, relative=relative)
    assert list(a.keys()) == ["Another_Name"]
    assert (
        a.__repr__()
        == f'SOMACollection(uri="{a.uri}"):\n  "Another_Name": SOMACollection(uri="{b.uri}")'
    )
    assert a["Another_Name"].uri == (tmp_path / "A" / "B").as_uri()
    del a

    # re-open, reconfirm
    aPrime = soma.SOMACollection(uri=(tmp_path / "A").as_uri())
    assert list(aPrime.keys()) == ["Another_Name"]
    assert (
        aPrime.__repr__()
        == f'SOMACollection(uri="{aPrime.uri}"):\n  "Another_Name": SOMACollection(uri="{b.uri}")'
    )
    assert aPrime["Another_Name"].uri == (tmp_path / "A" / "B").as_uri()
    del aPrime

    # move container, re-confirm wrt "relative" value
    os.rename((tmp_path / "A"), (tmp_path / "A_moved"))
    aMoved = soma.SOMACollection(uri=(tmp_path / "A_moved").as_uri())
    assert list(aMoved.keys()) == ["Another_Name"]
    if relative:
        assert aMoved["Another_Name"].uri == (tmp_path / "A_moved" / "B").as_uri()
        assert (
            aMoved.__repr__()
            == f'SOMACollection(uri="{aMoved.uri}"):\n  "Another_Name": SOMACollection(uri="{aMoved["Another_Name"].uri}")'
        )
    else:
        with pytest.raises(KeyError):
            aMoved["Another_Name"].uri

    del aMoved


def test_soma_collection_update_on_set(tmp_path):
    """
    SOMACollection.__setattr__ (and .set) have update semantics. Underlying
    tiledb.Group only has add/del. Verify.
    """

    sc = soma.SOMACollection(tmp_path.as_uri()).create()
    A = soma.SOMADenseNdArray(uri=(tmp_path / "A").as_uri()).create(
        type=pa.float64(), shape=(100, 10, 1)
    )
    B = soma.SOMADenseNdArray(uri=(tmp_path / "B").as_uri()).create(
        type=pa.float64(), shape=(100, 10, 1)
    )
    assert sc.exists()
    assert set(sc.keys()) == set([])

    sc["A"] = A
    assert set(sc.keys()) == set(["A"])
    assert sc["A"] == A

    sc["A"] = B
    assert set(sc.keys()) == set(["A"])
    assert sc["A"] == B


def test_exceptions_on_not_created(tmp_path):
    """
    When the collection has not been created, we should get
    a meaningful error
    """
    sc = soma.SOMACollection(tmp_path.as_uri())

    A = soma.SOMADenseNdArray(uri=(tmp_path / "A").as_uri()).create(
        type=pa.float64(), shape=(100, 10, 1)
    )
    with pytest.raises(SOMADoesNotExistError):
        sc["A"] = A
    with pytest.raises(SOMADoesNotExistError):
        del sc["A"]
    with pytest.raises(SOMADoesNotExistError):
        assert "A" not in sc
    with pytest.raises(SOMADoesNotExistError):
        list(sc)
    with pytest.raises(SOMADoesNotExistError):
        len(sc)
