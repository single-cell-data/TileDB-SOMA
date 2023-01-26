import os

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma.exception import DoesNotExistError
from tiledbsoma.options import SOMATileDBContext


# ----------------------------------------------------------------
def create_and_populate_dataframe(dataframe: soma.DataFrame) -> None:

    arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )

    dataframe.create_legacy(schema=arrow_schema)

    pydict = {}
    pydict["soma_joinid"] = [0, 1, 2, 3, 4]
    pydict["foo"] = [10, 20, 30, 40, 50]
    pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
    pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
    rb = pa.Table.from_pydict(pydict)
    dataframe.write(rb)


# ----------------------------------------------------------------
def create_and_populate_sparse_nd_array(
    sparse_nd_array: soma.SparseNDArray,
) -> None:
    nr = 10
    nc = 20
    sparse_nd_array.create_legacy(pa.int64(), [nr, nc])

    tensor = pa.SparseCOOTensor.from_numpy(
        data=np.asarray([7, 8, 9]),
        coords=[[0, 1], [2, 3], [3, 4]],
        shape=(nr, nc),
    )
    sparse_nd_array.write(tensor)


# ----------------------------------------------------------------
def test_collection_basic(tmp_path):
    basedir = tmp_path.as_uri()
    collection = soma.Collection(basedir)

    # ----------------------------------------------------------------
    assert not collection.exists()
    with pytest.raises(soma.DoesNotExistError):
        assert "foobar" not in collection

    collection.create_legacy()
    assert collection.exists()
    assert collection.uri == basedir
    assert "foobar" not in collection

    dataframe = soma.DataFrame(os.path.join(basedir, "sdf"))
    create_and_populate_dataframe(dataframe)

    sparse_nd_array = soma.SparseNDArray(os.path.join(basedir, "snda"))
    create_and_populate_sparse_nd_array(sparse_nd_array)

    collection.set("sdf", dataframe)
    collection.set("snda", sparse_nd_array)

    # ----------------------------------------------------------------
    readback_collection = soma.Collection(collection.uri)
    assert len(readback_collection) == 2

    with readback_collection.get("sdf").open_legacy() as sdf:
        assert len(sdf._tiledb_obj.df[:]) == 5

    with readback_collection.get("snda").open_legacy() as snda:
        assert len(snda._tiledb_obj.df[:]) == 3


@pytest.fixture(
    scope="function",
    params=[
        "Collection",
        "DataFrame",
        "DenseNDArray",
        "SparseNDArray",
    ],
)
def soma_object(request, tmp_path):
    """
    Make an empty test object of the given foundational class name.
    """
    uri = tmp_path.joinpath("object").as_uri()
    class_name = request.param

    if class_name == "Collection":
        so = soma.Collection(uri=uri)
        so.create_legacy()

    elif class_name == "DataFrame":
        so = soma.DataFrame(uri=uri)
        so.create_legacy(
            schema=pa.schema([("C", pa.float32()), ("D", pa.uint32())]),
            index_column_names=["D"],
        )

    elif class_name == "DenseNDArray":
        so = soma.DenseNDArray(uri=uri)
        so.create_legacy(type=pa.float64(), shape=(100, 10, 1))

    elif class_name == "SparseNDArray":
        so = soma.SparseNDArray(uri=uri)
        so.create_legacy(type=pa.int8(), shape=(11,))

    assert so is not None, f"Unknown class name: {class_name}"
    yield so
    so.delete()


def test_collection_mapping(soma_object, tmp_path):
    c = soma.Collection(uri=(tmp_path / "collection").as_uri())
    assert not c.exists()
    with pytest.raises(soma.DoesNotExistError):
        assert "foobar" not in c

    c.create_legacy()
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
    a = soma.Collection(uri=(tmp_path / "A").as_uri())
    a.create_legacy()
    assert a.exists()
    assert a.uri == (tmp_path / "A").as_uri()

    b = soma.Collection(uri=(tmp_path / "A" / "B").as_uri())
    b.create_legacy()
    assert b.exists()
    assert b.uri == (tmp_path / "A" / "B").as_uri()

    a.set("Another_Name", b, use_relative_uri=relative)
    assert list(a.keys()) == ["Another_Name"]
    assert (
        a.__repr__()
        == f'SOMACollection(uri="{a.uri}"):\n  "Another_Name": SOMACollection(uri="{b.uri}")'
    )
    assert a["Another_Name"].uri == (tmp_path / "A" / "B").as_uri()
    del a

    # re-open, reconfirm
    aPrime = soma.Collection(uri=(tmp_path / "A").as_uri())
    assert list(aPrime.keys()) == ["Another_Name"]
    assert (
        aPrime.__repr__()
        == f'SOMACollection(uri="{aPrime.uri}"):\n  "Another_Name": SOMACollection(uri="{b.uri}")'
    )
    assert aPrime["Another_Name"].uri == (tmp_path / "A" / "B").as_uri()
    del aPrime

    # move container, re-confirm wrt "relative" value
    os.rename((tmp_path / "A"), (tmp_path / "A_moved"))
    aMoved = soma.Collection(uri=(tmp_path / "A_moved").as_uri())
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


def test_collection_update_on_set(tmp_path):
    """
    Collection.__setattr__ (and .set) have update semantics. Underlying
    tiledb.Group only has add/del. Verify.
    """

    sc = soma.Collection(tmp_path.as_uri()).create_legacy()
    A = soma.DenseNDArray(uri=(tmp_path / "A").as_uri()).create_legacy(
        type=pa.float64(), shape=(100, 10, 1)
    )
    B = soma.DenseNDArray(uri=(tmp_path / "B").as_uri()).create_legacy(
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
    sc = soma.Collection(tmp_path.as_uri())

    A = soma.DenseNDArray(uri=(tmp_path / "A").as_uri()).create_legacy(
        type=pa.float64(), shape=(100, 10, 1)
    )
    with pytest.raises(DoesNotExistError):
        sc["A"] = A
    with pytest.raises(DoesNotExistError):
        del sc["A"]
    with pytest.raises(DoesNotExistError):
        assert "A" not in sc
    with pytest.raises(DoesNotExistError):
        list(sc)
    with pytest.raises(DoesNotExistError):
        len(sc)


def test_timestamp_inheritance(tmp_path):
    """
    When we specify read/write timestamps in SomaTileDBContext supplied to collection, those are
    inherited by elements accessed via the collection.
    """

    # create collection & A @ t=10
    ctx_write10 = SOMATileDBContext(_write_timestamp=10)
    sc = soma.Collection(tmp_path.as_uri(), context=ctx_write10).create_legacy()
    A = soma.DenseNDArray(uri=(tmp_path / "A").as_posix(), context=ctx_write10)
    A.create_legacy(type=pa.uint8(), shape=(2, 2))
    A.write(
        (slice(0, 2), slice(0, 2)),
        pa.Tensor.from_numpy(np.zeros((2, 2), dtype=np.uint8)),
    )
    sc["A"] = A

    # access A via collection @ t=20 and write more into it
    sc = soma.Collection(
        tmp_path.as_uri(), context=SOMATileDBContext(_write_timestamp=20)
    )
    sc["A"].write(
        (slice(0, 1), slice(0, 1)),
        pa.Tensor.from_numpy(np.ones((1, 1), dtype=np.uint8)),
    )

    # open A via collection with no or late timestamp => A should reflect both writes
    sc = soma.Collection(tmp_path.as_uri())
    assert sc["A"].read((slice(None), slice(None))).to_numpy().tolist() == [
        [1, 0],
        [0, 0],
    ]
    sc = soma.Collection(
        tmp_path.as_uri(), context=SOMATileDBContext(read_timestamp=25)
    )
    assert sc["A"].read((slice(None), slice(None))).to_numpy().tolist() == [
        [1, 0],
        [0, 0],
    ]

    # open A via collection with read_timestamp=15 => A should reflect first write only
    sc = soma.Collection(
        tmp_path.as_uri(), context=SOMATileDBContext(read_timestamp=15)
    )
    assert sc["A"].read((slice(None), slice(None))).to_numpy().tolist() == [
        [0, 0],
        [0, 0],
    ]
