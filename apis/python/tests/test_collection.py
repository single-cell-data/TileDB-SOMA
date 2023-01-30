import os

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma import factory
from tiledbsoma.exception import DoesNotExistError


# ----------------------------------------------------------------
def create_and_populate_dataframe(path: str) -> soma.DataFrame:

    arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )

    with soma.DataFrame.create(path, schema=arrow_schema) as df:
        pydict = {}
        pydict["soma_joinid"] = [0, 1, 2, 3, 4]
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.Table.from_pydict(pydict)
        df.write(rb)

    return factory.open(path)


# ----------------------------------------------------------------
def create_and_populate_sparse_nd_array(path: str) -> soma.SparseNDArray:
    nr = 10
    nc = 20
    with soma.SparseNDArray.create(
        path, type=pa.int64(), shape=(nr, nc)
    ) as sparse_nd_array:

        tensor = pa.SparseCOOTensor.from_numpy(
            data=np.asarray([7, 8, 9]),
            coords=[[0, 1], [2, 3], [3, 4]],
            shape=(nr, nc),
        )
        sparse_nd_array.write(tensor)
    return factory.open(path)


# ----------------------------------------------------------------
def test_collection_basic(tmp_path):
    basedir = tmp_path.as_uri()
    # ----------------------------------------------------------------
    with pytest.raises(soma.DoesNotExistError):
        soma.Collection.open(basedir)

    with soma.Collection.create(basedir) as collection:
        assert collection.uri == basedir
        assert "foobar" not in collection

        with create_and_populate_dataframe(os.path.join(basedir, "sdf")) as sdf:
            collection["sdf"] = sdf
        with create_and_populate_sparse_nd_array(os.path.join(basedir, "snda")) as snda:
            collection["snda"] = snda

    # ----------------------------------------------------------------
    readback_collection = soma.Collection.open(collection.uri)
    assert len(readback_collection) == 2

    with readback_collection["sdf"] as sdf:
        assert len(sdf._handle.reader.df[:]) == 5

    with readback_collection["snda"] as snda:
        assert len(snda._handle.reader.df[:]) == 3


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
        so = soma.Collection.create(uri)

    elif class_name == "DataFrame":
        so = soma.DataFrame.create(
            uri,
            schema=pa.schema([("C", pa.float32()), ("D", pa.uint32())]),
            index_column_names=["D"],
        )

    elif class_name == "DenseNDArray":
        so = soma.DenseNDArray.create(uri, type=pa.float64(), shape=(100, 10, 1))

    elif class_name == "SparseNDArray":
        so = soma.SparseNDArray.create(uri, type=pa.int8(), shape=(11,))
    else:
        raise ValueError(f"don't know how to make {class_name}")

    yield so
    so.close()


def test_collection_mapping(soma_object, tmp_path):
    collection_path = (tmp_path / "collection").as_uri()
    with pytest.raises(soma.DoesNotExistError):
        soma.Collection.open(collection_path)

    c = soma.Collection.create(collection_path)
    assert "foobar" not in c

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

    # XXX create delete global function
    # assert not (tmp_path / "collection").exists()


@pytest.mark.parametrize("relative", [False, True])
def test_collection_repr(tmp_path, relative):
    a = soma.Collection.create((tmp_path / "A").as_uri())
    assert a.uri == (tmp_path / "A").as_uri()

    b = soma.Collection.create((tmp_path / "A" / "B").as_uri())
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
    aPrime = soma.Collection.open((tmp_path / "A").as_uri())
    assert list(aPrime.keys()) == ["Another_Name"]
    assert (
        aPrime.__repr__()
        == f'SOMACollection(uri="{aPrime.uri}"):\n  "Another_Name": SOMACollection(uri="{b.uri}")'
    )
    assert aPrime["Another_Name"].uri == (tmp_path / "A" / "B").as_uri()
    del aPrime

    # move container, re-confirm wrt "relative" value
    os.rename((tmp_path / "A"), (tmp_path / "A_moved"))
    aMoved = soma.Collection.open((tmp_path / "A_moved").as_uri())
    assert list(aMoved.keys()) == ["Another_Name"]
    if relative:
        assert aMoved["Another_Name"].uri == (tmp_path / "A_moved" / "B").as_uri()
        assert (
            aMoved.__repr__()
            == f'SOMACollection(uri="{aMoved.uri}"):\n  "Another_Name": SOMACollection(uri="{aMoved["Another_Name"].uri}")'
        )
    else:
        with pytest.raises(DoesNotExistError):
            aMoved["Another_Name"].uri

    del aMoved


def test_collection_update_on_set(tmp_path):
    """
    Collection.__setattr__ (and .set) have update semantics. Underlying
    tiledb.Group only has add/del. Verify.
    """

    sc = soma.Collection.create(tmp_path.as_uri())
    A = soma.DenseNDArray.create(
        (tmp_path / "A").as_uri(), type=pa.float64(), shape=(100, 10, 1)
    )
    B = soma.DenseNDArray.create(
        uri=(tmp_path / "B").as_uri(), type=pa.float64(), shape=(100, 10, 1)
    )
    assert set(sc.keys()) == set([])

    sc["A"] = A
    assert set(sc.keys()) == set(["A"])
    assert sc["A"] == A

    sc["A"] = B
    assert set(sc.keys()) == set(["A"])
    assert sc["A"] == B
