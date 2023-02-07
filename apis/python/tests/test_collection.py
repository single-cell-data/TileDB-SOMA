import os
import pathlib
from typing import List, TypeVar, Union

import numpy as np
import pyarrow as pa
import pytest
from typing_extensions import Literal

import tiledbsoma as soma
from tiledbsoma import collection, factory, tiledb_object
from tiledbsoma.exception import DoesNotExistError
from tiledbsoma.options import SOMATileDBContext


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

    assert soma.Collection.exists(basedir)
    assert not soma.Experiment.exists(basedir)
    assert not soma.DenseNDArray.exists(basedir)

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


@pytest.mark.parametrize("relative", [False, True])
def test_collection_repr(tmp_path, relative) -> None:
    a = soma.Collection.create((tmp_path / "A").as_uri())
    assert a.uri == (tmp_path / "A").as_uri()

    b_uri = "B" if relative else (tmp_path / "A" / "B").as_uri()
    b = a.add_new_collection("Another_Name", uri=b_uri)

    assert b.uri == (tmp_path / "A" / "B").as_uri()

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


def test_cascading_close(tmp_path: pathlib.Path):
    owned_path = tmp_path / "owned"
    with soma.Collection.create(owned_path.as_uri()) as outer:
        # A tree of collections fully owned by us
        dog = outer.add_new_collection("dog")
        spitz = dog.add_new_collection("spitz")
        akita = spitz.add_new_collection("akita")
        hachiko = akita.add_new_dense_ndarray(
            "hachiko", type=pa.float64(), shape=(1, 2, 3)
        )
        shiba = spitz.add_new_collection("shiba")
        kabosu = shiba.add_new_sparse_ndarray("kabosu", type=pa.uint8(), shape=(10,))
        mutt = dog.add_new_collection("mutt")
        louis = mutt.add_new_dataframe(
            "louis",
            schema=pa.schema(
                (("soma_joinid", pa.int64()), ("stripes", pa.large_string()))
            ),
        )

        # A mix of collections we own and collections we don't own
        unowned_path = tmp_path / "unowned"
        unowned_path.mkdir()
        bird = outer.add_new_collection("bird")
        raptor = bird.add_new_collection("raptor")
        un_eagle = soma.Collection.create((unowned_path / "eagle").as_uri())
        raptor["eagle"] = un_eagle
        un_eagle_golden = un_eagle.add_new_collection("golden")
        un_corvid = soma.Collection.create((unowned_path / "corvid").as_uri())
        bird["corvid"] = un_corvid

        for elem in (
            dog,
            spitz,
            akita,
            hachiko,
            shiba,
            kabosu,
            mutt,
            louis,
            bird,
            raptor,
            un_eagle,
            un_eagle_golden,
            un_corvid,
        ):
            assert not elem.closed

    # Owned children should be closed
    for elem in (dog, spitz, akita, hachiko, shiba, kabosu, mutt, louis, bird, raptor):
        assert elem.closed
    # Unowned children should not be closed
    for elem in (un_eagle, un_eagle_golden, un_corvid):
        assert not elem.closed
    # Eagle cascading close
    un_eagle.close()
    assert un_eagle.closed
    assert un_eagle_golden.closed
    assert not un_corvid.closed

    # Corvid close
    un_corvid.close()
    assert un_corvid.closed

    all_elements: List[tiledb_object.AnyTileDBObject] = []

    def crawl(obj: tiledb_object.AnyTileDBObject):
        all_elements.append(obj)
        if isinstance(obj, collection.CollectionBase):
            for val in obj.values():
                crawl(val)

    with soma.Collection.open(owned_path.as_uri()) as reopened:
        # Accessing all of these results in reifying a SOMA object.
        # All of these are owned.
        crawl(reopened)
        assert len(all_elements) == 14
        assert not any(elem.closed for elem in all_elements)

        # Closing part of the subtree is fine (though not typical).
        reopened["bird"].close()
        assert reopened["bird"].closed
        assert reopened["bird"]["raptor"].closed
        # Doing so will not affect anything but that subtree.
        assert not reopened.closed
        assert not reopened["dog"].closed
    # Closing the reopened collection closes everything.
    assert all(elem.closed for elem in all_elements)


# Helper tests


@pytest.mark.parametrize(
    ("in_type", "want"),
    [
        (List[int], list),
        (set, set),
        (soma.Collection[object], soma.Collection),
    ],
)
def test_real_class(in_type, want):
    assert collection._real_class(in_type) is want


@pytest.mark.parametrize(
    "in_type", (Union[int, str], Literal["bacon"], TypeVar("_T", bound=List))
)
def test_real_class_fail(in_type):
    with pytest.raises(TypeError):
        collection._real_class(in_type)


@pytest.mark.parametrize(
    ("key", "want"),
    [
        ("hello", "hello"),
        ("good bye", "good_bye"),
        ("../beas/tie@boyz", "_beas_tie_boyz"),
        ("g0nna~let-the.BEAT", "g0nna_let_the_BEAT"),
        ("____DROP", "_DROP"),
    ],
)
def test_sanitize_for_path(key, want):
    assert collection._sanitize_for_path(key) == want


def test_timestamped_ops(tmp_path):
    """
    When we specify read/write timestamps in SOMATileDBContext supplied to collection, those are
    inherited by elements accessed via the collection.
    """

    # create collection @ t=10
    with soma.Collection.create(
        tmp_path.as_uri(), context=SOMATileDBContext(write_timestamp=10)
    ):
        pass

    # add array A to it @ t=20
    with soma.Collection.open(
        tmp_path.as_uri(), mode="w", context=SOMATileDBContext(write_timestamp=20)
    ) as sc:
        sc.add_new_dense_ndarray("A", type=pa.uint8(), shape=(2, 2)).write(
            (slice(0, 2), slice(0, 2)),
            pa.Tensor.from_numpy(np.zeros((2, 2), dtype=np.uint8)),
        )

    # access A via collection @ t=30 and write something into it
    with soma.Collection.open(
        tmp_path.as_uri(), mode="w", context=SOMATileDBContext(write_timestamp=30)
    ) as sc:
        sc["A"].write(
            (slice(0, 1), slice(0, 1)),
            pa.Tensor.from_numpy(np.ones((1, 1), dtype=np.uint8)),
        )

    # open A via collection with no timestamp => A should reflect both writes
    with soma.Collection.open(tmp_path.as_uri()) as sc:
        assert sc["A"].read((slice(None), slice(None))).to_numpy().tolist() == [
            [1, 0],
            [0, 0],
        ]

    # open A via collection @ t=25 => A should reflect first write only
    with soma.Collection.open(
        tmp_path.as_uri(), context=SOMATileDBContext(read_timestamp=25)
    ) as sc:
        assert sc["A"].read((slice(None), slice(None))).to_numpy().tolist() == [
            [0, 0],
            [0, 0],
        ]

    # open collection @ t=15 => A should not even be there
    with soma.Collection.open(
        tmp_path.as_uri(), context=SOMATileDBContext(read_timestamp=15)
    ) as sc:
        assert "A" not in sc
