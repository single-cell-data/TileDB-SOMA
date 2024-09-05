import datetime
import os
import pathlib
import textwrap
from typing import List, TypeVar, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from typeguard import suppress_type_checks
from typing_extensions import Literal

import tiledbsoma as soma
from tiledbsoma import _collection, _factory, _soma_group, _soma_object
from tiledbsoma._exception import DoesNotExistError, SOMAError
from tiledbsoma.options import SOMATileDBContext

from tests._util import raises_no_typeguard


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

    return _factory.open(path)


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
    return _factory.open(path)


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
        assert len(sdf.read().concat()) == 5

    with readback_collection["snda"] as snda:
        assert len(snda.read().tables().concat()) == 3


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
    assert c.members() == {"mumble": (soma_object.uri, "SOMACollection")}

    # TEMPORARY: This should no longer raise an error once TileDB supports
    # replacing an existing group member.
    with pytest.raises(soma.SOMAError):
        del c["mumble"]
        assert "mumble" not in c
        assert not c.get("mumble", False)


def test_delete_add(soma_object, tmp_path: pathlib.Path):
    tmp_uri = tmp_path.as_uri()
    with soma.Collection.create(tmp_uri) as create:
        create["porkchop sandwiches"] = soma_object

    with soma.open(tmp_uri, "w", soma_type=soma.Collection) as update:
        del update["porkchop sandwiches"]
        # TEMPORARY: This should no longer raise once TileDB supports replacing
        # an existing group member.
        with pytest.raises(soma.SOMAError):
            update["porkchop sandwiches"] = soma_object


@pytest.mark.parametrize("relative", [False, True])
def test_collection_repr(tmp_path: pathlib.Path, relative: bool) -> None:
    a_path = tmp_path / "A"
    a_uri = a_path.as_uri()
    a = soma.Collection.create(a_uri)
    assert a.uri == a_uri

    b_path = a_path / "B"
    b_uri = b_path.as_uri()
    b_uri_to_add = "B" if relative else b_uri
    b = a.add_new_collection("Another_Name", kind=soma.Experiment, uri=b_uri_to_add)

    assert b.uri == b_uri

    assert list(a.keys()) == ["Another_Name"]
    assert (
        repr(a)
        == textwrap.dedent(
            f"""
            <Collection {a_uri!r} (open for 'w') (1 item)
                'Another_Name': Experiment {b_uri!r} (open for 'w') (empty)>
            """
        ).strip()
    )
    assert a["Another_Name"] is b
    b.close()
    assert (
        repr(a)
        == textwrap.dedent(
            f"""
            <Collection {a_uri!r} (open for 'w') (1 item)
                'Another_Name': Experiment {b_uri!r} (CLOSED for 'w')>
            """
        ).strip()
    )
    a.close()
    assert repr(a) == f"<Collection {a_uri!r} (CLOSED for 'w')>"
    del a

    # re-open, reconfirm
    a_reopened = soma.Collection.open(a_uri)
    assert list(a_reopened.keys()) == ["Another_Name"]
    assert (
        repr(a_reopened)
        == textwrap.dedent(
            f"""
            <Collection {a_uri!r} (open for 'r') (1 item)
                'Another_Name': {b_uri!r} (unopened)>
            """
        ).strip()
    )
    assert a_reopened["Another_Name"].uri == b_uri
    assert (
        repr(a_reopened)
        == textwrap.dedent(
            f"""
            <Collection {a_uri!r} (open for 'r') (1 item)
                'Another_Name': Experiment {b_uri!r} (open for 'r') (empty)>
            """
        ).strip()
    )
    del a_reopened

    # move container, re-confirm wrt "relative" value
    a_moved_path = tmp_path / "A_moved"
    a_path.rename(a_moved_path)
    a_moved_uri = a_moved_path.as_uri()
    a_moved = soma.Collection.open(a_moved_uri)
    assert list(a_moved.keys()) == ["Another_Name"]
    if relative:
        new_b_uri = (a_moved_path / "B").as_uri()
        assert a_moved["Another_Name"].uri == new_b_uri
        assert (
            repr(a_moved)
            == textwrap.dedent(
                f"""
                <Collection {a_moved_uri!r} (open for 'r') (1 item)
                    'Another_Name': Experiment {new_b_uri!r} (open for 'r') (empty)>
                """
            ).strip()
        )
    else:
        with pytest.raises(DoesNotExistError):
            a_moved["Another_Name"]

        assert (
            repr(a_moved)
            == textwrap.dedent(
                f"""
                <Collection {a_moved_uri!r} (open for 'r') (1 item)
                    'Another_Name': {b_uri!r} (unopened)>
                """
            ).strip()
        )


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
    assert set(sc.keys()) == {"A"}
    assert sc["A"] == A

    with pytest.raises(SOMAError):
        sc["A"] = B
    assert set(sc.keys()) == {"A"}
    assert sc["A"] == A


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

    all_elements: List[_soma_object.AnySOMAObject] = []

    def crawl(obj: _soma_object.AnySOMAObject):
        all_elements.append(obj)
        if isinstance(obj, _collection.CollectionBase):
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
    with suppress_type_checks():
        assert _collection._real_class(in_type) is want


@pytest.mark.parametrize(
    "in_type", (Union[int, str], Literal["bacon"], TypeVar("_T", bound=List))
)
def test_real_class_fail(in_type):
    with raises_no_typeguard(TypeError):
        _collection._real_class(in_type)


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
    assert _soma_group._sanitize_for_path(key) == want


def test_timestamped_ops(tmp_path):
    """
    When we specify read/write timestamps in SOMATileDBContext supplied to collection, those are
    inherited by elements accessed via the collection.
    """

    # create collection @ t=10
    with soma.Collection.create(tmp_path.as_uri(), tiledb_timestamp=10):
        pass

    # add array A to it @ t=20
    with soma.Collection.open(tmp_path.as_uri(), mode="w", tiledb_timestamp=20) as sc:
        darr = sc.add_new_dense_ndarray("A", type=pa.uint8(), shape=(2, 2)).write(
            (slice(0, 2), slice(0, 2)),
            pa.Tensor.from_numpy(np.zeros((2, 2), dtype=np.uint8)),
        )
        assert darr.tiledb_timestamp_ms == 20

    # access A via collection @ t=30 and write something into it
    with soma.Collection.open(tmp_path.as_uri(), mode="w", tiledb_timestamp=30) as sc:
        darr = sc["A"].write(
            (slice(0, 1), slice(0, 1)),
            pa.Tensor.from_numpy(np.ones((1, 1), dtype=np.uint8)),
        )
        assert darr.tiledb_timestamp.isoformat() == "1970-01-01T00:00:00.030000+00:00"

    # open A via collection with no timestamp => A should reflect both writes
    with soma.Collection.open(tmp_path.as_uri()) as sc:
        darr = sc["A"]
        assert sc.tiledb_timestamp_ms == darr.tiledb_timestamp_ms
        assert darr.read((slice(None), slice(None))).to_numpy().tolist() == [
            [1, 0],
            [0, 0],
        ]

    # open A via collection @ t=25 => A should reflect first write only
    with soma.Collection.open(
        tmp_path.as_uri(), context=SOMATileDBContext(timestamp=25)
    ) as sc:
        assert sc["A"].read((slice(None), slice(None))).to_numpy().tolist() == [
            [0, 0],
            [0, 0],
        ]

    # open collection @ t=15 => A should not even be there
    with soma.Collection.open(
        tmp_path.as_uri(), context=SOMATileDBContext(timestamp=15)
    ) as sc:
        assert "A" not in sc

    # confirm timestamp validation in SOMATileDBContext
    with pytest.raises(ValueError):
        SOMATileDBContext(timestamp=-1)


def test_issue919(tmp_path):
    # Regression test for https://github.com/single-cell-data/TileDB-SOMA/issues/919
    # With write timestamp set, the final collection membership after adding multiple items was
    # non-deterministic, due to (i) reopening the TileDB write handle after each add operation,
    # (ii) TileDB reordering those writes randomly since they have the same timestamp, and (iii)
    # the add operations not being commutative under such randomization.
    # Fix was to eliminate (i).

    pdf = pd.DataFrame([(1, 1), (2, 2), (3, 3)], columns=["soma_joinid", "value"])
    schema = pa.Schema.from_pandas(pdf, preserve_index=False)

    for i in range(25):
        uri = str(tmp_path / str(i))

        context = SOMATileDBContext(timestamp=100)
        with soma.Collection.create(uri, context=context) as c:
            expt = c.add_new_collection("expt", soma.Experiment)
            expt.add_new_collection("causes_bug")
            expt.add_new_dataframe(
                "df", schema=schema, index_column_names=["soma_joinid"]
            )

        with soma.Collection.open(uri, context=context) as c:
            assert "df" in c["expt"] and "causes_bug" in c["expt"]
            df = c["expt"]["df"].read().concat().to_pandas()
            assert len(df) == 0


def test_context_timestamp(tmp_path: pathlib.Path):
    """Verifies that timestamps are inherited by collections."""
    fixed_time = SOMATileDBContext(timestamp=123)
    with soma.Collection.create(tmp_path.as_uri(), context=fixed_time) as coll:
        assert coll.tiledb_timestamp_ms == 123
        sub = coll.add_new_collection("sub_1")
        assert sub.tiledb_timestamp_ms == 123
        sub_sub = sub.add_new_collection("sub_sub")
        assert sub_sub.tiledb_timestamp_ms == 123

    with soma.Collection.open(tmp_path.as_uri(), context=fixed_time) as coll:
        assert coll.tiledb_timestamp_ms == 123
        sub_1 = coll["sub_1"]
        assert sub_1.tiledb_timestamp_ms == 123
        assert sub_1["sub_sub"].tiledb_timestamp_ms == 123

    with pytest.raises(soma.SOMAError):
        soma.open(tmp_path.as_uri(), context=fixed_time, tiledb_timestamp=100)

    with soma.Collection.open(
        tmp_path.as_uri(), context=fixed_time, tiledb_timestamp=234
    ) as coll:
        assert coll.tiledb_timestamp_ms == 234
        sub_1 = coll["sub_1"]
        assert sub_1.tiledb_timestamp_ms == 234
        assert sub_1["sub_sub"].tiledb_timestamp_ms == 234


def test_collection_reopen(tmp_path):
    # Ensure that reopen uses the correct mode
    soma.Collection.create(tmp_path.as_uri(), tiledb_timestamp=1)

    with soma.Collection.open(tmp_path.as_posix(), "r", tiledb_timestamp=1) as col1:
        with raises_no_typeguard(ValueError):
            col1.reopen("invalid")

        with col1.reopen("w", tiledb_timestamp=2) as col2:
            with col2.reopen("r", tiledb_timestamp=3) as col3:
                assert col1.mode == "r"
                assert col2.mode == "w"
                assert col3.mode == "r"
                assert col1.tiledb_timestamp_ms == 1
                assert col2.tiledb_timestamp_ms == 2
                assert col3.tiledb_timestamp_ms == 3

    ts1 = datetime.datetime(2023, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    ts2 = datetime.datetime(2024, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    with soma.Collection.open(tmp_path.as_posix(), "r", tiledb_timestamp=ts1) as col1:
        with col1.reopen("r", tiledb_timestamp=ts2) as col2:
            assert col1.mode == "r"
            assert col2.mode == "r"
            assert col1.tiledb_timestamp == ts1
            assert col2.tiledb_timestamp == ts2

    with soma.Collection.open(tmp_path.as_posix(), "w") as col1:
        with col1.reopen("w", tiledb_timestamp=None) as col2:
            with col2.reopen("w") as col3:
                assert col1.mode == "w"
                assert col2.mode == "w"
                assert col3.mode == "w"
                now = datetime.datetime.now(datetime.timezone.utc)
                assert col1.tiledb_timestamp <= now
                assert col2.tiledb_timestamp <= now
                assert col2.tiledb_timestamp <= now
