"""
Basic object creation with Carrara URIs
"""

from __future__ import annotations

from collections.abc import Generator
from uuid import uuid4

import pyarrow as pa
import pytest

import tiledbsoma as soma
import tiledb

from ._util import s3_from_tiledb_uri
from .conftest import BASE_URI


@pytest.mark.carrara
def test_dataframe_create(carrara_array_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    schema = pa.schema([pa.field("soma_joinid", pa.int64(), nullable=False), ("A", pa.int32())])
    domain = ((0, 100),)
    soma.DataFrame.create(carrara_array_path, schema=schema, domain=domain, context=carrara_context).close()
    with soma.open(carrara_array_path, context=carrara_context) as A:
        assert A.soma_type == "SOMADataFrame"
        assert A.domain == domain
        assert A.schema == schema

    tbl = pa.Table.from_pydict({"soma_joinid": [0, 2], "A": [100, 200]}).cast(schema)
    with soma.open(carrara_array_path, mode="w") as A:
        A.write(tbl)

    with soma.open(carrara_array_path) as A:
        assert tbl == A.read().concat()

    with soma.open(s3_from_tiledb_uri(carrara_array_path)) as A:
        assert A.soma_type == "SOMADataFrame"
        assert A.domain == domain
        assert A.schema == schema
        assert tbl == A.read().concat()


@pytest.mark.carrara
def test_sparsendarray_create(carrara_array_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    type = pa.float32()
    shape = (999, 101)
    soma.SparseNDArray.create(carrara_array_path, type=type, shape=shape, context=carrara_context).close()
    with soma.open(carrara_array_path, context=carrara_context) as A:
        assert A.soma_type == "SOMASparseNDArray"
        assert A.shape == shape
        assert A.type == type

    with soma.open(s3_from_tiledb_uri(carrara_array_path)) as A:
        assert A.soma_type == "SOMASparseNDArray"
        assert A.shape == shape
        assert A.type == type


@pytest.mark.parametrize("soma_cls", [soma.Collection, soma.Experiment, soma.Measurement])
@pytest.mark.carrara
def test_collection_create(soma_cls, carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    soma_cls.create(carrara_group_path, context=carrara_context).close()
    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert C.soma_type == soma_cls.soma_type
        assert len(C) == 0

    with soma.open(s3_from_tiledb_uri(carrara_group_path)) as C:
        assert C.soma_type == soma_cls.soma_type
        assert len(C) == 0


@pytest.mark.carrara
def test_collection_add_new(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    schema = pa.schema([pa.field("soma_joinid", pa.int64(), nullable=False), ("A", pa.int32())])
    domain = ((0, 100),)
    type_ = pa.float32()
    shape = (99, 101)

    soma.Collection.create(carrara_group_path, context=carrara_context).close()

    with soma.open(carrara_group_path, mode="w", context=carrara_context) as C:
        C.add_new_collection("child1")
        C.add_new_collection("child2", kind=soma.Experiment)
        C.add_new_collection("child3", kind=soma.Measurement)
        C.add_new_collection("child4", kind=soma.Collection)
        C.add_new_dataframe("child5", schema=schema, domain=domain)
        C.add_new_sparse_ndarray("child6", type=type_, shape=shape)
        C.add_new_dense_ndarray("child7", type=type_, shape=shape)

    children = set(f"child{i}" for i in range(1, 8))
    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert len(C) == len(children)
        assert set(C) == children

        assert all(C._handle.is_relative(chld) for chld in children)
        assert all(C[key].uri == f"{carrara_group_path}/{key}" for key in children)

        assert C["child1"].soma_type == "SOMACollection"
        assert len(C["child1"]) == 0

        assert C["child2"].soma_type == "SOMAExperiment"
        assert len(C["child2"]) == 0

        assert C["child3"].soma_type == "SOMAMeasurement"
        assert len(C["child3"]) == 0

        assert C["child4"].soma_type == "SOMACollection"
        assert len(C["child5"]) == 0

        assert C["child5"].soma_type == "SOMADataFrame"
        assert C["child5"].schema == schema
        assert C["child5"].domain == domain

        assert C["child6"].soma_type == "SOMASparseNDArray"
        assert C["child6"].type == type_
        assert C["child6"].shape == shape

        assert C["child7"].soma_type == "SOMADenseNDArray"
        assert C["child7"].type == type_
        assert C["child7"].shape == shape

    with tiledb.Group(carrara_group_path, ctx=carrara_context.tiledb_ctx) as G:
        for chld in children:
            assert G.is_relative(chld)

    with soma.open(s3_from_tiledb_uri(carrara_group_path), context=carrara_context) as C:
        assert set(C) == children

    with soma.open(f"{s3_from_tiledb_uri(carrara_group_path)}/child5", context=carrara_context) as df:
        assert df.soma_type == "SOMADataFrame"
        assert df.schema == schema
        assert df.domain == domain


@pytest.mark.xfail(reason="Not yet working correctly on Carrara")
@pytest.mark.carrara
def test_collection_setitem(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    type_ = pa.int8()
    shape = (1, 1, 1)

    soma.Collection.create(carrara_group_path, context=carrara_context).close()

    with soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C:
        A = soma.SparseNDArray.create(f"{carrara_group_path}/sparse_ndarray", type=type_, shape=shape)
        C["sparse_ndarray"] = A

        M = soma.Measurement.create(f"{carrara_group_path}/measurement")
        C["measurement"] = M

    # the underlying asset should always have the same path, regardless of member rename
    assert soma.SparseNDArray.exists(f"{s3_from_tiledb_uri(carrara_group_path)}/sparse_ndarray")

    children = {"sparse_ndarray", "measurement"}
    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert len(C) == len(children)
        assert set(C) == children

        assert C["sparse_ndarray"].soma_type == "SOMASparseNDArray"
        assert C["sparse_ndarray"].type == type_
        assert C["sparse_ndarray"].shape == shape

        assert C["measurement"].soma_type == "SOMAMeasurement"
        assert len(C["measurement"]) == 0

    with soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C:
        C["sparse_ndarray_newname"] = C["sparse_ndarray"]  # set is rename in Carrara

    # the underlying asset should always have the same path, regardless of member rename
    assert soma.SparseNDArray.exists(f"{s3_from_tiledb_uri(carrara_group_path)}/sparse_ndarray")

    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert C["sparse_ndarray_newname"].uri == A.uri

    # Creates fail because Carrara requires uniqueness on both name and path
    with pytest.raises(soma.SOMAError):
        _ = soma.SparseNDArray.create(f"{carrara_group_path}/sparse_ndarray", type=type_, shape=shape)  # uses A's path
    with pytest.raises(soma.SOMAError):
        _ = soma.SparseNDArray.create(
            f"{carrara_group_path}/sparse_ndarray_newname", type=type_, shape=shape
        )  # uses A's member name

    # creating another object with the same path should fail
    with soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C:
        with pytest.raises(soma.SOMAError):
            # collide's with A's path name
            _ = soma.SparseNDArray.create(f"{carrara_group_path}/sparse_ndarray", type=type_, shape=shape)

        with pytest.raises(soma.SOMAError):
            # collide's with A's member name
            _ = soma.SparseNDArray.create(f"{carrara_group_path}/sparse_ndarray_newname", type=type_, shape=shape)


@pytest.mark.carrara
def test_abs_path_is_error(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    soma.Collection.create(carrara_group_path, context=carrara_context).close()

    with (
        pytest.raises((tiledb.TileDBError, soma.SOMAError), match="URI does not start with tiledb://:"),
        soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C,
        soma.open(
            "s3://cellxgene-census-public-us-west-2/cell-census/2024-07-01/soma/",
            context=soma.SOMATileDBContext({"vfs.s3.region": "us-west-2"}),
        ) as census,
    ):
        C.set("CZI_census", census, use_relative_uri=False)


@pytest.mark.xfail(reason="Not yet working correctly on Carrara")
@pytest.mark.carrara
def test_collection_set(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    type_ = pa.int8()
    shape = (1, 1, 1)

    soma.Collection.create(carrara_group_path, context=carrara_context).close()

    with soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C:
        A = soma.SparseNDArray.create(f"{carrara_group_path}/sparse_ndarray_A", type=type_, shape=shape)
        C.set("sparse_ndarray_newname", A, use_relative_uri=True)  # add member is rename in Carrara

    with soma.Collection.open(carrara_group_path, mode="r", context=carrara_context) as C:
        assert "sparse_ndarray_newname" in C
        assert C["sparse_ndarray_newname"].soma_type == "SOMASparseNDArray"
        assert C["sparse_ndarray_newname"].uri == A.uri

    with soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C:
        with pytest.raises(Exception):
            # same array with two names not allowed in Carrara
            C.set("another_name", A)

        B = soma.SparseNDArray.create(f"{carrara_group_path}/sparse_ndarray_A", type=type_, shape=shape)
        with pytest.raises(Exception):
            # absolute URLs are not allowed in Carrara
            C.set("sparse_array_B", B, use_relative_uri=False)


@pytest.mark.carrara
def test_storage_path_generates_error(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    with pytest.raises(soma.SOMAError, match="Unsupported URI format"):
        soma.Collection.create(f"{carrara_group_path}/s3://tiledb-bruce/tmp/foobar")
    with pytest.raises(soma.SOMAError, match="Unsupported URI format"):
        soma.SparseNDArray.create(
            f"{carrara_group_path}/s3://tiledb-bruce/tmp/foobar", type=pa.float32(), shape=(10, 11)
        )


URI_REASON = "CLOUD-2323: Carrara URIs must contain only unaccented alphanumeric character, '.', '-', or '_'"
NAME_REASON = (
    "CLOUD-2323: Carrara Group member names must contain only unaccented alphanumeric character, '.', '-', or '_'"
)


class TestPathEncoding:
    """
    Carrara is currently very restrictive in legal path segments and group member names:

        > name must begin with an unaccented alphanumeric character, '.', '-', or '_'

    This rules out large segments of non-Latin-1 character sets, and many punctuation supported in
    most object stores and file systems.

    Add tests, but skip them for now.
    """

    @pytest.fixture(scope="class", autouse=True)
    def a_group(self, carrara_context: soma.SOMATileDBContext) -> Generator[soma.SOMACollection, None, None, None]:
        path = f"{BASE_URI}/{uuid4()}"

        soma.Collection.create(path, context=carrara_context).close()
        yield path

        try:
            with tiledb.Group(path, mode="m") as G:
                G.delete(recursive=True)
        except tiledb.TileDBError:
            pass

    @pytest.mark.parametrize(
        "key,idx",
        [
            ("hello", 0),
            pytest.param("Γειά σου", 1, marks=pytest.mark.skip(reason=URI_REASON)),  # noqa RUF001
            pytest.param("안녕하세요", 2, marks=pytest.mark.skip(reason=URI_REASON)),
            ("h-e-l-l-o", 3),
            ("hell.o", 4),
            pytest.param("!hello", 5, marks=pytest.mark.skip(reason=URI_REASON)),
            pytest.param("[hello]", 6, marks=pytest.mark.skip(reason=URI_REASON)),
            pytest.param("hel//lo", 7, marks=pytest.mark.skip(reason=URI_REASON)),
            pytest.param("hels3://bkt/lo", 8, marks=pytest.mark.skip(reason=URI_REASON)),
            ("_", 9),
            ("-", 10),
            (".test", 11),
        ],
    )
    @pytest.mark.carrara
    def test_path_encoding(
        self, key: str, idx: int, a_group: soma.SOMACollection, carrara_context: soma.SOMATileDBContext
    ) -> None:
        with soma.Collection.open(a_group, mode="w", context=carrara_context) as C:
            C.add_new_sparse_ndarray(key, type=pa.int16(), shape=(11, 3, idx + 1))

        with soma.open(a_group, context=carrara_context) as C:
            assert C[key].shape == (11, 3, idx + 1), f"Mismatch on key={key}"


class TestMemberNameEncoding:
    """
    Verify we can use non-ASCII member names.
    """

    @pytest.fixture(scope="class", autouse=True)
    def a_group(self, carrara_context: soma.SOMATileDBContext) -> Generator[soma.SOMACollection, None, None, None]:
        path = f"{BASE_URI}/{uuid4()}"

        soma.Collection.create(path, context=carrara_context).close()
        yield path

        try:
            with tiledb.Group(path, mode="m") as G:
                G.delete(recursive=True)
        except tiledb.TileDBError:
            pass

    @pytest.mark.parametrize(
        "key,idx",
        [
            ("hello", 0),
            pytest.param("Γειά σου", 1, marks=pytest.mark.skip(reason=NAME_REASON)),  # noqa RUF001
            pytest.param("안녕하세요", 2, marks=pytest.mark.skip(reason=NAME_REASON)),
            ("h-e-l-l-o", 3),
            ("hell.o", 4),
            pytest.param("!hello", 5, marks=pytest.mark.skip(reason=NAME_REASON)),
            pytest.param("[hello]", 6, marks=pytest.mark.skip(reason=NAME_REASON)),
            pytest.param("hel//lo", 7, marks=pytest.mark.skip(reason=NAME_REASON)),
            pytest.param("hels3://bkt/lo", 8, marks=pytest.mark.skip(reason=NAME_REASON)),
            ("_", 9),
            ("-", 10),
            (".test", 11),
            pytest.param("test!", 12, marks=pytest.mark.skip(reason=NAME_REASON)),
            pytest.param("test[1]", 13, marks=pytest.mark.skip(reason=NAME_REASON)),
            pytest.param("test(2)", 14, marks=pytest.mark.skip(reason=NAME_REASON)),
        ],
    )
    @pytest.mark.carrara
    def test_member_name_encoding(
        self, key: str, idx: int, a_group: soma.SOMACollection, carrara_context: soma.SOMATileDBContext
    ) -> None:
        with soma.Collection.open(a_group, mode="w", context=carrara_context) as C:
            A = soma.SparseNDArray.create(f"{a_group}/array_{idx}", type=pa.int16(), shape=(11, 3, idx + 1))

        with soma.Collection.open(a_group, mode="w", context=carrara_context) as C:
            C[key] = A

        with soma.open(a_group, context=carrara_context) as C:
            assert key in C
            assert C[key].shape == (11, 3, idx + 1), f"Mismatch on key={key}"


@pytest.mark.xfail(reason="CLOUD-2332")
@pytest.mark.carrara
def test_name_not_eq_path(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    """Carrara group member name must be EQ to the final path segment. This restriction
    is unique to Carrara -- on all other storage subsystems, the URI and the Collection
    member name are independent.

    There are many APIs in SOMA which allow these parameters to be separately specified.
    Check that we catch them all with a reasonable error.
    """

    soma.Collection.create(carrara_group_path, context=carrara_context).close()

    with pytest.raises(soma.SOMAError), soma.open(carrara_group_path, mode="w") as C:
        C.add_new_collection(key="test1", uri=f"{carrara_group_path}/not_test1")

    with pytest.raises(soma.SOMAError), soma.open(carrara_group_path, mode="w") as C:
        C.add_new_sparse_ndarray(
            key="test2", uri=f"{carrara_group_path}/not_test2", type=pa.float64(), shape=(100, 100)
        )

    with pytest.raises(soma.SOMAError), soma.open(carrara_group_path, mode="w") as C:
        A = soma.SparseNDArray.create(f"{carrara_group_path}/test3", type=pa.int8(), shape=(100, 101))
        C["not_test3"] = A

    # Verify there were no side-effects
    with soma.open(carrara_group_path, mode="r") as C:
        assert len(C.keys()) == 0
