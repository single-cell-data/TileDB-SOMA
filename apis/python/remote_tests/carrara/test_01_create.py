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


@pytest.mark.carrara
def test_sparsendarray_create(carrara_array_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    type = pa.float32()
    shape = (999, 101)
    soma.SparseNDArray.create(carrara_array_path, type=type, shape=shape, context=carrara_context).close()
    with soma.open(carrara_array_path, context=carrara_context) as A:
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


@pytest.mark.carrara
def test_collection_set(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    type_ = pa.int8()
    shape = (1, 1, 1)

    soma.Collection.create(carrara_group_path, context=carrara_context).close()

    with soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C:
        A = soma.SparseNDArray.create(f"{carrara_group_path}/sparse_ndarray", type=type_, shape=shape)
        C["sparse_ndarray"] = A

        M = soma.Measurement.create(f"{carrara_group_path}/measurement")
        C["measurement"] = M

    children = {"sparse_ndarray", "measurement"}
    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert len(C) == len(children)
        assert set(C) == children

        assert C["sparse_ndarray"].soma_type == "SOMASparseNDArray"
        assert C["sparse_ndarray"].type == type_
        assert C["sparse_ndarray"].shape == shape

        assert C["measurement"].soma_type == "SOMAMeasurement"
        assert len(C["measurement"]) == 0


@pytest.mark.carrara
def test_storage_path_generates_error(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    with pytest.raises(soma.SOMAError, match="Unsupported URI format"):
        soma.Collection.create(f"{carrara_group_path}/s3://tiledb-bruce/tmp/foobar")
    with pytest.raises(soma.SOMAError, match="Unsupported URI format"):
        soma.SparseNDArray.create(
            f"{carrara_group_path}/s3://tiledb-bruce/tmp/foobar", type=pa.float32(), shape=(10, 11)
        )


NAME_REASON = "CLOUD-2323: Carrara names must contain only unaccented alphanumeric character, '.', '-', or '_'"


class TestPathEncoding:
    """
    Carrara is currently very restritive in legal path segments and group member names:

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


@pytest.mark.xfail(reason="CLOUD-2332")
def test_name_noteq_path(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
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
