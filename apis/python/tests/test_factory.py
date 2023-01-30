from typing import Type

import numpy as np
import pytest
import tiledb

import tiledbsoma as soma
from tiledbsoma import constants, factory

UNKNOWN_ENCODING_VERSION = "3141596"


@pytest.fixture
def soma_collection(tmp_path):
    return soma.Collection.create(tmp_path.as_posix())


@pytest.fixture
def tiledb_factory(soma_collection, object_type, metadata_key, encoding_version):
    """Create a parent group and an object with specified metadata"""
    object_uri = f"{soma_collection.uri}/object"

    # create object
    if object_type == "array":
        object_cls = tiledb.Array
        schema = tiledb.ArraySchema(
            domain=tiledb.Domain(
                tiledb.Dim(name="rows", domain=(0, 100), dtype=np.int64)
            ),
            attrs=[
                tiledb.Attr(name="a", dtype=np.int32),
                tiledb.Attr(name="b", dtype=np.float32),
            ],
        )
        tiledb.Array.create(object_uri, schema)
        with tiledb.open(object_uri, mode="w") as A:
            _setmetadata(A, metadata_key, encoding_version)
    else:
        object_cls = tiledb.Group
        tiledb.group_create(object_uri)
        with tiledb.Group(object_uri, mode="w") as G:
            _setmetadata(G, metadata_key, encoding_version)

    return object_uri, soma_collection.context, object_cls


@pytest.mark.parametrize(
    "object_type,metadata_key,encoding_version,expected_soma_type",
    [
        ("group", "SOMAExperiment", constants.SOMA_ENCODING_VERSION, soma.Experiment),
        ("group", "SOMAMeasurement", constants.SOMA_ENCODING_VERSION, soma.Measurement),
        ("group", "SOMACollection", constants.SOMA_ENCODING_VERSION, soma.Collection),
        ("array", "SOMADataFrame", constants.SOMA_ENCODING_VERSION, soma.DataFrame),
        (
            "array",
            "SOMADenseNDArray",
            constants.SOMA_ENCODING_VERSION,
            soma.DenseNDArray,
        ),
        pytest.param(
            "array",
            "SOMADenseNdArray",
            constants.SOMA_ENCODING_VERSION,
            soma.DenseNDArray,
            marks=pytest.mark.xfail(reason="TileDB-SOMA bug #800"),
        ),
        ("array", "SOMASparseNDArray", "1", soma.SparseNDArray),
        pytest.param(
            "array",
            "SOMASparseNdArray",
            constants.SOMA_ENCODING_VERSION,
            soma.SparseNDArray,
            marks=pytest.mark.xfail(reason="TileDB-SOMA bug #800"),
        ),
    ],
)
def test_factory(tiledb_factory, expected_soma_type: Type):
    """Happy path tests"""
    uri, *_ = tiledb_factory
    soma_obj = factory.open(uri)
    assert isinstance(soma_obj, expected_soma_type)


@pytest.mark.parametrize(
    "object_type,metadata_key,encoding_version",
    [
        ("group", "SOMAExperiment", UNKNOWN_ENCODING_VERSION),
        ("group", "SOMAMeasurement", UNKNOWN_ENCODING_VERSION),
        ("group", "SOMACollection", UNKNOWN_ENCODING_VERSION),
        ("array", "SOMADataFrame", UNKNOWN_ENCODING_VERSION),
        ("array", "SOMADenseNDArray", UNKNOWN_ENCODING_VERSION),
        ("array", "SOMASparseNDArray", UNKNOWN_ENCODING_VERSION),
    ],
)
def test_factory_unsupported_version(tiledb_factory):
    """All of these should raise, as they are encoding formats from the future"""
    with pytest.raises(ValueError):
        factory.open(tiledb_factory[0])


@pytest.mark.parametrize(
    "object_type,metadata_key,encoding_version",
    [
        ("array", "AnUnknownTypeName", constants.SOMA_ENCODING_VERSION),
        ("group", "AnUnknownTypeName", constants.SOMA_ENCODING_VERSION),
        ("array", "AnUnknownTypeName", None),
        ("group", "AnUnknownTypeName", None),
        ("array", None, constants.SOMA_ENCODING_VERSION),
        ("group", None, constants.SOMA_ENCODING_VERSION),
        ("array", None, None),
        ("group", None, None),
        (
            "array",
            "SOMACollection",
            constants.SOMA_ENCODING_VERSION,
        ),  # Collections can't be arrays
        (
            "group",
            "SOMADataFrame",
            constants.SOMA_ENCODING_VERSION,
        ),  # DataFrame can't be a group
    ],
)
def test_factory_unsupported_types(tiledb_factory):
    """Illegal or non-existant metadata"""
    with pytest.raises(soma.SOMAError):
        factory.open(tiledb_factory[0])


def test_factory_unknown_files():
    """Test with non-TileDB files or other weirdness"""
    with pytest.raises(soma.SOMAError):
        factory.open("/tmp/no/such/file/exists/")


def _setmetadata(open_tdb_object, metadata_key, encoding_version):
    """set only those values which are not None"""
    changes = {}
    if metadata_key is not None:
        changes[constants.SOMA_OBJECT_TYPE_METADATA_KEY] = metadata_key
    if encoding_version is not None:
        changes[constants.SOMA_ENCODING_VERSION_METADATA_KEY] = encoding_version
    if changes:
        open_tdb_object.meta.update(changes)
