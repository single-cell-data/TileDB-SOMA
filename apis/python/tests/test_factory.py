from time import sleep
from typing import Type

import numpy as np
import pytest

import tiledbsoma as soma
from tiledbsoma import _constants
import tiledb

UNKNOWN_ENCODING_VERSION = "3141596"


@pytest.fixture
def tiledb_object_uri(tmp_path, object_type, metadata_typename, encoding_version):
    """Create an object with specified metadata"""
    object_uri = f"{tmp_path}/object"

    # create object
    if object_type == "array":
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
            _setmetadata(A, metadata_typename, encoding_version)
    else:
        tiledb.group_create(object_uri)
        with tiledb.Group(object_uri, mode="w") as G:
            _setmetadata(G, metadata_typename, encoding_version)

    return object_uri


@pytest.mark.parametrize(
    "object_type,metadata_typename,encoding_version,expected_soma_type",
    [
        ("group", "SOMAExperiment", _constants.SOMA_ENCODING_VERSION, soma.Experiment),
        (
            "group",
            "SOMAMeasurement",
            _constants.SOMA_ENCODING_VERSION,
            soma.Measurement,
        ),
        ("group", "SOMACollection", _constants.SOMA_ENCODING_VERSION, soma.Collection),
        ("array", "SOMADataFrame", _constants.SOMA_ENCODING_VERSION, soma.DataFrame),
        (
            "array",
            "SOMADenseNDArray",
            _constants.SOMA_ENCODING_VERSION,
            soma.DenseNDArray,
        ),
        (
            "array",
            "SOMADenseNdArray",
            _constants.SOMA_ENCODING_VERSION,
            soma.DenseNDArray,
        ),
        (
            "array",
            "SOMASparseNDArray",
            _constants.SOMA_ENCODING_VERSION,
            soma.SparseNDArray,
        ),
        (
            "array",
            "SOMASparseNdArray",
            _constants.SOMA_ENCODING_VERSION,
            soma.SparseNDArray,
        ),
    ],
)
def test_open(tiledb_object_uri, expected_soma_type: Type):
    """Happy path tests"""
    # TODO: Fix Windows test failures without the following.
    sleep(0.01)
    soma_obj = soma.open(tiledb_object_uri)
    assert isinstance(soma_obj, expected_soma_type)
    typed_soma_obj = soma.open(tiledb_object_uri, soma_type=expected_soma_type)
    assert isinstance(typed_soma_obj, expected_soma_type)
    str_typed_soma_obj = soma.open(
        tiledb_object_uri, soma_type=expected_soma_type.soma_type
    )
    assert isinstance(str_typed_soma_obj, expected_soma_type)
    assert expected_soma_type.exists(tiledb_object_uri)


@pytest.mark.parametrize(
    ("object_type", "metadata_typename", "encoding_version", "wrong_type"),
    [
        ("group", "SOMAExperiment", _constants.SOMA_ENCODING_VERSION, soma.Measurement),
        ("group", "SOMAMeasurement", _constants.SOMA_ENCODING_VERSION, soma.DataFrame),
        (
            "group",
            "SOMAMeasurement",
            _constants.SOMA_ENCODING_VERSION,
            "SOMACollection",
        ),
        (
            "array",
            "SOMADenseNDArray",
            _constants.SOMA_ENCODING_VERSION,
            soma.Collection,
        ),
        (
            "array",
            "SOMADenseNdArray",
            _constants.SOMA_ENCODING_VERSION,
            soma.SparseNDArray,
        ),
        (
            "array",
            "SOMASparseNDArray",
            _constants.SOMA_ENCODING_VERSION,
            "SOMADenseNDArray",
        ),
    ],
)
def test_open_wrong_type(tiledb_object_uri, wrong_type):
    with pytest.raises((soma.SOMAError, TypeError)):
        soma.open(tiledb_object_uri, soma_type=wrong_type)


@pytest.mark.parametrize(
    "object_type,metadata_typename,encoding_version",
    [
        ("group", "SOMAExperiment", UNKNOWN_ENCODING_VERSION),
        ("group", "SOMAMeasurement", UNKNOWN_ENCODING_VERSION),
        ("group", "SOMACollection", UNKNOWN_ENCODING_VERSION),
        ("array", "SOMADataFrame", UNKNOWN_ENCODING_VERSION),
        ("array", "SOMADenseNDArray", UNKNOWN_ENCODING_VERSION),
        ("array", "SOMASparseNDArray", UNKNOWN_ENCODING_VERSION),
    ],
)
def test_factory_unsupported_version(tiledb_object_uri):
    """All of these should raise, as they are encoding formats from the future"""
    # TODO: Fix Windows test failures without the following.
    sleep(0.01)
    with pytest.raises(ValueError):
        soma.open(tiledb_object_uri)


@pytest.mark.parametrize(
    "object_type,metadata_typename,encoding_version",
    [
        ("array", "AnUnknownTypeName", _constants.SOMA_ENCODING_VERSION),
        ("group", "AnUnknownTypeName", _constants.SOMA_ENCODING_VERSION),
        ("array", "AnUnknownTypeName", None),
        ("group", "AnUnknownTypeName", None),
        ("array", None, _constants.SOMA_ENCODING_VERSION),
        ("group", None, _constants.SOMA_ENCODING_VERSION),
        ("array", None, None),
        ("group", None, None),
        (
            "array",
            "SOMACollection",
            _constants.SOMA_ENCODING_VERSION,
        ),  # Collections can't be arrays
        (
            "group",
            "SOMADataFrame",
            _constants.SOMA_ENCODING_VERSION,
        ),  # DataFrame can't be a group
    ],
)
def test_factory_unsupported_types(tiledb_object_uri):
    """Illegal or non-existant metadata"""
    with pytest.raises(soma.SOMAError):
        soma.open(tiledb_object_uri)


def test_factory_unknown_files():
    """Test with non-TileDB files or other weirdness"""
    with pytest.raises(soma.SOMAError):
        soma.open("/tmp/no/such/file/exists/")


def _setmetadata(open_tdb_object, metadata_typename, encoding_version):
    """set only those values which are not None"""
    changes = {}
    if metadata_typename is not None:
        changes[_constants.SOMA_OBJECT_TYPE_METADATA_KEY] = metadata_typename
    if encoding_version is not None:
        changes[_constants.SOMA_ENCODING_VERSION_METADATA_KEY] = encoding_version
    if changes:
        open_tdb_object.meta.update(changes)
