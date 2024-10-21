from time import sleep
from typing import Type

import numpy as np
import pyarrow as pa
import pytest

import tiledbsoma as soma
from tiledbsoma import _constants

UNKNOWN_ENCODING_VERSION = "3141596"


@pytest.fixture
def tiledb_object_uri(tmp_path, metadata_typename, encoding_version, soma_type):
    """Create an object with specified metadata"""
    object_uri = f"{tmp_path}/object"
    kwargs = {}

    if issubclass(soma_type, (soma.DenseNDArray, soma.SparseNDArray)):
        kwargs["type"] = pa.int64()
        kwargs["shape"] = (100,)
    elif issubclass(soma_type, soma.DataFrame):
        kwargs["schema"] = pa.schema(
            [("rows", pa.int64()), ("a", pa.int32()), ("b", pa.float32())]
        )

    soma_type.create(object_uri, tiledb_timestamp=1, **kwargs).close()

    with soma_type.open(object_uri, "w", tiledb_timestamp=2) as soma_obj:
        _setmetadata(soma_obj, metadata_typename, encoding_version)

    return object_uri


@pytest.mark.parametrize(
    "metadata_typename, encoding_version, soma_type",
    [
        ("SOMAExperiment", _constants.SOMA_ENCODING_VERSION, soma.Experiment),
        ("SOMAMeasurement", _constants.SOMA_ENCODING_VERSION, soma.Measurement),
        ("SOMACollection", _constants.SOMA_ENCODING_VERSION, soma.Collection),
        ("SOMADataFrame", _constants.SOMA_ENCODING_VERSION, soma.DataFrame),
        ("SOMADenseNDArray", _constants.SOMA_ENCODING_VERSION, soma.DenseNDArray),
        ("SOMADenseNdArray", _constants.SOMA_ENCODING_VERSION, soma.DenseNDArray),
        ("SOMASparseNDArray", _constants.SOMA_ENCODING_VERSION, soma.SparseNDArray),
        ("SOMASparseNdArray", _constants.SOMA_ENCODING_VERSION, soma.SparseNDArray),
    ],
)
def test_open(tiledb_object_uri, soma_type: Type):
    """Happy path tests"""
    # TODO: Fix Windows test failures without the following.
    sleep(0.01)
    soma_obj = soma.open(tiledb_object_uri, tiledb_timestamp=2)
    assert isinstance(soma_obj, soma_type)
    typed_soma_obj = soma.open(
        tiledb_object_uri, soma_type=soma_type, tiledb_timestamp=2
    )
    assert isinstance(typed_soma_obj, soma_type)
    str_typed_soma_obj = soma.open(
        tiledb_object_uri, soma_type=soma_type.soma_type, tiledb_timestamp=2
    )
    assert isinstance(str_typed_soma_obj, soma_type)
    assert soma_type.exists(tiledb_object_uri)


@pytest.mark.parametrize(
    ("metadata_typename", "encoding_version", "soma_type"),
    [
        ("SOMAExperiment", _constants.SOMA_ENCODING_VERSION, soma.Measurement),
        ("SOMAMeasurement", _constants.SOMA_ENCODING_VERSION, soma.DataFrame),
        ("SOMAMeasurement", _constants.SOMA_ENCODING_VERSION, soma.Collection),
        ("SOMADenseNDArray", _constants.SOMA_ENCODING_VERSION, soma.Collection),
        ("SOMADenseNdArray", _constants.SOMA_ENCODING_VERSION, soma.SparseNDArray),
        ("SOMASparseNDArray", _constants.SOMA_ENCODING_VERSION, soma.DenseNDArray),
    ],
)
def test_open_wrong_type(tiledb_object_uri, soma_type):
    with pytest.raises((soma.SOMAError, TypeError)):
        soma.open(tiledb_object_uri, soma_type=soma_type, tiledb_timestamp=2)


@pytest.mark.parametrize(
    "metadata_typename, encoding_version, soma_type",
    [
        ("SOMAExperiment", UNKNOWN_ENCODING_VERSION, soma.Experiment),
        ("SOMAMeasurement", UNKNOWN_ENCODING_VERSION, soma.Measurement),
        ("SOMACollection", UNKNOWN_ENCODING_VERSION, soma.Collection),
        ("SOMADataFrame", UNKNOWN_ENCODING_VERSION, soma.DataFrame),
        ("SOMADenseNDArray", UNKNOWN_ENCODING_VERSION, soma.DenseNDArray),
        ("SOMASparseNDArray", UNKNOWN_ENCODING_VERSION, soma.SparseNDArray),
    ],
)
def test_factory_unsupported_version(tiledb_object_uri):
    """All of these should raise, as they are encoding formats from the future"""
    # TODO: Fix Windows test failures without the following.
    sleep(0.01)
    with pytest.raises(ValueError):
        soma.open(tiledb_object_uri, tiledb_timestamp=2)


@pytest.mark.parametrize(
    "metadata_typename, encoding_version, soma_type",
    [
        (
            "AnUnknownTypeName",
            _constants.SOMA_ENCODING_VERSION,
            soma.DataFrame,
        ),  # Invalid type
        (
            "AnUnknownTypeName",
            _constants.SOMA_ENCODING_VERSION,
            soma.Collection,
        ),  # Invalid type
        ("AnUnknownTypeName", None, soma.DataFrame),  # Invalid type and no version
        ("AnUnknownTypeName", None, soma.Collection),  # Invalid type and no version
        (None, _constants.SOMA_ENCODING_VERSION, soma.DataFrame),  # No type given
        (None, _constants.SOMA_ENCODING_VERSION, soma.Collection),  # No type give
        (None, None, soma.DataFrame),  # Neither type nor version filled
        (None, None, soma.Collection),  # Neither type nor version filled
        (
            "SOMACollection",
            _constants.SOMA_ENCODING_VERSION,
            soma.DataFrame,
        ),  # Collections can't be an array
        (
            "SOMADataFrame",
            _constants.SOMA_ENCODING_VERSION,
            soma.Collection,
        ),  # DataFrame can't be a group
    ],
)
def test_factory_unsupported_types(tiledb_object_uri):
    """Illegal or non-existant metadata"""
    with pytest.raises(soma.SOMAError):
        soma.open(tiledb_object_uri, tiledb_timestamp=2)


def test_factory_unknown_files():
    """Test with non-TileDB files or other weirdness"""
    with pytest.raises(soma.SOMAError):
        soma.open("/tmp/no/such/file/exists/")


def _setmetadata(open_tdb_object, metadata_typename, encoding_version):
    """force modify the metadata values"""
    set_metadata = open_tdb_object._handle._handle.set_metadata
    del_metadata = open_tdb_object._handle._handle.delete_metadata

    if metadata_typename is not None:
        val = np.array([metadata_typename], "S")
        set_metadata(_constants.SOMA_OBJECT_TYPE_METADATA_KEY, val, True)
    else:
        del_metadata(_constants.SOMA_OBJECT_TYPE_METADATA_KEY, True)

    if encoding_version is not None:
        val = np.array([encoding_version], "S")
        set_metadata(_constants.SOMA_ENCODING_VERSION_METADATA_KEY, val, True)
    else:
        del_metadata(_constants.SOMA_ENCODING_VERSION_METADATA_KEY, True)
