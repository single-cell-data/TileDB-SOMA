import math
from typing import Any, Dict

import numpy as np
import pyarrow as pa
import pytest
from typeguard import suppress_type_checks

import tiledbsoma as soma
from tiledbsoma import _factory

from tests._util import raises_no_typeguard

""""
Metadata handling tests for all SOMA foundational datatypes.
"""


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


def test_metadata(soma_object):
    """Basic API endpoints"""
    # Verify the metadata is empty to start. "Empty" defined as no keys
    # other than soma_ keys.
    uri = soma_object.uri
    timestamp = soma_object.tiledb_timestamp_ms
    with soma_object:
        assert non_soma_metadata(soma_object) == {}
        non_soma_keys = [k for k in soma_object.metadata if not k.startswith("soma_")]
        assert non_soma_keys == []
        as_dict = dict(soma_object.metadata)
        assert frozenset(soma_object.metadata) == frozenset(as_dict)
        assert "foobar" not in soma_object.metadata

        soma_object.metadata["foobar"] = True
        assert non_soma_metadata(soma_object) == {"foobar": True}

    with pytest.raises(soma.SOMAError):
        soma_object.metadata["x"] = "y"

    with _factory.open(uri, "r", tiledb_timestamp=timestamp) as read_obj:
        assert non_soma_metadata(read_obj) == {"foobar": True}
        assert "foobar" in read_obj.metadata
        # Double-check the various getter methods
        for k, v in read_obj.metadata.items():
            assert k in read_obj.metadata
            assert read_obj.metadata.get(k) == v
            assert read_obj.metadata[k] == v
        assert read_obj.metadata.get("freeble", "blorp") == "blorp"
        with pytest.raises(soma.SOMAError):
            read_obj.metadata["x"] = "y"

    with _factory.open(uri, "w", tiledb_timestamp=timestamp + 1) as second_write:
        second_write.metadata.update(stay="frosty", my="friends")
        assert non_soma_metadata(second_write) == {
            "foobar": True,
            "stay": "frosty",
            "my": "friends",
        }

    with _factory.open(uri, "w", tiledb_timestamp=timestamp + 2) as third_write:
        del third_write.metadata["stay"]
        third_write.metadata["my"] = "enemies"
        assert non_soma_metadata(third_write) == {"foobar": True, "my": "enemies"}
        assert "stay" not in third_write.metadata
        assert third_write.metadata.get("stay", False) is False

    with _factory.open(uri, "r") as second_read:
        assert non_soma_metadata(second_read) == {"foobar": True, "my": "enemies"}
        # We don't want to test the exact metadata format,
        # just that it includes the key–value pairs.
        meta_repr = repr(second_read.metadata)
        # 'True' might get turned into '1', so only check the key.
        assert "'foobar': " in meta_repr
        assert "'my': 'enemies'" in meta_repr
    # ...but closed metadata does not.
    with suppress_type_checks():  # type checking eagerly evaluates properties including `len`, which fails on a closed object
        meta_repr_closed = repr(second_read.metadata)
    assert "foobar" not in meta_repr_closed


def test_add_delete_metadata(soma_object):
    uri = soma_object.uri
    with soma_object:
        soma_object.metadata["heres"] = "johnny"
        assert non_soma_metadata(soma_object) == {"heres": "johnny"}
        del soma_object.metadata["heres"]
        assert non_soma_metadata(soma_object) == {}

    with _factory.open(uri) as reader:
        assert non_soma_metadata(reader) == {}


def test_delete_add_metadata(soma_object):
    uri = soma_object.uri
    timestamp = soma_object.tiledb_timestamp_ms
    with soma_object:
        soma_object.metadata["hdyfn"] = "destruction"
        assert non_soma_metadata(soma_object) == {"hdyfn": "destruction"}

    with _factory.open(uri, "w", tiledb_timestamp=timestamp + 1) as second_write:
        assert non_soma_metadata(second_write) == {"hdyfn": "destruction"}
        del second_write.metadata["hdyfn"]
        assert non_soma_metadata(second_write) == {}
        second_write.metadata["hdyfn"] = "somebody new"
        assert non_soma_metadata(second_write) == {"hdyfn": "somebody new"}

    with _factory.open(uri, "r", tiledb_timestamp=timestamp + 1) as reader:
        assert non_soma_metadata(reader) == {"hdyfn": "somebody new"}


def test_set_set_metadata(soma_object):
    uri = soma_object.uri
    timestamp = soma_object.tiledb_timestamp_ms

    with soma_object:
        soma_object.metadata["content"] = "content"

    with _factory.open(uri, "w", tiledb_timestamp=timestamp + 1) as second_write:
        second_write.metadata["content"] = "confidence"
        second_write.metadata["content"] = "doubt"

    with _factory.open(uri, "r", tiledb_timestamp=timestamp + 1) as reader:
        assert non_soma_metadata(reader) == {"content": "doubt"}


def test_set_delete_metadata(soma_object):
    uri = soma_object.uri
    timestamp = soma_object.tiledb_timestamp_ms

    with soma_object:
        soma_object.metadata["possession"] = "obsession"

    with _factory.open(uri, "w", tiledb_timestamp=timestamp + 1) as second_write:
        second_write.metadata["possession"] = "funny thing about opinions"
        del second_write.metadata["possession"]

    with _factory.open(uri, "r", tiledb_timestamp=timestamp + 1) as reader:
        assert non_soma_metadata(reader) == {}


def non_soma_metadata(obj) -> Dict[str, Any]:
    return {k: v for (k, v) in obj.metadata.items() if not k.startswith("soma_")}


@pytest.mark.parametrize(
    "test_value",
    [
        True,
        False,
        0,
        1.00000001,
        -3.1415,
        "",
        "\x00",
        "\U00000000",  # get's casted to \x00
        "\x10abc",
        "\U00081a63×\x84\x94𘪩a\U000a4f44Î\x10m",
        np.str_("foo"),
        "a string",
        math.nan,
        math.inf,
        -math.inf,
    ],
)
def test_metadata_marshalling_OK(soma_object, test_value):
    """
    Test the various data type marshalling we expect to work,
    which is any Arrow primitive and Arrow strings
    """
    uri = soma_object.uri

    with soma_object:
        soma_object.metadata["test_value"] = test_value

    with _factory.open(uri, "r") as read_soma_object:
        assert "test_value" in read_soma_object.metadata

        val = read_soma_object.metadata["test_value"]
        if isinstance(test_value, float) and math.isnan(test_value):
            # By definition, NaN != NaN, so we can't just compare
            assert math.isnan(val)
        else:
            # Since an empty string is transformed to a NULL byte by numpy, passing a NULL byte is treated accordingly.
            if test_value == "\x00":
                assert val == ""
            else:
                assert val == test_value


@pytest.mark.parametrize(
    "bad_value",
    [["a", "b", "c"], {"a": False}, [1, 2, 3], np.arange(10)],
)
def test_metadata_marshalling_FAIL(soma_object, bad_value):
    """Verify that unsupported metadata types raise an error immediately."""

    with raises_no_typeguard(TypeError):
        soma_object.metadata["test_value"] = bad_value

    assert "test_value" not in soma_object.metadata


@pytest.mark.parametrize(
    "good_key",
    ["", "\x10abc", "\U00081a63×\x84\x94𘪩a\U000a4f44Î\x10m", "a string"],
)
def test_metadata_good_key(soma_object, good_key):
    """Verify that unsupported metadata types raise an error immediately."""

    soma_object.metadata[good_key] = "test_value"


@pytest.mark.parametrize(
    "bad_key",
    ["\x00", "AA\x00BB", "AA\U00000000BB"],
)
def test_metadata_bad_key(soma_object, bad_key):
    """Verify that unsupported metadata types raise an error immediately."""

    soma_object.metadata[bad_key] = "test_value"

    with pytest.raises(soma.SOMAError):
        soma_object._handle.metadata._write()


@pytest.mark.parametrize(
    "bad_value",
    [
        "AA\x00BB",
        "AA\U00000000BB",
        np.str_("foo\x00bar"),
        b"foo",
        b"\xc2",
        b"\x00",
        np.bytes_("foo"),
    ],
)
def test_metadata_bad_string_value(soma_object, bad_value):
    """Verify that unsupported metadata types raise an error immediately."""

    soma_object.metadata["test_key"] = bad_value

    with pytest.raises(soma.SOMAError):
        soma_object._handle.metadata._write()
