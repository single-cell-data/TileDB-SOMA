import math

import pyarrow as pa
import pytest

import tiledbsoma as soma

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
    uri = tmp_path.as_posix()
    class_name = request.param

    if class_name == "Collection":
        so = soma.Collection(uri=uri)
        so.create()

    elif class_name == "DataFrame":
        so = soma.DataFrame(uri=uri)
        so.create(
            schema=pa.schema([("C", pa.float32()), ("D", pa.uint32())]),
            index_column_names=["D"],
        )

    elif class_name == "DenseNDArray":
        so = soma.DenseNDArray(uri=uri)
        so.create(type=pa.float64(), shape=(100, 10, 1))

    elif class_name == "SparseNDArray":
        so = soma.SparseNDArray(uri=uri)
        so.create(type=pa.int8(), shape=(11,))

    assert so is not None, f"Unknown class name: {class_name}"
    yield so
    so.delete()


def test_metadata(soma_object):
    """Basic API endpoints"""
    assert soma_object.exists()

    # Verify the metadata is empty to start. "Empty" defined as no keys
    # other than soma_ keys.
    keys = list(
        filter(lambda s: not s.startswith("soma_"), soma_object.metadata().keys())
    )
    assert keys == []
    assert len(soma_object.metadata()) == len(soma_object.metadata().keys())
    assert len(soma_object.metadata()) == len(soma_object.metadata().as_dict())
    assert "foobar" not in soma_object.metadata()

    soma_object.metadata()["foobar"] = True
    assert "foobar" in soma_object.metadata()
    for k, v in soma_object.metadata().as_dict().items():
        assert k in soma_object.metadata()
        assert soma_object.metadata().get(k) == v
        assert soma_object.metadata()[k] == v

    # also check set()
    soma_object.metadata()["stay"] = "frosty"
    assert "stay" in soma_object.metadata()
    assert soma_object.metadata()["stay"] == "frosty"

    del soma_object.metadata()["stay"]
    assert "stay" not in soma_object.metadata()
    assert soma_object.metadata().get("stay", False) is False


@pytest.mark.parametrize(
    "test_value",
    [
        True,
        False,
        0,
        1.00000001,
        -3.1415,
        "",
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
    soma_object.metadata()["test_value"] = test_value
    assert "test_value" in soma_object.metadata()

    val = soma_object.metadata()["test_value"]
    if type(test_value) is float and math.isnan(test_value):
        # By definition, NaN != NaN, so we can't just compare
        assert math.isnan(val)
    else:
        assert val == test_value


@pytest.mark.parametrize(
    "test_value",
    [["a", "b", "c"], {"a": False}],
)
def test_metadata_marshalling_FAIL(soma_object, test_value):
    """Test the various data type marshalling we expect to FAIL"""

    with pytest.raises(TypeError):
        soma_object.metadata()["test_value"] = test_value

    assert "test_value" not in soma_object.metadata()
