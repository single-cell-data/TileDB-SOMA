import pytest
from somacore import ResultOrder

from tiledbsoma._util import (
    dense_index_to_shape,
    dense_indices_to_shape,
    slice_to_numeric_range,
    uri_joinpath,
)


def test_uri_joinpath_file():
    assert uri_joinpath("/A/", "B") == "/A/B"
    assert uri_joinpath("/A/", "/B") == "/B"
    assert uri_joinpath("/A/B", "C") == "/A/B/C"
    assert uri_joinpath("/A/B/", "C") == "/A/B/C"
    assert uri_joinpath("/A/B/", "../C/") == "/A/B/../C"


def test_uri_joinpath_s3():
    assert uri_joinpath("s3://bucket/", "A") == "s3://bucket/A"
    assert uri_joinpath("s3://bucket/A", "B/") == "s3://bucket/A/B/"
    assert uri_joinpath("s3://bucket/A/", "B") == "s3://bucket/A/B"
    assert uri_joinpath("s3://bucket/A/", "B/") == "s3://bucket/A/B/"
    assert uri_joinpath("s3://bucket/A", "/B/") == "s3://bucket/B/"

    with pytest.raises(ValueError):
        assert uri_joinpath("s3://A/B/", "../C/")


def test_uri_joinpath_tiledb():
    assert uri_joinpath("tiledb://acct/", "A") == "tiledb://acct/A"
    assert (
        uri_joinpath("tiledb://acct/s3://bucket/C", "D")
        == "tiledb://acct/s3://bucket/C/D"
    )

    with pytest.raises(ValueError):
        assert uri_joinpath("tiledb://acct/A/", "../B/")


@pytest.mark.parametrize(
    "io",
    [
        # Note: SOMA slices are doubly inclusive
        {"coord": 1, "array_length": 10, "expected_shape": 1},
        {"coord": slice(None), "array_length": 10, "expected_shape": 10},
        {"coord": slice(0, 10), "array_length": 10, "expected_shape": 10},
        {"coord": slice(0, 9), "array_length": 10, "expected_shape": 10},
        {"coord": slice(0, 8), "array_length": 10, "expected_shape": 9},
        {"coord": slice(0, 8), "array_length": 10, "expected_shape": 9},
        {"coord": slice(3, 3), "array_length": 10, "expected_shape": 1},
        {"coord": slice(None, 3), "array_length": 10, "expected_shape": 4},
        {"coord": slice(3, None), "array_length": 10, "expected_shape": 7},
    ],
)
def test_dense_index_to_shape(io):
    assert dense_index_to_shape(io["coord"], io["array_length"]) == io["expected_shape"]


@pytest.mark.parametrize(
    "io",
    [
        # Note: SOMA slices are doubly inclusive
        {
            "coord": (1,),
            "input_shape": (10,),
            "result_order": ResultOrder.ROW_MAJOR,
            "output_shape": (1,),
        },
        {
            "coord": (1,),
            "input_shape": (10,),
            "result_order": ResultOrder.COLUMN_MAJOR,
            "output_shape": (1,),
        },
        {
            "coord": (1, 2),
            "input_shape": (10, 20, 30),
            "result_order": ResultOrder.ROW_MAJOR,
            "output_shape": (1, 1, 30),
        },
        {
            "coord": (1, 2),
            "input_shape": (10, 20, 30),
            "result_order": ResultOrder.COLUMN_MAJOR,
            "output_shape": (30, 1, 1),
        },
        {
            "coord": (1, 2, 3),
            "input_shape": (10, 20, 30),
            "result_order": ResultOrder.ROW_MAJOR,
            "output_shape": (1, 1, 1),
        },
        {
            "coord": (1, 2, 3),
            "input_shape": (10, 20, 30),
            "result_order": ResultOrder.COLUMN_MAJOR,
            "output_shape": (1, 1, 1),
        },
    ],
)
def test_dense_indices_to_shape(io):
    assert (
        dense_indices_to_shape(io["coord"], io["input_shape"], io["result_order"])
        == io["output_shape"]
    )


@pytest.mark.parametrize(
    "io",
    [
        {
            "coord": (1, 2, 3, 4),
            "input_shape": (10, 20, 30),
            "result_order": ResultOrder.ROW_MAJOR,
            "throws": ValueError,
        },
        {
            "coord": (1, 2, 3, 4),
            "input_shape": (10, 20, 30),
            "result_order": ResultOrder.COLUMN_MAJOR,
            "throws": ValueError,
        },
    ],
)
def test_dense_indices_to_shape_error_cases(io):
    with pytest.raises(io["throws"]):
        dense_indices_to_shape(io["coord"], io["input_shape"], io["result_order"])


@pytest.mark.parametrize(
    ("start_stop", "domain", "want"),
    [
        ((None, None), (-10, 10), None),
        ((5, None), (-10, 10), (5, 10)),
        ((None, 5), (-10, 10), (-10, 5)),
        ((-3, 0), (-10, 10), (-3, 0)),
        ((None, None), ("", ""), None),
    ],
)
def test_slice_to_range(start_stop, domain, want):
    slc = slice(*start_stop)
    assert want == slice_to_numeric_range(slc, domain)


@pytest.mark.parametrize(
    ("start_stop", "domain", "exc"),
    [
        ((-20, -15), (-10, 10), ValueError),
        ((15, 20), (-10, 10), ValueError),
        ((5, 10), ("", ""), TypeError),
        ((5, -5), (-10, 10), ValueError),
        (("yes", "no"), (5, 10), TypeError),
        (("START", None), ("", ""), TypeError),
    ],
)
def test_slice_to_range_bad(start_stop, domain, exc):
    with pytest.raises(exc):
        slc = slice(*start_stop)
        slice_to_numeric_range(slc, domain)
