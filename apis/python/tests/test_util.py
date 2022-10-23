import pyarrow as pa
import pytest

from tiledbsoma.util import ids_to_list, uri_joinpath


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


def test_ids_to_list():
    assert ids_to_list([]) == pa.chunked_array(pa.array([]))
    assert ids_to_list([]) == pa.chunked_array(pa.array([]))
    assert ids_to_list([1, 2, 3]) == pa.chunked_array(pa.array([1, 2, 3]))
    assert ids_to_list(pa.array([])) == pa.chunked_array(pa.array([]))
    assert ids_to_list(pa.array([1, 2, 3])) == pa.chunked_array(pa.array([1, 2, 3]))
    assert ids_to_list(slice(None)) == pa.chunked_array(pa.array([]))
    assert ids_to_list(slice(1, 5)) == pa.chunked_array(pa.array([1, 2, 3, 4, 5]))
    assert ids_to_list(slice(5, 1)) == pa.chunked_array(pa.array([5, 4, 3, 2, 1]))
