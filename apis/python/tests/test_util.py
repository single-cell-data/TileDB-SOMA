import pytest
from somacore import ResultOrder

from tiledbsoma._util import (
    dense_index_to_shape,
    dense_indices_to_shape,
    sanitize_key,
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
    assert uri_joinpath("tiledb://acct/s3://bucket/C", "D") == "tiledb://acct/s3://bucket/C/D"

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
    assert dense_indices_to_shape(io["coord"], io["input_shape"], io["result_order"]) == io["output_shape"]


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


@pytest.mark.parametrize(
    ("key", "sanitized"),
    (
        ("<>", "%3C%3E"),
        ("#%&*", "%23%25%26%2A"),
        ("CONFIG$", "CONFIG%24"),
        ("name_with_trailing_space_ ", "name_with_trailing_space_%20"),
        (" name_with_leading_space", "%20name_with_leading_space"),
        ("无效的文件名", "%E6%97%A0%E6%95%88%E7%9A%84%E6%96%87%E4%BB%B6%E5%90%8D"),
        (
            "path/无效的文件名",
            "path%2F%E6%97%A0%E6%95%88%E7%9A%84%E6%96%87%E4%BB%B6%E5%90%8D",
        ),
        ("path/with/ space-before-filename", "path%2Fwith%2F%20space-before-filename"),
        ("path/with/space-after-filename ", "path%2Fwith%2Fspace-after-filename%20"),
        ("%%%%%%%%%%%", "%25%25%25%25%25%25%25%25%25%25%25"),
        ("valid/path/with/slashes", "valid%2Fpath%2Fwith%2Fslashes"),
        (
            "nested/path/with_underscores/with-dashes",
            "nested%2Fpath%2Fwith_underscores%2Fwith-dashes",
        ),
        ("path/with+special-characters!", "path%2Fwith+special-characters!"),
        ("name%20with%20encoded%20spaces", "name%2520with%2520encoded%2520spaces"),
        ("name%2Fwith%2Fencoded%2Fslashes", "name%252Fwith%252Fencoded%252Fslashes"),
        (
            "path/name%20with%20encoded%20spaces",
            "path%2Fname%2520with%2520encoded%2520spaces",
        ),
        (
            "path/name%2Fwith%2Fencoded%2Fslashes",
            "path%2Fname%252Fwith%252Fencoded%252Fslashes",
        ),
        (
            "%20%20%20%20%20%20%20%20%20",
            "%2520%2520%2520%2520%2520%2520%2520%2520%2520",
        ),
        (
            "/path/with/mixed/slashes\\and\\backslashes",
            "%2Fpath%2Fwith%2Fmixed%2Fslashes%5Cand%5Cbackslashes",
        ),
        (
            "path//with///multiple////slashes",
            "path%2F%2Fwith%2F%2F%2Fmultiple%2F%2F%2F%2Fslashes",
        ),
        ("/./root_parent", "%2F.%2Froot_parent"),
        ("path/./dot_as_directory", "path%2F.%2Fdot_as_directory"),
        ("file.with..dot_segments", "file.with..dot_segments"),
        ("~/user_home_dir", "~%2Fuser_home_dir"),
        ("path.with../dot_segments/subdir", "path.with..%2Fdot_segments%2Fsubdir"),
        ("path/with/.dot-before-filename", "path%2Fwith%2F.dot-before-filename"),
        ("path/with/dot-after-filename.", "path%2Fwith%2Fdot-after-filename."),
        ("CON", "CON"),
        ("~", "~"),
    ),
)
def test_sanitize_paths(key, sanitized):
    assert sanitized == sanitize_key(key)


@pytest.mark.parametrize("key", ("..", "."))
def test_invalid_sanitize_paths(key):
    with pytest.raises(ValueError):
        assert sanitize_key(key)
