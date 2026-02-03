import pytest

from tiledbsoma import ResultOrder
from tiledbsoma._util import (
    dense_index_to_shape,
    dense_indices_to_shape,
    is_relative_uri,
    make_relative_path,
    sanitize_key,
    slice_to_numeric_range,
    uri_joinpath,
)


def test_is_relative_uri() -> None:
    assert is_relative_uri("A")
    assert is_relative_uri("A/B")
    assert is_relative_uri("~")
    assert is_relative_uri("./A/B")
    assert is_relative_uri("../A/B")

    assert not is_relative_uri("file:///A")  # file: only supports absolute paths
    assert not is_relative_uri("file://./A")  # file: only supports absolute paths

    assert not is_relative_uri("s3://foo/A")
    assert not is_relative_uri("gs://foo/A")
    assert not is_relative_uri("tiledb://foo/")


def test_uri_joinpath_posix():
    # absolute base
    assert uri_joinpath("/A/", "B") == "/A/B"
    assert uri_joinpath("/A/", "/B") == "/B"
    assert uri_joinpath("/A/B", "C") == "/A/B/C"
    assert uri_joinpath("/A/B/", "C") == "/A/B/C"
    assert uri_joinpath("/A/B/", "../C/") == "/A/B/../C"

    # relative base
    assert uri_joinpath("A/", "B") == "A/B"
    assert uri_joinpath("A/", "/B") == "/B"
    assert uri_joinpath("A/B", "C") == "A/B/C"
    assert uri_joinpath("A/B/", "C") == "A/B/C"
    assert uri_joinpath("A/B/", "../C/") == "A/B/../C"


@pytest.mark.parametrize("scheme", ["s3", "gs"])
def test_uri_joinpath_object_store(scheme):
    assert uri_joinpath(f"{scheme}://bucket/", "A") == f"{scheme}://bucket/A"
    assert uri_joinpath(f"{scheme}://bucket/A", "B/") == f"{scheme}://bucket/A/B/"
    assert uri_joinpath(f"{scheme}://bucket/A/", "B") == f"{scheme}://bucket/A/B"
    assert uri_joinpath(f"{scheme}://bucket/A/", "B/") == f"{scheme}://bucket/A/B/"
    assert uri_joinpath(f"{scheme}://bucket/A", "/B/") == f"{scheme}://bucket/B/"

    with pytest.raises(ValueError):
        assert uri_joinpath(f"{scheme}://A/B/", "../C/")


def test_make_relative_path_posix():
    assert make_relative_path("/A/B/C", "/A/B/") == "C"
    assert make_relative_path("/A/B/C/", "/A/B/") == "C"
    assert make_relative_path("/A/B/C", "/A/B") == "C"
    assert make_relative_path("/A/B/C/", "/A/B") == "C"

    with pytest.raises(ValueError, match="different scheme"):
        make_relative_path("/A/B/C/", "s3://A/B/")
    with pytest.raises(ValueError, match="is not in the subpath of"):
        make_relative_path("A/B/C/", "/A/B/")


@pytest.mark.parametrize("scheme", ["s3", "gs"])
def test_make_relative_path_object_store(scheme):
    assert make_relative_path(f"{scheme}://A/B/C", f"{scheme}://A/B/") == "C"
    assert make_relative_path(f"{scheme}://A/B/C/", f"{scheme}://A/B/") == "C"
    assert make_relative_path(f"{scheme}://A/B/C", f"{scheme}://A/B") == "C"
    assert make_relative_path(f"{scheme}://A/B/C/", f"{scheme}://A/B") == "C"

    with pytest.raises(ValueError, match="different scheme"):
        make_relative_path("/A/B/C/", f"{scheme}://A/B/")
    with pytest.raises(ValueError, match="different scheme"):
        make_relative_path(f"{scheme}:/A/B/C/", "/A/B/")
    with pytest.raises(ValueError, match="is not in the subpath of"):
        make_relative_path(f"{scheme}://A/C/D/", f"{scheme}://A/B/")


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
    assert sanitized == sanitize_key(key, "tiledbv2")


@pytest.mark.parametrize("key", ("..", "."))
def test_invalid_sanitize_paths(key):
    with pytest.raises(ValueError):
        assert sanitize_key(key, "tiledbv2")
