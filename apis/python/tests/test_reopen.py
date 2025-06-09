import datetime

import pytest

import tiledbsoma

from ._util import create_basic_object, raises_no_typeguard


@pytest.mark.parametrize(
    "soma_type",
    [
        "SOMAExperiment",
        "SOMAMeasurement",
        "SOMACollection",
        "SOMAScene",
        "SOMADataFrame",
        "SOMASparseNDArray",
        "SOMADenseNDArray",
        "SOMAScene",
        "SOMAPointCloudDataFrame",
        "SOMAGeometryDataFrame",
        "SOMAMultiscaleImage",
    ],
)
def test_experiment_reopen(tmp_path, soma_type):
    # Create object and an error is thrown for reopen with invalid mode.
    uri = f"{tmp_path}/{soma_type}"
    with create_basic_object(soma_type, uri, tiledb_timestamp=1) as x1:
        with raises_no_typeguard(ValueError):
            x1.reopen("invalid")

    # Check reopen without specifying timestamp.
    with tiledbsoma.open(uri, "r") as x1:
        assert x1.mode == "r"
        x2 = x1.reopen("w")
        assert x2 is x1
        assert x1.mode == "w"
    assert x1.closed

    with tiledbsoma.open(uri, "w") as x1:
        assert x1.mode == "w"
        x2 = x1.reopen("r")
        assert x2 is x1
        assert x1.mode == "r"
    assert x1.closed

    for mode in ["r", "w"]:
        with tiledbsoma.open(uri, mode) as x1:
            assert x1.mode == mode
            x2 = x1.reopen(mode)
            assert x2 is x1
            assert x1.mode == mode
    assert x1.closed

    # Check modifying the timestamp.
    with tiledbsoma.open(uri, "r") as x1:
        x2 = x1.reopen("w", tiledb_timestamp=2)
        assert x2 is x1
        assert x1.mode == "w"
        assert x1.tiledb_timestamp_ms == 2
        x2 = x1.reopen("r", tiledb_timestamp=3)
        assert x2 is x1
        assert x1.mode == "r"
        assert x1.tiledb_timestamp_ms == 3
    assert x1.closed

    ts1 = datetime.datetime(2023, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    ts2 = datetime.datetime(2024, 1, 1, 1, 0, tzinfo=datetime.timezone.utc)
    for mode in ["r", "w"]:
        with tiledbsoma.open(uri, "r", tiledb_timestamp=ts1) as x1:
            x2 = x1.reopen("r", tiledb_timestamp=ts2)
            assert x2 is x1
            assert x1.mode == "r"
            assert x1.tiledb_timestamp == ts2
        assert x1.closed

        with tiledbsoma.open(uri, "w", tiledb_timestamp=1) as x1:
            x2 = x1.reopen("w", tiledb_timestamp=None)
            assert x1 is x1
            assert x1.mode == "w"
            now = datetime.datetime.now(datetime.timezone.utc)
            assert x1.tiledb_timestamp <= now
        assert x1.closed
