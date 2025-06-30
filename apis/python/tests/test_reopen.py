import datetime

import pyarrow as pa
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
    with create_basic_object(soma_type, uri, tiledb_timestamp=1) as x1, raises_no_typeguard(ValueError):
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
    for _mode in ["r", "w"]:
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


def test_collection_add_elements(tmp_path):
    uri = f"{tmp_path}/coll_add_and_reopen"
    with tiledbsoma.Collection.create(uri) as coll:
        assert "subcoll" not in coll
        assert "array1" not in coll
        with tiledbsoma.Collection.open(uri, mode="w") as coll2:
            coll2.add_new_collection("subcoll1")
            coll2.add_new_sparse_ndarray("array1", type=pa.float32(), shape=(100, 100))
            assert "subcoll1" in coll2
            assert "array1" in coll2

        assert "subcoll1" not in coll
        assert "array1" not in coll

        coll.reopen(mode="r")

        assert "subcoll1" in coll
        assert "array1" in coll


def test_reopen_scene_add_coord_space(tmp_path):
    uri = f"{tmp_path}/scene_reopen_add_coord_space"
    coords = tiledbsoma.CoordinateSpace([tiledbsoma.Axis("x", "meters"), tiledbsoma.Axis("y", "meters")])
    with create_basic_object("SOMAScene", uri) as scene:
        assert scene.coordinate_space is None
        with tiledbsoma.Scene.open(uri, mode="w") as scene2:
            scene2.coordinate_space = coords
            assert scene2.coordinate_space == coords
        assert scene.coordinate_space is None
        scene.reopen(mode="r")
        assert scene.coordinate_space == coords


def test_reopen_multiscale_image_add_level(tmp_path):
    uri = f"{tmp_path}/multiscal_image_reopen_add_level"
    # Need to create and close Multiscale image to fully flush creation metadata.
    image = tiledbsoma.MultiscaleImage.create(uri, type=pa.uint8(), level_shape=(3, 16, 32))
    image.close()

    with tiledbsoma.MultiscaleImage.open(uri) as image:
        assert len(image.levels()) == 1

        with tiledbsoma.MultiscaleImage.open(uri, mode="w") as image2:
            image2.add_new_level("level2", shape=(3, 8, 16))
            assert len(image2.levels()) == 2
        assert len(image.levels()) == 1
        image.reopen(mode="r")
        assert len(image.levels()) == 2
