import pytest

from ._util import create_basic_object

somacore = pytest.importorskip("somacore")

# Note: using `isinstance` instead of `issubclass` for testing since `issubclass` doesn't support Protocols with properties.


def test_dataframe_protocol(tmp_path):
    df = create_basic_object("SOMADataFrame", tmp_path.as_posix())
    assert isinstance(df, somacore.SOMAObject)
    assert isinstance(df, somacore.DataFrame)


def test_sparse_nd_protocol(tmp_path):
    array = create_basic_object("SOMASparseNDArray", tmp_path.as_posix())
    assert isinstance(array, somacore.SOMAObject)
    assert isinstance(array, somacore.SparseNDArray)


def test_dense_nd_protocol(tmp_path):
    array = create_basic_object("SOMADenseNDArray", tmp_path.as_posix())
    assert isinstance(array, somacore.SOMAObject)
    assert isinstance(array, somacore.DenseNDArray)


def test_collection_protocol(tmp_path):
    collection = create_basic_object("SOMACollection", tmp_path.as_posix())
    assert isinstance(collection, somacore.SOMAObject)
    assert isinstance(collection, somacore.Collection)


def test_experiment_protocol(tmp_path):
    experiment = create_basic_object("SOMAExperiment", tmp_path.as_posix())
    assert isinstance(experiment, somacore.SOMAObject)
    assert isinstance(experiment, somacore.Experiment)


def test_measurement_protocol(tmp_path):
    measurement = create_basic_object("SOMAMeasurement", tmp_path.as_posix())
    assert isinstance(measurement, somacore.SOMAObject)
    assert isinstance(measurement, somacore.Measurement)


def test_scene_protocol(tmp_path):
    scene = create_basic_object("SOMAScene", tmp_path.as_posix())
    assert isinstance(scene, somacore.SOMAObject)
    assert isinstance(scene, somacore.Scene)


def test_multiscale_image_protocol(tmp_path):
    image = create_basic_object("SOMAMultiscaleImage", tmp_path.as_posix())
    assert isinstance(image, somacore.SOMAObject)
    assert isinstance(image, somacore.MultiscaleImage)


def test_point_cloud_dataframe_protocol(tmp_path):
    dataframe = create_basic_object("SOMAPointCloudDataFrame", tmp_path.as_posix())
    assert isinstance(dataframe, somacore.SOMAObject)
    assert isinstance(dataframe, somacore.PointCloudDataFrame)


def test_geometry_dataframe_protocol(tmp_path):
    dataframe = create_basic_object("SOMAGeometryDataFrame", tmp_path.as_posix())
    assert isinstance(dataframe, somacore.SOMAObject)
    assert isinstance(dataframe, somacore.GeometryDataFrame)
