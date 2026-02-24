import pytest

from ._util import create_basic_object

somacore = pytest.importorskip("somacore")

# Note: using `isinstance` instead of `issubclass` for testing since `issubclass` doesn't support Protocols with properties.


def test_datafram_protocol(tmp_path):
    df = create_basic_object("SOMADataFrame", tmp_path.as_posix())
    isinstance(df, somacore.SOMAObject)
    isinstance(df, somacore.DataFrame)


def test_sparse_nd_protocol(tmp_path):
    array = create_basic_object("SOMASparseNDArray", tmp_path.as_posix())
    isinstance(array, somacore.SOMAObject)
    isinstance(array, somacore.SparseNDArray)


def test_dense_nd_protocol(tmp_path):
    array = create_basic_object("SOMADenseNDArray", tmp_path.as_posix())
    isinstance(array, somacore.SOMAObject)
    isinstance(array, somacore.DenseNDArray)


def test_collection_protocol(tmp_path):
    collection = create_basic_object("SOMACollection", tmp_path.as_posix())
    isinstance(collection, somacore.SOMAObject)
    isinstance(collection, somacore.Collection)


def test_experiment_protocol(tmp_path):
    experiment = create_basic_object("SOMAExperiment", tmp_path.as_posix())
    isinstance(experiment, somacore.SOMAObject)
    isinstance(experiment, somacore.Experiment)


def test_measurement_protocol(tmp_path):
    measurement = create_basic_object("SOMAMeasurement", tmp_path.as_posix())
    isinstance(measurement, somacore.SOMAObject)
    isinstance(measurement, somacore.Measurement)


def test_scene_protocol(tmp_path):
    scene = create_basic_object("SOMAScene", tmp_path.as_posix())
    isinstance(scene, somacore.SOMAObject)
    isinstance(scene, somacore.Scene)


def test_multiscale_image_protocol(tmp_path):
    image = create_basic_object("SOMAMultiscaleImage", tmp_path.as_posix())
    isinstance(image, somacore.SOMAObject)
    isinstance(image, somacore.MultiscaleImage)


def test_point_cloud_dataframe_protocol(tmp_path):
    dataframe = create_basic_object("SOMAPointCloudDataFrame", tmp_path.as_posix())
    isinstance(dataframe, somacore.SOMAObject)
    isinstance(dataframe, somacore.PointCloudDataFrame)


def test_geometry_dataframe_protocol(tmp_path):
    dataframe = create_basic_object("SOMAGeometryDataFrame", tmp_path.as_posix())
    isinstance(dataframe, somacore.SOMAObject)
    isinstance(dataframe, somacore.GeometryDataFrame)
