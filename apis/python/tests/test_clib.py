import pytest

import tiledbsoma
import tiledbsoma.pytiledbsoma as somaclib

from tests._util import create_basic_object


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
        "SOMAPointCloudDataFrame",
        "SOMAGeometryDataFrame",
        "SOMAMultiscaleImage",
    ],
)
def test_get_soma_type_metadata_value(tmp_path_factory, soma_type):
    uri = str(tmp_path_factory.mktemp(f"simple_{soma_type.lower()}"))

    with create_basic_object(soma_type, uri) as soma_obj:
        soma_obj.close()

    ctx = tiledbsoma.SOMAContext()
    actual_type = somaclib.get_soma_type_metadata_value(uri, ctx.native_context, None)
    assert actual_type == soma_type


@pytest.mark.parametrize(
    "soma_type",
    [
        "SOMADataFrame",
        "SOMASparseNDArray",
        "SOMADenseNDArray",
        "SOMAPointCloudDataFrame",
        "SOMAGeometryDataFrame",
    ],
)
def test_get_soma_type_metadata_value_from_array(tmp_path_factory, soma_type):
    uri = str(tmp_path_factory.mktemp(f"simple_{soma_type.lower()}"))
    with create_basic_object(soma_type, uri) as soma_obj:
        soma_obj.close()

    ctx = tiledbsoma.SOMAContext()
    actual_type = somaclib.get_soma_type_metadata_value_from_array(uri, ctx.native_context, None)
    assert actual_type == soma_type


@pytest.mark.parametrize(
    "soma_type",
    [
        "SOMAExperiment",
        "SOMAMeasurement",
        "SOMACollection",
        "SOMAScene",
        "SOMAMultiscaleImage",
    ],
)
def test_get_soma_type_metadata_value_from_group(tmp_path_factory, soma_type):
    uri = str(tmp_path_factory.mktemp(f"simple_{soma_type.lower()}"))

    with create_basic_object(soma_type, uri) as soma_obj:
        soma_obj.close()

    ctx = tiledbsoma.SOMAContext()
    actual_type = somaclib.get_soma_type_metadata_value_from_group(uri, ctx.native_context, None)
    assert actual_type == soma_type


def test_get_soma_type_metadata_value_from_array_error(tmp_path_factory):
    uri = str(tmp_path_factory.mktemp("test_collection"))
    with create_basic_object("SOMACollection", uri) as coll:
        coll.close()

    ctx = tiledbsoma.SOMAContext()
    with pytest.raises(Exception, match=r".* Array does not exist."):
        somaclib.get_soma_type_metadata_value_from_array(uri, ctx.native_context, None)


def test_get_soma_type_metadata_value_from_group_error(tmp_path_factory):
    uri = str(tmp_path_factory.mktemp("test_dataframe"))
    with create_basic_object("SOMADataFrame", uri) as df:
        df.close()

    ctx = tiledbsoma.SOMAContext()
    with pytest.raises(Exception, match=r".* Group does not exist."):
        somaclib.get_soma_type_metadata_value_from_group(uri, ctx.native_context, None)
