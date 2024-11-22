from urllib.parse import urljoin

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
import shapely
import somacore

import tiledbsoma as soma

gpd = pytest.importorskip("geopandas")
soma_outgest = pytest.importorskip("tiledbsoma.experimental.outgest")
spatialdata = pytest.importorskip("spatialdata")


@pytest.fixture(scope="module")
def sample_point_cloud_dataframe_2d(tmp_path_factory):
    baseuri = tmp_path_factory.mktemp("export_point_cloud_dataframe").as_uri()
    uri = urljoin(baseuri, "sample_2d")
    with soma.PointCloudDataFrame.create(
        uri,
        schema=pa.schema([("x", pa.float64()), ("y", pa.float64())]),
        domain=[[0, 1], [0, 1]],
    ) as point_cloud:
        x_data = np.array([0, 0, 0.5, 0.5], dtype=np.float64)
        y_data = np.array([0, 0.5, 0, 0.5], dtype=np.float64)
        data = pa.Table.from_pydict(
            {"soma_joinid": np.arange(4), "x": x_data, "y": y_data}
        )
        point_cloud.write(data)
        point_cloud.metadata["soma_geometry"] = 2.0
        point_cloud.metadata["soma_geometry_type"] = "radius"
    point_cloud = soma.PointCloudDataFrame.open(uri)
    return point_cloud


def test_export_to_shapes_2d(sample_point_cloud_dataframe_2d):
    """Test exporting a simple point cloud to a SpatialData shape model."""
    # Export PointCloudDataFrame to shapes.
    shape = soma_outgest.to_spatial_data_shapes(
        sample_point_cloud_dataframe_2d,
        scene_id="scene0",
        scene_dim_map={"x_scene": "x", "y_scene": "y"},
        soma_joinid_name="obs_id",
        transform=somacore.IdentityTransform(
            ("x_scene", "y_scene"), ("x_points", "y_points")
        ),
    )

    # Check this is valid storage for the SpatialData "Shapes" model.
    spatialdata.models.ShapesModel.validate(shape)

    # Check the dataframe.
    expected = gpd.GeoDataFrame.from_dict(
        {
            "obs_id": np.arange(4),
            "radius": 2.0 * np.ones((4,), dtype=np.float64),
            "geometry": shapely.points(
                [[0, 0], [0, 0.5], [0.5, 0], [0.5, 0.5]]
            ).tolist(),
        }
    )
    assert all(expected == shape)

    # Check the metadata.
    metadata = dict(shape.attrs)
    for key, val in metadata.items():
        print(f"{key}: {val}")
    assert len(metadata) == 1
    assert metadata["transform"] == {"scene0": spatialdata.transformations.Identity()}


def test_export_to_points_2d(sample_point_cloud_dataframe_2d):
    """Test exporting a simple point cloud to a SpatialData shape model."""
    # Export PointCloudDataFrame to shapes.
    points = soma_outgest.to_spatial_data_points(
        sample_point_cloud_dataframe_2d,
        scene_id="scene0",
        scene_dim_map={"x_scene": "x", "y_scene": "y"},
        soma_joinid_name="obs_id",
        transform=somacore.IdentityTransform(
            ("x_scene", "y_scene"), ("x_points", "y_points")
        ),
    )

    # Check this is valid storage for the SpatialData "Points" model.
    spatialdata.models.PointsModel.validate(points)

    # Check the dataframe.
    expected = pd.DataFrame.from_dict(
        {
            "x": [0, 0, 0.5, 0.5],
            "y": [0, 0.5, 0, 0.5],
            "obs_id": np.arange(4),
        }
    )
    assert all(points == expected)

    # Check the metadata.
    metadata = dict(points.attrs)
    for key, val in metadata.items():
        print(f"{key}: {val}")
    assert len(metadata) == 1
    assert metadata["transform"] == {"scene0": spatialdata.transformations.Identity()}
