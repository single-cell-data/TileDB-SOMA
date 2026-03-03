from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
import shapely
import typeguard

import tiledbsoma as soma
from tiledbsoma._core_options import BatchSize


def test_point_cloud_bad_create(tmp_path):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "bad_create")

    # all spatial axis must have the same type
    asch = pa.schema([("x", pa.float64()), ("y", pa.int64())])
    with pytest.raises(ValueError):
        soma.PointCloudDataFrame.create(
            urljoin(baseuri, "different_types"),
            schema=asch,
            domain=[[0, 9]],
        )

    # type must be integral or floating-point
    asch = pa.schema([("x", pa.large_string()), ("y", pa.large_string())])
    with pytest.raises(ValueError):
        soma.PointCloudDataFrame.create(urljoin(baseuri, "bad_type"), schema=asch, domain=[[0, 9]])


def test_point_cloud_basic_read(tmp_path):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "basic_read")

    asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])

    # defaults
    with soma.PointCloudDataFrame.create(
        urljoin(baseuri, "default"),
        schema=asch,
        domain=[[-10000, 10000], [-10000, 10000], [0, 10]],
    ) as ptc:
        pydict = {}
        pydict["soma_joinid"] = [1, 2, 3, 4, 5]
        pydict["x"] = [10, 20, 30, 40, 50]
        pydict["y"] = [4.1, 5.2, 6.3, 7.4, 8.5]

        rb = pa.Table.from_pydict(pydict, schema=asch.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    with soma.PointCloudDataFrame.open(urljoin(baseuri, "default"), "r") as ptc:
        assert set(ptc.schema.names) == {"soma_joinid", "x", "y"}
        assert ptc.index_column_names == ("x", "y", "soma_joinid")
        assert ptc.axis_names == ("x", "y")

        table = ptc.read().concat()
        assert ptc.count == len(ptc) == table.num_rows == 5
        assert table.num_columns == 3
        assert [e.as_py() for e in table["soma_joinid"]] == pydict["soma_joinid"]
        assert [e.as_py() for e in table["x"]] == pydict["x"]
        assert [e.as_py() for e in table["y"]] == pydict["y"]

    # with user defined values
    with soma.PointCloudDataFrame.create(
        urljoin(baseuri, "user_defined"),
        schema=asch,
        coordinate_space="x",
        domain=((1, 10), (0, 10)),
    ) as ptc:
        pydict = {}
        pydict["soma_joinid"] = [1, 2, 3, 4, 5]
        pydict["x"] = [1, 2, 3, 4, 5]
        pydict["y"] = [4.1, 5.2, 6.3, 7.4, 8.5]

        rb = pa.Table.from_pydict(pydict, schema=asch.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    with soma.PointCloudDataFrame.open(urljoin(baseuri, "user_defined"), "r") as ptc:
        assert set(ptc.schema.names) == set(["soma_joinid", "x", "y"])
        assert ptc.index_column_names == ("x", "soma_joinid")
        assert ptc.axis_names == ("x",)
        assert ptc.domain == ((1, 10), (0, 10))

        table = ptc.read().concat()
        assert ptc.count == len(ptc) == table.num_rows == 5
        assert table.num_columns == 3
        assert [e.as_py() for e in table["soma_joinid"]] == pydict["soma_joinid"]
        assert [e.as_py() for e in table["x"]] == pydict["x"]
        assert [e.as_py() for e in table["y"]] == pydict["y"]

        assert (
            ptc.metadata[soma._constants.SOMA_SPATIAL_VERSION_METADATA_KEY]
            == soma._constants.SOMA_SPATIAL_ENCODING_VERSION
        )

    # Ensure it cannot be opened by another type
    with pytest.raises(soma.SOMAError):
        soma.DataFrame.open(urljoin(baseuri, "default"))

    with pytest.raises(soma.SOMAError):
        soma.SparseNDArray.open(urljoin(baseuri, "default"))

    with pytest.raises(soma.SOMAError):
        soma.DenseNDArray.open(urljoin(baseuri, "default"))

    with pytest.raises(soma.SOMAError):
        soma.Collection.open(urljoin(baseuri, "default"))

    with pytest.raises(soma.SOMAError):
        soma.Experiment.open(urljoin(baseuri, "default"))

    with pytest.raises(soma.SOMAError):
        soma.Measurement.open(urljoin(baseuri, "default"))

    with pytest.raises(soma.SOMAError):
        soma.Scene.open(urljoin(baseuri, "default"))

    with pytest.raises(soma.SOMAError):
        soma.MultiscaleImage.open(urljoin(baseuri, "default"))


def test_point_cloud_coordinate_space(tmp_path):
    uri = tmp_path.as_uri()

    asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])

    with soma.PointCloudDataFrame.create(uri, schema=asch, domain=((-10, 10), (-10, 10), (0, 1000))) as ptc:
        assert len(ptc.coordinate_space) == 2
        assert ptc.coordinate_space.axis_names == ("x", "y")
        assert ptc.coordinate_space.axes == (soma.Axis(name="x"), soma.Axis(name="y"))

        # Axis names do not match
        with pytest.raises(ValueError):
            ptc.coordinate_space = soma.CoordinateSpace([soma.Axis(name="a"), soma.Axis(name="y")])

        ptc.coordinate_space = soma.CoordinateSpace([soma.Axis(name="x", unit="m"), soma.Axis(name="y", unit="in")])
        assert ptc.coordinate_space[0] == soma.Axis(name="x", unit="m")
        assert ptc.coordinate_space[1] == soma.Axis(name="y", unit="in")


def test_point_cloud_bad_read_spatial_region(tmp_path):
    uri = tmp_path.as_uri()

    schema = pa.schema([("x", pa.float64()), ("y", pa.float64())])

    with soma.PointCloudDataFrame.create(uri, schema=schema, domain=[[-10000, 10000], [-10000, 10000], [0, 10]]) as ptc:
        pydict = {
            "soma_joinid": [1, 2, 3, 4, 5],
            "x": [10, 20, 30, 40, 50],
            "y": [4.1, 5.2, 6.3, 7.4, 8.5],
        }
        rb = pa.Table.from_pydict(pydict, schema=schema.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    with soma.PointCloudDataFrame.open(uri, "r") as ptc:
        # Cannot specify the output coordinate space when transform is None
        with pytest.raises(ValueError):
            ptc.read_spatial_region(
                region_transform=None,
                region_coord_space=soma.CoordinateSpace((soma.Axis(name="x"), soma.Axis(name="y"))),
            )

        # The input axes of the transform must match the axes of the coordinate
        # space the requested region is defined in
        with pytest.raises(ValueError):
            ptc.read_spatial_region(
                region_transform=soma.IdentityTransform(input_axes=("x", "y"), output_axes=("x", "y")),
                region_coord_space=soma.CoordinateSpace([soma.Axis(name="y"), soma.Axis(name="x")]),
            )


@pytest.mark.parametrize(
    "name, region, value_filter, expected_output",
    [
        (
            "Bounding box",
            [10, 4.1, 50, 8.5],
            None,
            {
                "soma_joinid": [1, 2, 3, 4, 5],
                "x": [10, 20, 30, 40, 50],
                "y": [4.1, 5.2, 6.3, 7.4, 8.5],
            },
        ),
        (
            "Full region",
            [0, 0, 100, 100],
            None,
            {
                "soma_joinid": [1, 2, 3, 4, 5],
                "x": [10, 20, 30, 40, 50],
                "y": [4.1, 5.2, 6.3, 7.4, 8.5],
            },
        ),
        (
            "Partial region",
            [15, 4, 35, 7],
            None,
            {
                "soma_joinid": [2, 3],
                "x": [20, 30],
                "y": [5.2, 6.3],
            },
        ),
        (
            "Single point",
            [10, 4.1, 10, 4.1],
            None,
            {
                "soma_joinid": [1],
                "x": [10],
                "y": [4.1],
            },
        ),
        (
            "No points",
            [60, 10, 80, 20],
            None,
            {
                "soma_joinid": [],
                "x": [],
                "y": [],
            },
        ),
        (
            "Flipped box",
            [50, 10, 20, 0],
            None,
            {
                "soma_joinid": [2, 3, 4, 5],
                "x": [20, 30, 40, 50],
                "y": [5.2, 6.3, 7.4, 8.5],
            },
        ),
        (
            "Bounding box with value filter",
            [10, 4.1, 50, 8.5],
            "soma_joinid in [1, 2, 3, 4, 5]",
            {
                "soma_joinid": [1, 2, 3, 4, 5],
                "x": [10, 20, 30, 40, 50],
                "y": [4.1, 5.2, 6.3, 7.4, 8.5],
            },
        ),
        (
            "Partial region with value filter",
            [15, 4, 35, 7],
            "soma_joinid in [2, 3]",
            {
                "soma_joinid": [2, 3],
                "x": [20, 30],
                "y": [5.2, 6.3],
            },
        ),
        (
            "Single point with value filter",
            [10, 4.1, 10, 4.1],
            "soma_joinid == 1",
            {
                "soma_joinid": [1],
                "x": [10],
                "y": [4.1],
            },
        ),
        (
            "No points with value filter",
            [60, 10, 80, 20],
            "soma_joinid in [1, 2, 3, 4, 5]",
            {
                "soma_joinid": [],
                "x": [],
                "y": [],
            },
        ),
        (
            "Flipped box with value filter",
            [50, 10, 20, 0],
            "soma_joinid in [2, 3, 4, 5]",
            {
                "soma_joinid": [2, 3, 4, 5],
                "x": [20, 30, 40, 50],
                "y": [5.2, 6.3, 7.4, 8.5],
            },
        ),
    ],
)
def test_point_cloud_read_spatial_region_basic_2d(tmp_path, name, region, value_filter, expected_output):
    uri = tmp_path.as_uri()

    schema = pa.schema([("x", pa.float64()), ("y", pa.float64())])

    with soma.PointCloudDataFrame.create(uri, schema=schema, domain=[[-10000, 10000], [-10000, 10000], [0, 10]]) as ptc:
        pydict = {
            "soma_joinid": [1, 2, 3, 4, 5],
            "x": [10, 20, 30, 40, 50],
            "y": [4.1, 5.2, 6.3, 7.4, 8.5],
        }
        rb = pa.Table.from_pydict(pydict, schema=schema.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    with soma.PointCloudDataFrame.open(uri, "r") as ptc:
        actual_output = ptc.read_spatial_region(region=region, value_filter=value_filter)
        assert actual_output.data.concat().to_pydict() == expected_output

        actual_output = ptc.read_spatial_region(region=shapely.box(*region), value_filter=value_filter)
        assert actual_output.data.concat().to_pydict() == expected_output

        actual_output = ptc.read_spatial_region(region=region, column_names=["x"], value_filter=value_filter)
        if value_filter:
            assert actual_output.data.concat().to_pydict() == {
                "soma_joinid": expected_output["soma_joinid"],
                "x": expected_output["x"],
            }
        else:
            assert actual_output.data.concat().to_pydict() == {
                "x": expected_output["x"],
            }

        actual_output = ptc.read_spatial_region(region=region, column_names=["y"], value_filter=value_filter)
        if value_filter:
            assert actual_output.data.concat().to_pydict() == {
                "soma_joinid": expected_output["soma_joinid"],
                "y": expected_output["y"],
            }
        else:
            assert actual_output.data.concat().to_pydict() == {
                "y": expected_output["y"],
            }


@pytest.mark.skip("3D regions not supported yet")
@pytest.mark.parametrize(
    "name, region, expected_output",
    [
        (
            "Bounding box",
            [10, 4.1, 50, 8.5, 0, 100],
            {
                "soma_joinid": [1, 2, 3, 4, 5],
                "x": [10, 20, 30, 40, 50],
                "y": [4.1, 5.2, 6.3, 7.4, 8.5],
                "z": [0, 1, 2, 3, 4],
            },
        ),
        (
            "Full region",
            [0, 0, 100, 100, 0, 10],
            {
                "soma_joinid": [1, 2, 3, 4, 5],
                "x": [10, 20, 30, 40, 50],
                "y": [4.1, 5.2, 6.3, 7.4, 8.5],
                "z": [0, 1, 2, 3, 4],
            },
        ),
        (
            "Partial region",
            [15, 4, 35, 7, 0, 5],
            {
                "soma_joinid": [2, 3],
                "x": [20, 30],
                "y": [5.2, 6.3],
                "z": [1, 2],
            },
        ),
        (
            "Single point",
            [10, 4.1, 10, 4.1, 0, 0],
            {
                "soma_joinid": [1],
                "x": [10],
                "y": [4.1],
                "z": [0],
            },
        ),
        (
            "No points",
            [60, 10, 80, 20, 0, 1],
            {
                "soma_joinid": [],
                "x": [],
                "y": [],
                "z": [],
            },
        ),
        (
            "Flipped box",
            [50, 10, 20, 0, 0, 5],
            {
                "soma_joinid": [2, 3, 4, 5],
                "x": [20, 30, 40, 50],
                "y": [5.2, 6.3, 7.4, 8.5],
                "z": [1, 2, 3, 4],
            },
        ),
    ],
)
def test_point_cloud_read_spatial_region_basic_3d(tmp_path, name, region, expected_output):
    uri = tmp_path.as_uri()

    schema = pa.schema([("x", pa.float64()), ("y", pa.float64()), ("z", pa.float64())])

    with soma.PointCloudDataFrame.create(
        uri,
        schema=schema,
        coordinate_space=("x", "y", "z"),
        domain=[[-10000, 10000], [-10000, 10000], [-10000, 10000], [0, 10]],
    ) as ptc:
        pydict = {
            "soma_joinid": [1, 2, 3, 4, 5],
            "x": [10, 20, 30, 40, 50],
            "y": [4.1, 5.2, 6.3, 7.4, 8.5],
            "z": [0, 1, 2, 3, 4],
        }
        rb = pa.Table.from_pydict(pydict, schema=schema.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    with soma.PointCloudDataFrame.open(uri, "r") as ptc:
        actual_output = ptc.read_spatial_region(region=region)
        assert actual_output.data.concat().to_pydict() == expected_output


@pytest.mark.parametrize(
    "name, region, exc_type",
    [
        ("Empty", [], ValueError),
        ("Too few", [10, 4.1], ValueError),
        ("Too many", [0, 0, 100, 100, 50], ValueError),
        ("Bad type", ["a", "b", "c", "d"], typeguard.TypeCheckError),
    ],
)
def test_point_cloud_read_spatial_region_2d_bad(tmp_path, name, region, exc_type):
    uri = tmp_path.as_uri()

    schema = pa.schema([("x", pa.float64()), ("y", pa.float64())])

    with soma.PointCloudDataFrame.create(uri, schema=schema, domain=[[-10000, 10000], [-10000, 10000], [0, 10]]) as ptc:
        pydict = {
            "soma_joinid": [1, 2, 3, 4, 5],
            "x": [10, 20, 30, 40, 50],
            "y": [4.1, 5.2, 6.3, 7.4, 8.5],
        }
        rb = pa.Table.from_pydict(pydict, schema=schema.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    with soma.PointCloudDataFrame.open(uri, "r") as ptc, pytest.raises(exc_type):
        ptc.read_spatial_region(region=region)


@pytest.mark.skip("3D regions not supported yet")
@pytest.mark.parametrize(
    "name, region, exc_type",
    [
        ("Empty", [], ValueError),
        ("Too few", [10, 4.1, 20, 6.3], ValueError),
        ("Too many", [0, 0, 100, 100, 50, 10, 0], ValueError),
        ("Bad type", ["a", "b", "c", "d", "e", "f"], typeguard.TypeCheckError),
    ],
)
def test_point_cloud_read_spatial_region_3d_bad(tmp_path, name, region, exc_type):
    uri = tmp_path.as_uri()

    schema = pa.schema([("x", pa.float64()), ("y", pa.float64()), ("z", pa.float64())])

    with soma.PointCloudDataFrame.create(uri, schema=schema, domain=[[0, 9], [0, 10]]) as ptc:
        pydict = {
            "soma_joinid": [1, 2, 3, 4, 5],
            "x": [10, 20, 30, 40, 50],
            "y": [4.1, 5.2, 6.3, 7.4, 8.5],
            "z": [0, 1, 2, 3, 4],
        }
        rb = pa.Table.from_pydict(pydict, schema=schema.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    with soma.PointCloudDataFrame.open(uri, "r") as ptc, pytest.raises(exc_type):
        ptc.read_spatial_region(region=region)


def point_cloud_read_spatial_region_transform_setup(uri, transform, input_axes, kwargs):
    schema = pa.schema([("x", pa.float64()), ("y", pa.float64())])

    with soma.PointCloudDataFrame.create(uri, schema=schema, domain=[[-10000, 10000], [-10000, 10000], [0, 10]]) as ptc:
        pydict = {
            "soma_joinid": [1, 2, 3, 4, 5],
            "x": [10, 20, 30, 40, 50],
            "y": [4.1, 5.2, 6.3, 7.4, 8.5],
        }
        rb = pa.Table.from_pydict(pydict, schema=schema.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    output_names = ("x", "y")
    input_names = tuple(axis.name for axis in input_axes)

    with soma.PointCloudDataFrame.open(uri, "r") as ptc:
        read_spatial_region = ptc.read_spatial_region(
            region_transform=transform(input_axes=input_names, output_axes=output_names, **kwargs),
        )

        assert read_spatial_region.data.concat() == ptc.read().concat()

        data_coordinate_space = read_spatial_region.data_coordinate_space
        assert data_coordinate_space.axis_names == ptc.coordinate_space.axis_names == output_names
        assert len(data_coordinate_space) == len(output_names)

        output_coordinate_space = read_spatial_region.output_coordinate_space
        assert output_coordinate_space.axes == input_axes
        assert output_coordinate_space.axis_names == input_names
        assert len(output_coordinate_space) == len(input_names)

        assert isinstance(read_spatial_region.coordinate_transform, transform)

    return read_spatial_region.coordinate_transform


def test_point_cloud_read_spatial_region_identity_transform(tmp_path):
    coordinate_transform = point_cloud_read_spatial_region_transform_setup(
        uri=tmp_path.as_uri(),
        transform=soma.IdentityTransform,
        input_axes=(soma.Axis(name="x"), soma.Axis(name="y")),
        kwargs={},
    )

    # Identity will always have a scale of 1.0
    assert list(coordinate_transform.scale_factors) == [1, 1]
    assert coordinate_transform.scale == 1.0
    assert list(coordinate_transform.inverse_transform().scale_factors) == [1, 1]


@pytest.mark.parametrize(
    "name, scale, expected_scale_factors, expected_scale",
    [
        ("Identity", 1, [1, 1], 1.0),
        ("Scale up", 2, [0.5, 0.5], 0.5),
        ("Negative scale up", -2, [-0.5, -0.5], -0.5),
        ("Very large", 1e6, [1e-6, 1e-6], 1e-6),
        ("Very small", 1e-6, [1e6, 1e6], 1e6),
        ("Scale down", 0.5, [2, 2], 2.0),
        ("Negative scale down", -0.5, [-2, -2], -2.0),
        ("Inversion", -1, [-1, -1], -1.0),
        ("Fractional scale", 0.75, [4 / 3, 4 / 3], 4 / 3),
    ],
)
def test_point_cloud_read_spatial_region_uniform_scale_transform(
    tmp_path,
    name,
    scale,
    expected_scale_factors,
    expected_scale,
):
    coordinate_transform = point_cloud_read_spatial_region_transform_setup(
        uri=tmp_path.as_uri(),
        transform=soma.UniformScaleTransform,
        input_axes=(soma.Axis(name="x2"), soma.Axis(name="y2")),
        kwargs={"scale": scale},
    )

    assert list(coordinate_transform.scale_factors) == expected_scale_factors
    assert coordinate_transform.scale == expected_scale
    assert coordinate_transform.inverse_transform().scale == scale


def test_point_cloud_read_spatial_region_uniform_scale_transform_bad(tmp_path):
    # Error out for zero scale
    with pytest.raises(ZeroDivisionError):
        point_cloud_read_spatial_region_transform_setup(
            uri=tmp_path.as_uri(),
            transform=soma.UniformScaleTransform,
            input_axes=(soma.Axis(name="x"), soma.Axis(name="y")),
            kwargs={"scale": 0},
        )


@pytest.mark.parametrize(
    "name, scale_factors, expected_scale_factors",
    [
        ["Identity", (1, 1), (1, 1)],
        ["Scale up and down", (2, 0.5), (0.5, 2.0)],
        ["Very large and small scale", (1e6, 1e-6), (1e-6, 1e6)],
        ["Single inversion", (-1, 2), (-1, 0.5)],
        ["Double inversion", (-2, -0.5), (-0.5, -2)],
        ["One dimension inverted", (1, -1), (1, -1)],
        ["Scale up in both dimensions", (2, 3), (0.5, (1 / 3))],
        ["Small scale factors", (0.1, 0.1), (10, 10)],
    ],
)
def test_point_cloud_read_spatial_region_scale_transform(tmp_path, name, scale_factors, expected_scale_factors):
    coordinate_transform = point_cloud_read_spatial_region_transform_setup(
        uri=tmp_path.as_uri(),
        transform=soma.ScaleTransform,
        input_axes=(soma.Axis(name="x2"), soma.Axis(name="y2")),
        kwargs={"scale_factors": scale_factors},
    )

    assert tuple(coordinate_transform.scale_factors) == expected_scale_factors
    assert tuple(coordinate_transform.inverse_transform().scale_factors) == scale_factors


@pytest.mark.parametrize(
    "name, matrix, expected_aug_matrix",
    [
        (
            "Identity",
            np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        ),
        (
            "Translation",
            np.array([[1, 0, 2], [0, 1, 3], [0, 0, 1]]),
            np.array([[1, 0, -2], [0, 1, -3], [0, 0, 1]]),
        ),
        (
            "Scaling",
            np.array([[2, 0, 0], [0, 3, 0], [0, 0, 1]]),
            np.array([[0.5, 0, 0], [0, 1 / 3, 0], [0, 0, 1]]),
        ),
        (
            "Rotation (90 degrees counter-clockwise)",
            np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
            np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
        ),
        (
            "Reflection over the x-axis",
            np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]]),
            np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]]),
        ),
        (
            "Reflection over the y-axis",
            np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        ),
        (
            "Shearing (x direction)",
            np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]]),
            np.array([[1, -1, 0], [0, 1, 0], [0, 0, 1]]),
        ),
        (
            "Shearing (y direction)",
            np.array([[1, 0, 0], [1, 1, 0], [0, 0, 1]]),
            np.array([[1, 0, 0], [-1, 1, 0], [0, 0, 1]]),
        ),
        (
            "Scaling and translation",
            np.array([[2, 0, 1], [0, 3, 2], [0, 0, 1]]),
            np.array([[0.5, 0, -0.5], [0, 1 / 3, -2 / 3], [0, 0, 1]]),
        ),
        (
            "Rotation and scaling",
            np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]]),
            np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]]),
        ),
    ],
)
def test_point_cloud_read_spatial_region_affine_transform(tmp_path, name, matrix, expected_aug_matrix):
    coordinate_transform = point_cloud_read_spatial_region_transform_setup(
        uri=tmp_path.as_uri(),
        transform=soma.AffineTransform,
        input_axes=(soma.Axis(name="x2"), soma.Axis(name="y2")),
        kwargs={"matrix": matrix},
    )

    aug_matrix = coordinate_transform.augmented_matrix
    assert np.array_equal(aug_matrix, expected_aug_matrix)

    aug_matrix_inv = coordinate_transform.inverse_transform().augmented_matrix
    assert np.array_equal(aug_matrix_inv, matrix)


def test_point_cloud_read_spatial_region_region_coord_space(tmp_path):
    uri = tmp_path.as_uri()

    schema = pa.schema([("x", pa.float64()), ("y", pa.float64())])

    with soma.PointCloudDataFrame.create(uri, schema=schema, domain=[[-10000, 10000], [-10000, 10000], [0, 10]]) as ptc:
        pydict = {
            "soma_joinid": [1, 2, 3, 4, 5],
            "x": [10, 20, 30, 40, 50],
            "y": [4.1, 5.2, 6.3, 7.4, 8.5],
        }
        rb = pa.Table.from_pydict(pydict, schema=schema.insert(0, pa.field("soma_joinid", pa.int64())))
        ptc.write(rb)

    with soma.PointCloudDataFrame.open(uri, "r") as ptc:
        output = ptc.read_spatial_region()
        assert output.output_coordinate_space.axis_names == ("x", "y")

        output = ptc.read_spatial_region(
            region_transform=soma.IdentityTransform(input_axes=("a", "b"), output_axes=("x", "y")),
            region_coord_space=soma.CoordinateSpace([soma.Axis(name="a"), soma.Axis(name="b")]),
        )
        assert output.output_coordinate_space.axis_names == ("a", "b")

        # Cannot specify the output coordinate space when transform is None
        with pytest.raises(ValueError):
            ptc.read_spatial_region(region_coord_space=soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")]))

        # The input axes of the transform must match the axes of the coordinate
        # space the requested region is defined in
        with pytest.raises(ValueError):
            ptc.read_spatial_region(
                region_transform=soma.IdentityTransform(input_axes=("a", "b"), output_axes=("x", "y")),
                region_coord_space=soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")]),
            )


@pytest.mark.parametrize(
    "dtype",
    ["int8", "int16", "int32", "uint8", "uint16", "uint32", "uint64", "float32", "float64"],
)
def test_fragments_in_writes(tmp_path, dtype):
    uri = tmp_path.as_posix()

    # --- three dataframes, all with identical schema
    df_0 = pd.DataFrame(
        {
            "x": pd.Series([0, 1, 2, 3], dtype=dtype),
            "y": pd.Series([0, 1, 2, 3], dtype=dtype),
            "soma_joinid": pd.Series([0, 1, 2, 3], dtype=np.int64),
        },
    )
    df_1 = pd.DataFrame(
        {
            "x": pd.Series([4, 5, 6, 7], dtype=dtype),
            "y": pd.Series([0, 1, 2, 3], dtype=dtype),
            "soma_joinid": pd.Series([4, 5, 6, 7], dtype=np.int64),
        },
    )
    df_2 = pd.DataFrame(
        {
            "x": pd.Series([8, 9, 10, 11], dtype=dtype),
            "y": pd.Series([0, 1, 2, 3], dtype=dtype),
            "soma_joinid": pd.Series([8, 9, 10, 11], dtype=np.int64),
        },
    )
    expected_df = pd.concat([df_0, df_1, df_2], ignore_index=True)

    soma.PointCloudDataFrame.create(
        uri,
        schema=pa.Schema.from_pandas(df_0, preserve_index=False),
        domain=[[0, 11], [0, 11], [0, 11]],
    ).close()

    with soma.PointCloudDataFrame.open(uri, mode="w") as A:
        # Three-chunk table
        A.write(
            pa.concat_tables(
                [
                    pa.Table.from_pandas(df_0, preserve_index=False),
                    pa.Table.from_pandas(df_1, preserve_index=False),
                    pa.Table.from_pandas(df_2, preserve_index=False),
                ],
            ),
            platform_config=soma.TileDBWriteOptions(sort_coords=False),
        )

    # There should be a single fragment even though there are three chunks (and
    # therefore three submits) in the array because we only finalize once at
    # the end
    assert len(list((Path(uri) / "__commits").iterdir())) == 1
    assert len(list((Path(uri) / "__fragments").iterdir())) == 1

    with soma.open(uri) as A:
        df = A.read().concat().to_pandas()
        np.testing.assert_array_equal(df, expected_df)


@pytest.mark.parametrize(
    "delete_coords,value_filter,expected_index",
    [
        pytest.param((slice(-10, 10),), None, [], id="delete all with slice on x"),
        pytest.param(
            (slice(None, None), slice(None, None), slice(0, 24)), None, [], id="delete all with slice on soma_joinid"
        ),
        pytest.param((slice(None, None), (0, -2, -1, 2)), None, [3, 8, 13, 18, 23], id="delete by points on y"),
        pytest.param(
            (slice(None, None), pa.chunked_array([[0, -2], [-1, 2]])),
            None,
            [3, 8, 13, 18, 23],
            id="delete by points on y with chunked array",
        ),
        pytest.param((slice(None, 1),), None, np.arange(20, 25), id="delete by slice on x"),
        pytest.param((slice(None, None),), "soma_joinid < 11", np.arange(11, 25), id="delete by value_filter only"),
    ],
)
def test_delete_cells(tmp_path, delete_coords, value_filter, expected_index):
    uri = str(tmp_path)
    join_data = np.arange(25, dtype=np.int64)
    x_data, y_data = np.meshgrid(np.arange(-2, 3, dtype=np.int64), np.arange(-2, 3, dtype=np.int64), indexing="ij")
    x_data = x_data.flatten()
    y_data = y_data.flatten()
    with soma.PointCloudDataFrame.create(
        uri, schema=pa.schema([("x", pa.int64()), ("y", pa.int64())]), domain=[[-2, 2], [-2, 2], [0, 24]]
    ) as point_cloud:
        data = pa.Table.from_pydict({"soma_joinid": join_data, "x": x_data, "y": y_data})
        point_cloud.write(data)

    with soma.PointCloudDataFrame.open(uri, mode="d") as points:
        points.delete_cells(delete_coords, value_filter=value_filter)

    with soma.PointCloudDataFrame.open(uri) as points:
        actual_table = points.read(result_order="row-major").concat()

    expected_table = pa.Table.from_pydict(
        {
            "x": pa.array([x_data[index] for index in expected_index], type=pa.int64()),
            "y": pa.array([y_data[index] for index in expected_index], type=pa.int64()),
            "soma_joinid": pa.array([join_data[index] for index in expected_index], type=pa.int64()),
        },
        schema=pa.schema([
            pa.field("x", pa.int64(), nullable=False),
            pa.field("y", pa.int64(), nullable=False),
            pa.field("soma_joinid", pa.int64(), nullable=False),
        ]),
    )
    assert actual_table == expected_table


@pytest.mark.parametrize(
    "index_type,delete_coords,value_filter,expected_index",
    [
        pytest.param(pa.float32(), (slice(-1.0, 0.0),), None, np.arange(15, 25), id="delete slice float"),
        pytest.param(pa.float64(), (slice(-1, 0),), None, np.arange(15, 25), id="delete slice double"),
        pytest.param(
            pa.float32(),
            ((-1, 1, 0),),
            None,
            np.hstack((np.arange(5, 10), np.arange(15, 20))),
            id="delete points float",
        ),
        pytest.param(
            pa.float64(),
            ((-1, 1, 0),),
            None,
            np.hstack((np.arange(5, 10), np.arange(15, 20))),
            id="delete points double",
        ),
    ],
)
def test_delete_cells_floating_point(tmp_path, index_type, delete_coords, value_filter, expected_index):
    uri = str(tmp_path)
    join_data = np.arange(25, dtype=np.int64)
    x_data, y_data = np.meshgrid(np.linspace(-1, 1, num=5), np.linspace(-1, 1, num=5), indexing="ij")
    x_data = x_data.flatten()
    y_data = y_data.flatten()
    with soma.PointCloudDataFrame.create(
        uri, schema=pa.schema([("x", index_type), ("y", index_type)]), domain=[[-1, 1], [-1, 1], [0, 24]]
    ) as point_cloud:
        data = pa.Table.from_pydict({"soma_joinid": join_data, "x": x_data, "y": y_data}, point_cloud.schema)
        point_cloud.write(data)

    with soma.PointCloudDataFrame.open(uri, mode="d") as points:
        points.delete_cells(delete_coords, value_filter=value_filter)

    with soma.PointCloudDataFrame.open(uri) as points:
        actual_table = points.read(result_order="row-major").concat()

    expected_table = pa.Table.from_pydict(
        {
            "x": pa.array([x_data[index] for index in expected_index], type=index_type),
            "y": pa.array([y_data[index] for index in expected_index], type=index_type),
            "soma_joinid": pa.array([join_data[index] for index in expected_index], type=pa.int64()),
        },
        schema=pa.schema([
            pa.field("x", index_type, nullable=False),
            pa.field("y", index_type, nullable=False),
            pa.field("soma_joinid", pa.int64(), nullable=False),
        ]),
    )
    assert actual_table == expected_table


def test_delete_cells_exceptions(tmp_path):
    with soma.PointCloudDataFrame.create(
        str(tmp_path),
        coordinate_space=("x", "y", "z"),
        schema=pa.schema({"x": pa.float64(), "y": pa.float64(), "z": pa.float64()}),
        domain=[[-100, 100], [-100, 100], [-100, 100], [0, 1000]],
    ) as points:
        points.close()
    with soma.PointCloudDataFrame.open(str(tmp_path), mode="w") as points:
        assert points.mode == "w"
        with pytest.raises((soma.SOMAError, RuntimeError)):
            points.delete_cells((slice(1, 4),))

    with soma.PointCloudDataFrame.open(str(tmp_path), mode="r") as points:
        assert points.mode == "r"
        with pytest.raises((soma.SOMAError, RuntimeError)):
            points.delete_cells((slice(1, 4),))

    with soma.PointCloudDataFrame.open(str(tmp_path), mode="d") as points:
        with pytest.raises(IndexError):
            points.delete_cells((slice(-1000, -900),))
        with pytest.raises(IndexError):
            points.delete_cells((slice(None, None), (10_000, 80, 900)))
        with pytest.raises(ValueError):
            points.delete_cells(tuple())
        with pytest.raises(ValueError):
            points.delete_cells((slice(3, 1),))
        with pytest.raises(TypeError):
            points.delete_cells((("one", "five"),))


def test_point_cloud_dataframe_batch_size_not_implemented(tmp_path):
    uri = (tmp_path / "pc_df").as_posix()
    asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])

    with soma.PointCloudDataFrame.create(
        uri,
        schema=asch,
        domain=[[-10000, 10000], [-10000, 10000], [0, 10]],
    ) as pc_df:
        pydict = {}
        pydict["soma_joinid"] = [1, 2, 3, 4, 5]
        pydict["x"] = [10, 20, 30, 40, 50]
        pydict["y"] = [4.1, 5.2, 6.3, 7.4, 8.5]

        rb = pa.Table.from_pydict(pydict)
        pc_df.write(rb)

    with soma.PointCloudDataFrame.open(uri) as pc_df:
        result = pc_df.read().concat()
        assert result.num_rows == 5

        with pytest.raises(NotImplementedError, match=r"batch_size.*not yet implemented"):
            list(pc_df.read(batch_size=BatchSize(count=10)))
