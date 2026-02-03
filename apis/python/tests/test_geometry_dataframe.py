import numpy as np
import pyarrow as pa
import pytest
import shapely

import tiledbsoma as soma
from tiledbsoma._core_options import BatchSize


@pytest.fixture
def geom_dataframe(tmp_path):
    """Factory fixture for creating GeometryDataFrame."""

    def _make(schema=None, domain=None, data=None, subdir="geom_df"):
        if schema is None:
            schema = pa.schema([("quality", pa.float32())])
        if domain is None:
            domain = [[(-10, 10), (-10, 10)], [0, 100]]

        uri = (tmp_path / subdir).as_uri()

        with soma.GeometryDataFrame.create(uri, schema=schema, domain=domain) as geom:
            if data is not None:
                rb = pa.Table.from_pydict(data)
                geom.from_outlines(rb)

        return uri

    return _make


def test_geometry_domain_deprecated(tmp_path):
    uri = tmp_path.as_uri()

    asch = pa.schema([("quality", pa.float32())])

    with pytest.deprecated_call(), soma.GeometryDataFrame.create(uri, schema=asch, domain=None) as geom:
        soma_domain = geom.domain

        assert len(soma_domain) == 2

        for idx, name in enumerate(geom.index_column_names):
            if name == "soma_geometry":
                f64info = np.finfo(np.float64)

                assert len(soma_domain[idx][0]) == 2
                assert len(soma_domain[idx][1]) == 2

                for axis_idx, axis_name in enumerate(geom.coordinate_space.axis_names):
                    assert soma_domain[idx][0][axis_name] == f64info.min

                for axis_idx, axis_name in enumerate(geom.coordinate_space.axis_names):
                    assert soma_domain[idx][1][axis_name] == f64info.max


def test_geometry_domain(tmp_path):
    uri = tmp_path.as_uri()
    domain = [[(-1000, 1000), (0, 100)], [0, 100]]

    asch = pa.schema([("quality", pa.float32())])

    with soma.GeometryDataFrame.create(uri, schema=asch, domain=domain) as geom:
        soma_domain = geom.domain

        assert len(soma_domain) == 2

        for idx, name in enumerate(geom.index_column_names):
            if name == "soma_geometry":
                assert len(soma_domain[idx][0]) == 2
                assert len(soma_domain[idx][1]) == 2

                for axis_idx, axis_name in enumerate(geom.coordinate_space.axis_names):
                    assert soma_domain[idx][0][axis_name] == domain[0][axis_idx][0]

                for axis_idx, axis_name in enumerate(geom.coordinate_space.axis_names):
                    assert soma_domain[idx][1][axis_name] == domain[0][axis_idx][1]


def test_geometry_coordinate_space(tmp_path):
    uri = tmp_path.as_uri()

    asch = pa.schema([("quality", pa.float32())])

    with soma.GeometryDataFrame.create(uri, schema=asch, domain=[[(-100, 100), (-100, 100)], [0, 100]]) as geom:
        assert len(geom.coordinate_space) == 2
        assert geom.coordinate_space.axis_names == ("x", "y")
        assert geom.coordinate_space.axes == (soma.Axis(name="x"), soma.Axis(name="y"))

        # Axis names do not match
        with pytest.raises(ValueError):
            geom.coordinate_space = soma.CoordinateSpace([soma.Axis(name="a"), soma.Axis(name="y")])

        geom.coordinate_space = soma.CoordinateSpace([soma.Axis(name="x", unit="m"), soma.Axis(name="y", unit="in")])
        assert geom.coordinate_space[0] == soma.Axis(name="x", unit="m")
        assert geom.coordinate_space[1] == soma.Axis(name="y", unit="in")


def test_geometry_basic_read(geom_dataframe):
    triangle = shapely.Polygon([(0, 0), (0, 1), (1, 0), (0, 0)])
    rect = shapely.Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])

    data = {
        "soma_geometry": [
            [0.0, 0, 0, 1, 1, 0, 0, 0],
            [0.0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
        ],
        "soma_joinid": [1, 2],
        "quality": [4.1, 5.2],
    }
    uri = geom_dataframe(data=data)

    with soma.GeometryDataFrame.open(uri) as geom:
        result = geom.read().concat()

        assert result[0].to_numpy()[0] == triangle.wkb
        assert result[0].to_numpy()[1] == rect.wkb

        assert shapely.from_wkb(result["soma_geometry"].to_numpy()[0]) == shapely.Polygon(
            [(0, 0), (0, 1), (1, 0), (0, 0)],
        )
        assert shapely.from_wkb(result["soma_geometry"].to_numpy()[1]) == shapely.Polygon(
            [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)],
        )


def test_geometry_basic_spatial_read(geom_dataframe):
    data = {
        "soma_geometry": [
            [0.0, 0, 0, 1, 1, 0, 0, 0],
            [2.0, 0, 2, 1, 3, 1, 3, 0, 2, 0],
        ],
        "soma_joinid": [1, 2],
        "quality": [4.1, 5.2],
    }
    uri = geom_dataframe(data=data)

    with soma.GeometryDataFrame.open(uri) as geom:
        result = geom.read_spatial_region(region=[0.5, 0.5, 1.5, 1.5]).data.concat()

        # Internal columns will be hidden in a subsequent PR
        assert len(result) == 1

        assert shapely.from_wkb(result["soma_geometry"].to_numpy()[0]) == shapely.Polygon(
            [(0, 0), (0, 1), (1, 0), (0, 0)],
        )


def test_delete_cells_not_implemented(tmp_path):
    uri = tmp_path.as_uri()
    with soma.GeometryDataFrame.create(
        uri, schema=pa.schema([("quality", pa.float32())]), domain=[[(0, 100), (0, 100)], [0, 10]]
    ) as geom:
        geom.close()

    with soma.GeometryDataFrame.open(uri, mode="d") as geom:
        assert geom.mode == "d"
        with pytest.raises(NotImplementedError):
            geom.delete_cells((slice(None, None),))


def test_geometry_dataframe_batch_size_not_implemented(geom_dataframe):
    data = {
        "soma_geometry": [
            [0.0, 0, 0, 1, 1, 0, 0, 0],
            [0.0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
        ],
        "soma_joinid": [1, 2],
        "quality": [4.1, 5.2],
    }

    uri = geom_dataframe(data=data)

    with soma.GeometryDataFrame.open(uri) as geom_df:
        result = geom_df.read().concat()
        assert result.num_rows == 2

        with pytest.raises(NotImplementedError, match=r"batch_size.*not yet implemented"):
            list(geom_df.read(batch_size=BatchSize(count=10)))
