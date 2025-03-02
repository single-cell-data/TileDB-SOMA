import numpy as np
import pyarrow as pa
import pytest
import shapely

import tiledbsoma as soma


@pytest.mark.parametrize("domain", [None, [[(-1000, 1000), (0, 100)], [0, 100]]])
def test_geometry_domain(tmp_path, domain):
    uri = tmp_path.as_uri()

    asch = pa.schema([("quality", pa.float32())])

    with soma.GeometryDataFrame.create(uri, schema=asch, domain=domain) as geom:
        soma_domain = geom.domain

        assert len(soma_domain) == 2

        for idx, name in enumerate(geom.index_column_names):
            if name == "soma_geometry":
                f64info = np.finfo(np.float64)

                assert len(soma_domain[idx][0]) == 2
                assert len(soma_domain[idx][1]) == 2

                for axis_idx, axis_name in enumerate(geom.coordinate_space.axis_names):
                    assert soma_domain[idx][0][axis_name] == (
                        f64info.min if domain is None else domain[0][axis_idx][0]
                    )

                for axis_idx, axis_name in enumerate(geom.coordinate_space.axis_names):
                    assert soma_domain[idx][1][axis_name] == (
                        f64info.max if domain is None else domain[0][axis_idx][1]
                    )


def test_geometry_coordinate_space(tmp_path):
    uri = tmp_path.as_uri()

    asch = pa.schema([("quality", pa.float32())])

    with soma.GeometryDataFrame.create(uri, schema=asch) as geom:
        assert len(geom.coordinate_space) == 2
        assert geom.coordinate_space.axis_names == ("x", "y")
        assert geom.coordinate_space.axes == (soma.Axis(name="x"), soma.Axis(name="y"))

        # Axis names do not match
        with pytest.raises(ValueError):
            geom.coordinate_space = soma.CoordinateSpace(
                [soma.Axis(name="a"), soma.Axis(name="y")]
            )

        geom.coordinate_space = soma.CoordinateSpace(
            [soma.Axis(name="x", unit="m"), soma.Axis(name="y", unit="in")]
        )
        assert geom.coordinate_space[0] == soma.Axis(name="x", unit="m")
        assert geom.coordinate_space[1] == soma.Axis(name="y", unit="in")


def test_geometry_basic_read(tmp_path):
    uri = tmp_path.as_uri()

    asch = pa.schema([("quality", pa.float32())])
    triangle = shapely.Polygon([(0, 0), (0, 1), (1, 0), (0, 0)])
    rect = shapely.Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])

    with soma.GeometryDataFrame.create(
        uri, schema=asch, domain=[[(-10, 10), (-10, 10)], [0, 100]]
    ) as geom:
        pydict = {}
        pydict["soma_geometry"] = [
            [0.0, 0, 0, 1, 1, 0, 0, 0],
            [0.0, 0, 0, 1, 1, 1, 1, 0, 0, 0],
        ]
        pydict["soma_joinid"] = [1, 2]
        pydict["quality"] = [4.1, 5.2]

        rb = pa.Table.from_pydict(pydict)
        geom.from_outlines(rb)

    with soma.GeometryDataFrame.open(uri) as geom:
        result = geom.read().concat()

        assert result[0].to_numpy()[0] == triangle.wkb
        assert result[0].to_numpy()[1] == rect.wkb

        assert shapely.from_wkb(
            result["soma_geometry"].to_numpy()[0]
        ) == shapely.Polygon([(0, 0), (0, 1), (1, 0), (0, 0)])
        assert shapely.from_wkb(
            result["soma_geometry"].to_numpy()[1]
        ) == shapely.Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])


def test_geometry_basic_spatial_read(tmp_path):
    uri = tmp_path.as_uri()

    asch = pa.schema([("quality", pa.float32())])

    with soma.GeometryDataFrame.create(
        uri, schema=asch, domain=[[(-10, 10), (-10, 10)], [0, 100]]
    ) as geom:
        pydict = {}
        pydict["soma_geometry"] = [
            [0.0, 0, 0, 1, 1, 0, 0, 0],
            [2.0, 0, 2, 1, 3, 1, 3, 0, 2, 0],
        ]
        pydict["soma_joinid"] = [1, 2]
        pydict["quality"] = [4.1, 5.2]

        rb = pa.Table.from_pydict(pydict)
        geom.from_outlines(rb)

    with soma.GeometryDataFrame.open(uri) as geom:

        result = geom.read_spatial_region(region=[0.5, 0.5, 1.5, 1.5]).data.concat()

        # Internal columns will be hidden in a subsequent PR
        assert len(result) == 1

        assert shapely.from_wkb(
            result["soma_geometry"].to_numpy()[0]
        ) == shapely.Polygon([(0, 0), (0, 1), (1, 0), (0, 0)])
