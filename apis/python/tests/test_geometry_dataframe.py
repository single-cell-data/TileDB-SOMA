import numpy as np
import pyarrow as pa
import pytest

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
