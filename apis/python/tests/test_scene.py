import json
from urllib.parse import urljoin

import numpy as np
import pyarrow as pa
import pytest
import typeguard

import tiledbsoma as soma


def create_and_populate_df(uri: str) -> soma.DataFrame:
    obs_arrow_schema = pa.schema(
        [
            ("foo", pa.int32()),
            ("bar", pa.float64()),
            ("baz", pa.large_string()),
        ]
    )

    with soma.DataFrame.create(uri, schema=obs_arrow_schema) as obs:
        pydict = {}
        pydict["soma_joinid"] = [0, 1, 2, 3, 4]
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        rb = pa.Table.from_pydict(pydict)
        obs.write(rb)

    return soma.DataFrame.open(uri)


def test_scene_basic(tmp_path):
    baseuri = tmp_path.as_uri()

    with soma.Scene.create(baseuri) as scene:
        assert scene.uri == baseuri

        with pytest.raises(TypeError):
            scene["obsl"] = soma.Experiment.create(urljoin(baseuri, "obs"))
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)
        scene["obsl"]["df"] = create_and_populate_df(urljoin(obsl_uri, "df"))

        with pytest.raises(TypeError):
            scene["varl"] = soma.Measurement.create(urljoin(baseuri, "var"))
        varl_uri = urljoin(baseuri, "varl")
        scene["varl"] = soma.Collection.create(varl_uri)
        scene["varl"]["sparse"] = soma.SparseNDArray.create(
            urljoin(varl_uri, "sparse"), type=pa.int64(), shape=(10,)
        )
        scene["varl"]["dense"] = soma.DenseNDArray.create(
            urljoin(varl_uri, "dense"), type=pa.int64(), shape=(10,)
        )

        img_uri = urljoin(baseuri, "img")
        scene["img"] = soma.Collection.create(img_uri)
        scene["img"]["col"] = soma.Collection.create(urljoin(img_uri, "col"))

    assert not soma.Collection.exists(baseuri)
    assert soma.Scene.exists(baseuri)
    assert soma.Collection.exists(obsl_uri)
    assert soma.Collection.exists(varl_uri)
    assert soma.Collection.exists(img_uri)
    assert soma.Measurement.exists(urljoin(baseuri, "var"))
    assert soma.SparseNDArray.exists(urljoin(varl_uri, "sparse"))
    assert soma.DenseNDArray.exists(urljoin(varl_uri, "dense"))
    assert soma.Collection.exists(urljoin(img_uri, "col"))

    with soma.Scene.open(baseuri) as scene:
        assert scene is not None
        assert scene.obsl is not None
        assert scene.obsl["df"] is not None
        assert scene.varl is not None
        assert scene.varl["sparse"] is not None
        assert scene.varl["dense"] is not None
        assert scene.img is not None
        assert scene.img["col"] is not None

        assert len(scene) == 3
        assert scene.soma_type == "SOMAScene"

        assert scene.obsl == scene["obsl"]
        assert len(scene.obsl) == 1
        assert scene.obsl["df"] == scene["obsl"]["df"]

        assert scene.varl == scene["varl"]
        assert len(scene.varl) == 2
        assert scene.varl["sparse"] == scene["varl"]["sparse"]
        assert scene.varl["dense"] == scene["varl"]["dense"]

        assert scene.img == scene["img"]
        assert len(scene.img) == 1
        assert scene.img["col"] == scene["img"]["col"]

    with pytest.raises(soma.DoesNotExistError):
        soma.Scene.open("bad uri")


def test_measurement_with_var_scene(tmp_path):
    baseuri = tmp_path.as_uri()
    obs_scene_uri = urljoin(baseuri, "obs_scene")

    with soma.Measurement.create(baseuri) as mea:
        with pytest.raises(TypeError):
            mea["obs_scene"] = soma.SparseNDArray.create(obs_scene_uri)
        mea["obs_scene"] = create_and_populate_df(obs_scene_uri)

    assert soma.Measurement.exists(baseuri)
    assert soma.DataFrame.exists(obs_scene_uri)


def test_scene_coord_space(tmp_path):
    uri = tmp_path.as_uri()

    coord_space = soma.CoordinateSpace(
        [
            soma.Axis(name="x"),
            soma.Axis(name="y"),
        ]
    )
    coord_space_json = """
    [
        {"name": "x", "unit": null},
        {"name": "y", "unit": null}
    ]
    """

    with soma.Scene.create(uri) as scene:
        assert scene.coordinate_space is None
        assert "soma_coordinate_space" not in scene.metadata

        # Setter only takes in CoordinateSpace
        with pytest.raises(typeguard.TypeCheckError):
            scene.coordinate_space = None
        with pytest.raises(typeguard.TypeCheckError):
            scene.coordinate_space = [soma.Axis(name="x"), soma.Axis(name="y")]

        # Reserved metadata key should not be settable?
        # with pytest.raises(soma.SOMAError):
        #     scene.metadata["soma_coordinate_space"] = coord_space_json

        scene.coordinate_space = coord_space
        assert scene.coordinate_space == coord_space
        assert json.loads(scene.metadata["soma_coordinate_space"]) == json.loads(
            coord_space_json
        )

    with soma.Scene.open(uri) as scene:
        assert scene.coordinate_space == coord_space


@pytest.mark.parametrize(
    "coord_transform, transform_kwargs",
    [
        (soma.AffineTransform, {"matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}),
        (soma.ScaleTransform, {"scale_factors": [1, 1]}),
        (soma.UniformScaleTransform, {"scale": 1}),
        (soma.IdentityTransform, {}),
    ],
)
def test_scene_point_cloud(tmp_path, coord_transform, transform_kwargs):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_scene_point_cloud")

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        ptc_uri = urljoin(obsl_uri, "ptc")
        asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])
        coord_space = soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")])

        # TODO replace with Scene.add_new_point_cloud_dataframe when implemented
        scene["obsl"]["ptc"] = soma.PointCloudDataFrame.create(ptc_uri, schema=asch)

        transform = coord_transform(
            input_axes=("x", "y"), output_axes=("x", "y"), **transform_kwargs
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_point_cloud_dataframe("ptc", transform)

        scene.coordinate_space = coord_space

        # No SOMAObject named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_point_cloud_dataframe("bad", transform)

        # Not a PointCloudDataFrame
        scene["obsl"]["col"] = soma.Collection.create(urljoin(obsl_uri, "col"))
        with pytest.raises(typeguard.TypeCheckError):
            scene.set_transform_to_point_cloud_dataframe("col", transform)

        # Transform not set
        with pytest.raises(KeyError):
            scene.get_transform_to_point_cloud_dataframe("ptc")

        scene.set_transform_to_point_cloud_dataframe("ptc", transform)

        ptc_transform = scene.get_transform_to_point_cloud_dataframe("ptc")
        if isinstance(coord_transform, soma.AffineTransform):
            assert np.array_equal(
                ptc_transform.augmented_matrix,
                transform.augmented_matrix,
            )
        elif isinstance(coord_transform, soma.ScaleTransform):
            assert np.array_equal(
                ptc_transform.scale_factors,
                transform.scale_factors,
            )
        elif isinstance(
            coord_transform, (soma.UniformScaleTransform, soma.IdentityTransform)
        ):
            assert ptc_transform.scale == transform.scale


@pytest.mark.parametrize(
    "coord_transform, transform_kwargs",
    [
        (soma.AffineTransform, {"matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}),
        (soma.ScaleTransform, {"scale_factors": [1, 1]}),
        (soma.UniformScaleTransform, {"scale": 1}),
        (soma.IdentityTransform, {}),
    ],
)
def test_scene_multiscale_image(tmp_path, coord_transform, transform_kwargs):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_scene_multiscale_image")

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        img_uri = urljoin(baseuri, "img")
        scene["img"] = soma.Collection.create(img_uri)

        msi_uri = urljoin(img_uri, "msi")
        coord_space = soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")])

        # TODO replace with Scene.add_multiscale_image when implemented
        scene["img"]["msi"] = soma.MultiscaleImage.create(
            msi_uri, type=pa.int64(), reference_level_shape=[1, 2, 3]
        )

        transform = coord_transform(
            input_axes=("x", "y"),
            output_axes=("x", "y"),
            **transform_kwargs,
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_multiscale_image("msi", transform)

        scene.coordinate_space = coord_space

        # No MultiscaleImage named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_multiscale_image("bad", transform)

        # Transform not set
        with pytest.raises(KeyError):
            scene.get_transform_to_multiscale_image("msi")

        # Not a MultiscaleImage
        scene["img"]["col"] = soma.Collection.create(urljoin(img_uri, "col"))
        with pytest.raises(typeguard.TypeCheckError):
            scene.set_transform_to_multiscale_image("col", transform)

        scene.set_transform_to_multiscale_image("msi", transform)

        msi_transform = scene.get_transform_to_multiscale_image("msi")
        if isinstance(coord_transform, soma.AffineTransform):
            assert np.array_equal(
                msi_transform.augmented_matrix,
                transform.augmented_matrix,
            )
        elif isinstance(coord_transform, soma.ScaleTransform):
            assert np.array_equal(
                msi_transform.scale_factors,
                transform.scale_factors,
            )
        elif isinstance(
            coord_transform, (soma.UniformScaleTransform, soma.IdentityTransform)
        ):
            assert msi_transform.scale == transform.scale


@pytest.mark.skip("GeometryDataFrame not supported yet")
@pytest.mark.parametrize(
    "coord_transform, transform_kwargs",
    [
        (soma.AffineTransform, {"matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}),
        (soma.ScaleTransform, {"scale_factors": [1, 1]}),
        (soma.UniformScaleTransform, {"scale": 1}),
        (soma.IdentityTransform, {}),
    ],
)
def test_scene_geometry_dataframe(tmp_path, coord_transform, transform_kwargs):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_scene_geometry_dataframe")

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        gdf_uri = urljoin(obsl_uri, "gdf")
        asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])
        coord_space = soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")])

        # TODO replace with Scene.add_new_geometry_dataframe when implemented
        scene["obsl"]["gdf"] = soma.GeometryDataFrame.create(gdf_uri, schema=asch)

        transform = coord_transform(
            input_axes=("x", "y"), output_axes=("x", "y"), **transform_kwargs
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_geometry_dataframe("gdf", transform)

        scene.coordinate_space = coord_space

        # No SOMAObject named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_geometry_dataframe("bad", transform)

        # Not a GeometryDataFrame
        scene["obsl"]["col"] = soma.Collection.create(urljoin(obsl_uri, "col"))
        with pytest.raises(typeguard.TypeCheckError):
            scene.set_transform_to_geometry_dataframe("col", transform)

        # Transform not set
        with pytest.raises(KeyError):
            scene.get_transform_to_geometry_dataframe("gdf")

        scene.set_transform_to_geometry_dataframe("gdf", transform)

        gdf_transform = scene.get_transform_to_geometry_dataframe("gdf")
        if isinstance(coord_transform, soma.AffineTransform):
            assert np.array_equal(
                gdf_transform.augmented_matrix,
                transform.augmented_matrix,
            )
        elif isinstance(coord_transform, soma.ScaleTransform):
            assert np.array_equal(
                gdf_transform.scale_factors,
                transform.scale_factors,
            )
        elif isinstance(
            coord_transform, (soma.UniformScaleTransform, soma.IdentityTransform)
        ):
            assert gdf_transform.scale == transform.scale
