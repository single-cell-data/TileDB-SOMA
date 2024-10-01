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


def test_scene_point_cloud_with_affine_transform(tmp_path):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/", "test_scene_point_cloud_with_affine_transform"
    )

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        ptc_uri = urljoin(obsl_uri, "ptc")
        asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])
        coord_space = soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")])

        # TODO replace with Scene.add_new_point_cloud_dataframe when implemented
        scene["obsl"]["ptc"] = soma.PointCloudDataFrame.create(ptc_uri, schema=asch)

        transform = soma.AffineTransform(
            input_axes=("x", "y"),
            output_axes=("x", "y"),
            matrix=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_point_cloud_dataframe("transcripts", transform)

        scene.coordinate_space = coord_space

        # No PointCloudDataFrame named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_point_cloud_dataframe("bad", transform)

        scene.set_transform_to_point_cloud_dataframe("ptc", transform)
        assert np.array_equal(
            scene.get_transform_to_point_cloud_dataframe("ptc").augmented_matrix,
            transform.augmented_matrix,
        )


def test_scene_point_cloud_with_scale_transform(tmp_path):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/", "test_scene_point_cloud_with_scale_transform"
    )

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        ptc_uri = urljoin(obsl_uri, "ptc")
        asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])
        coord_space = soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")])

        # TODO replace with Scene.add_new_point_cloud_dataframe when implemented
        scene["obsl"]["ptc"] = soma.PointCloudDataFrame.create(ptc_uri, schema=asch)

        transform = soma.ScaleTransform(
            input_axes=("x", "y"), output_axes=("x", "y"), scale_factors=[1, 1]
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_point_cloud_dataframe("transcripts", transform)

        scene.coordinate_space = coord_space

        # No PointCloudDataFrame named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_point_cloud_dataframe("bad", transform)

        scene.set_transform_to_point_cloud_dataframe("ptc", transform)
        assert np.array_equal(
            scene.get_transform_to_point_cloud_dataframe("ptc").scale_factors,
            transform.scale_factors,
        )


def test_scene_point_cloud_with_uniform_scale_transform(tmp_path):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/", "test_scene_point_cloud_with_uniform_scale_transform"
    )

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        ptc_uri = urljoin(obsl_uri, "ptc")
        asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])
        coord_space = soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")])

        # TODO replace with Scene.add_new_point_cloud_dataframe when implemented
        scene["obsl"]["ptc"] = soma.PointCloudDataFrame.create(ptc_uri, schema=asch)

        transform = soma.UniformScaleTransform(
            input_axes=("x", "y"), output_axes=("x", "y"), scale=1
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_point_cloud_dataframe("transcripts", transform)

        scene.coordinate_space = coord_space

        # No PointCloudDataFrame named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_point_cloud_dataframe("bad", transform)

        scene.set_transform_to_point_cloud_dataframe("ptc", transform)
        assert (
            scene.get_transform_to_point_cloud_dataframe("ptc").scale == transform.scale
        )


def test_scene_point_cloud_with_identity_scale_transform(tmp_path):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/", "test_scene_point_cloud_with_identity_scale_transform"
    )

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        ptc_uri = urljoin(obsl_uri, "ptc")
        asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])
        coord_space = soma.CoordinateSpace([soma.Axis(name="x"), soma.Axis(name="y")])

        # TODO replace with Scene.add_new_point_cloud_dataframe when implemented
        scene["obsl"]["ptc"] = soma.PointCloudDataFrame.create(ptc_uri, schema=asch)

        transform = soma.IdentityTransform(
            input_axes=("x", "y"), output_axes=("x", "y")
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_point_cloud_dataframe("transcripts", transform)

        scene.coordinate_space = coord_space

        # No PointCloudDataFrame named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_point_cloud_dataframe("bad", transform)

        scene.set_transform_to_point_cloud_dataframe("ptc", transform)
        assert (
            scene.get_transform_to_point_cloud_dataframe("ptc").scale == transform.scale
        )


def test_scene_multiscale_image_with_affine_transform(tmp_path):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/", "test_scene_multiscale_image_with_affine_transform"
    )

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

        transform = soma.AffineTransform(
            input_axes=("x", "y"),
            output_axes=("x", "y"),
            matrix=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_multiscale_image("transcripts", transform)

        scene.coordinate_space = coord_space

        # No PointCloudDataFrame named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_multiscale_image("bad", transform)

        scene.set_transform_to_multiscale_image("msi", transform)
        assert np.array_equal(
            scene.get_transform_to_multiscale_image("msi").augmented_matrix,
            transform.augmented_matrix,
        )


def test_scene_multiscale_image_with_scale_transform(tmp_path):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/", "test_scene_multiscale_image_with_scale_transform"
    )

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

        transform = soma.ScaleTransform(
            input_axes=("x", "y"), output_axes=("x", "y"), scale_factors=[1, 1]
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_multiscale_image("transcripts", transform)

        scene.coordinate_space = coord_space

        # No PointCloudDataFrame named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_multiscale_image("bad", transform)

        scene.set_transform_to_multiscale_image("msi", transform)
        assert np.array_equal(
            scene.get_transform_to_multiscale_image("msi").scale_factors,
            transform.scale_factors,
        )


def test_scene_multiscale_image_with_uniform_scale_transform(tmp_path):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/",
        "test_scene_multiscale_image_with_uniform_scale_transform",
    )

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

        transform = soma.UniformScaleTransform(
            input_axes=("x", "y"), output_axes=("x", "y"), scale=1
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_multiscale_image("transcripts", transform)

        scene.coordinate_space = coord_space

        # No PointCloudDataFrame named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_multiscale_image("bad", transform)

        scene.set_transform_to_multiscale_image("msi", transform)
        assert scene.get_transform_to_multiscale_image("msi").scale == transform.scale


def test_scene_multiscale_image_with_identity_scale_transform(tmp_path):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/",
        "test_scene_multiscale_image_with_identity_scale_transform",
    )

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

        transform = soma.IdentityTransform(
            input_axes=("x", "y"), output_axes=("x", "y")
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_multiscale_image("transcripts", transform)

        scene.coordinate_space = coord_space

        # No PointCloudDataFrame named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_multiscale_image("bad", transform)

        scene.set_transform_to_multiscale_image("msi", transform)
        assert scene.get_transform_to_multiscale_image("msi").scale == transform.scale
