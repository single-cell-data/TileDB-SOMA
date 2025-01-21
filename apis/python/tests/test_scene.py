import json
from urllib.parse import urljoin

import pyarrow as pa
import pytest
import typeguard

import tiledbsoma as soma

from ._util import assert_transform_equal


def test_scene_basic(tmp_path):
    uri = f"{tmp_path.as_uri()}/scene_basic"

    with soma.Scene.create(uri) as scene:
        assert scene.uri == uri

        # Test `obsl`, `varl`, and `img` throw type errors for invalid collection types.
        with pytest.raises(TypeError):
            scene["obsl"] = soma.Measurement.create(f"{uri}/obsl_")
        with pytest.raises(TypeError):
            scene["varl"] = soma.Measurement.create(f"{uri}/varl_")
        with pytest.raises(TypeError):
            scene["img"] = soma.Measurement.create(f"{uri}/img_")
        with pytest.raises(TypeError):
            scene.add_new_collection("obsl", soma.Experiment)
        with pytest.raises(TypeError):
            scene.add_new_collection("varl", soma.Experiment)
        with pytest.raises(TypeError):
            scene.add_new_collection("img", soma.Experiment)

        # Add obsl.
        with scene.add_new_collection("obsl") as obsl:
            obsl.add_new_dataframe(
                "df",
                schema=pa.schema(
                    [
                        ("foo", pa.int32()),
                        ("bar", pa.float64()),
                        ("baz", pa.large_string()),
                    ]
                ),
                domain=[[0, 9]],
            )

        # Add varl.
        with scene.add_new_collection("varl") as varl:
            varl.add_new_collection("col")

        # Add img
        with scene.add_new_collection("img") as img:
            img.add_new_sparse_ndarray("sparse", type=pa.int64(), shape=(10,))
            img.add_new_dense_ndarray("dense", type=pa.int64(), shape=(10,))

    assert not soma.Collection.exists(uri)
    assert soma.Scene.exists(uri)
    assert soma.Collection.exists(f"{uri}/obsl")
    assert soma.DataFrame.exists(f"{uri}/obsl/df")
    assert soma.Collection.exists(f"{uri}/varl")
    assert soma.Collection.exists(f"{uri}/varl/col")
    assert soma.Collection.exists(f"{uri}/img")
    assert soma.SparseNDArray.exists(f"{uri}/img/sparse")
    assert soma.DenseNDArray.exists(f"{uri}/img/dense")

    with soma.Scene.open(uri) as scene:
        assert scene is not None
        assert scene.obsl is not None
        assert scene.obsl["df"] is not None
        assert scene.varl is not None
        assert scene.varl["col"] is not None
        assert scene.img["sparse"] is not None
        assert scene.img["dense"] is not None

        assert len(scene) == 3
        assert scene.soma_type == "SOMAScene"

        assert (
            scene.metadata[soma._constants.SOMA_SPATIAL_VERSION_METADATA_KEY]
            == soma._constants.SOMA_SPATIAL_ENCODING_VERSION
        )

        assert scene.obsl == scene["obsl"]
        assert len(scene.obsl) == 1
        assert scene.obsl["df"] == scene["obsl"]["df"]

        assert scene.varl == scene["varl"]
        assert len(scene.varl) == 1
        assert scene.varl["col"] == scene["varl"]["col"]

        assert scene.img == scene["img"]
        assert len(scene.img) == 2
        assert scene.img["sparse"] == scene["img"]["sparse"]
        assert scene.img["dense"] == scene["img"]["dense"]

    with pytest.raises(soma.DoesNotExistError):
        soma.Scene.open("bad uri")

    # Ensure it cannot be opened by another type
    with pytest.raises(soma.SOMAError):
        soma.DataFrame.open(uri)

    with pytest.raises(soma.SOMAError):
        soma.SparseNDArray.open(uri)

    with pytest.raises(soma.SOMAError):
        soma.DenseNDArray.open(uri)

    with pytest.raises(soma.SOMAError):
        soma.PointCloudDataFrame.open(uri)

    with pytest.raises(soma.SOMAError):
        soma.Collection.open(uri)

    with pytest.raises(soma.SOMAError):
        soma.Experiment.open(uri)

    with pytest.raises(soma.SOMAError):
        soma.Measurement.open(uri)

    with pytest.raises(soma.SOMAError):
        soma.MultiscaleImage.open(uri)


@pytest.mark.parametrize(
    "input,expected",
    (
        (["x", "y"], soma.CoordinateSpace.from_axis_names(("x", "y"))),
        (
            soma.CoordinateSpace([soma.Axis("x", "meters"), soma.Axis("y", "meters")]),
            soma.CoordinateSpace([soma.Axis("x", "meters"), soma.Axis("y", "meters")]),
        ),
    ),
)
def test_scene_coord_space_at_create(tmp_path, input, expected):
    uri = tmp_path.as_uri()

    with soma.Scene.create(uri, coordinate_space=input) as scene:

        # Reserved metadata key should not be settable?
        # with pytest.raises(soma.SOMAError):
        #     scene.metadata["soma_coordinate_space"] = "user_metadata"

        assert scene.coordinate_space == expected

    with soma.Scene.open(uri) as scene:
        assert scene.coordinate_space == expected


def test_scene_coord_space_after_create(tmp_path):
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


class TestSceneDeepSubcollections:
    """Tests on a Scene with multiple layers of subcollections.

    Scene structure:

        scene
        |
        ├- obsl
        ├- varl
        |   └- RNA
        └- suns
            └- suns
                └- final
    """

    @pytest.fixture(scope="class")
    def scene(self, tmp_path_factory):
        """Creates and returns a scene for reading that"""
        baseuri = tmp_path_factory.mktemp("scene").as_uri()
        scene_uri = urljoin(baseuri, "multi-collection")

        # Create a scene with multi-level collections.
        with soma.Scene.create(scene_uri) as scene:

            obsl = scene.add_new_collection("obsl")
            obsl.metadata["name"] = "obsl"

            varl = scene.add_new_collection("varl")
            varl.metadata["name"] = "varl"

            rna = varl.add_new_collection("RNA")
            rna.metadata["name"] = "varl/RNA"

            # Add a collection that is not part of the set data model.
            # Using 'suns' for spatial-uns.
            suns = scene.add_new_collection("suns")
            suns.metadata["name"] = "suns"

            suns2 = suns.add_new_collection("suns")
            suns2.metadata["name"] = "suns/suns"

            fin = suns2.add_new_collection("final")
            fin.metadata["name"] = "suns/suns/final"

        scene = scene.open(scene_uri)
        yield scene
        scene.close()

    def test_open_subcollection_no_items(self, scene):
        with pytest.raises(ValueError):
            scene._open_subcollection([])

    @pytest.mark.parametrize(
        "subcollection",
        ["bad_name", ["obsl", "bad_name"], ["bad_name", "obsl"]],
    )
    def test_open_subcollection_keyerror(self, scene, subcollection):
        with pytest.raises(KeyError):
            scene._open_subcollection(subcollection)

    @pytest.mark.parametrize(
        "subcollection,expected_metadata",
        [
            ("obsl", "obsl"),
            (["obsl"], "obsl"),
            ("varl", "varl"),
            (["varl", "RNA"], "varl/RNA"),
            ("suns", "suns"),
            (["suns", "suns"], "suns/suns"),
            (["suns", "suns", "final"], "suns/suns/final"),
        ],
    )
    def test_open_subcolletion(self, scene, subcollection, expected_metadata):
        coll = scene._open_subcollection(subcollection)
        actual_metadata = coll.metadata["name"]
        assert actual_metadata == expected_metadata


def test_scene_point_cloud(tmp_path):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_scene_point_cloud")

    with soma.Scene.create(baseuri) as scene:
        # Create obsl.
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        # Add parameters for the point cloud.
        asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])
        elem_coord_space = soma.CoordinateSpace(
            [soma.Axis(name="x", unit="nm"), soma.Axis(name="y", unit="nm")]
        )
        transform = soma.ScaleTransform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("x", "y"),
            scale_factors=[-1, 1],
        )

        # Cannot set transform before the scene coordinate space is set.
        with pytest.raises(soma.SOMAError):
            scene.add_new_point_cloud_dataframe(
                "ptc",
                subcollection="obsl",
                transform=transform,
                schema=asch,
                coordinate_space=elem_coord_space,
            )

        # Set scene coordinate space.
        scene_coord_space = soma.CoordinateSpace(
            [soma.Axis(name="x_scene"), soma.Axis(name="y_scene")]
        )
        scene.coordinate_space = scene_coord_space

        # Mismatch in transform input axes and coordinate space axes.
        bad_transform = soma.ScaleTransform(
            input_axes=("xbad", "ybad"),
            output_axes=("x", "y"),
            scale_factors=[-1, 1],
        )
        with pytest.raises(ValueError):
            scene.add_new_point_cloud_dataframe(
                "ptc",
                subcollection="obsl",
                transform=bad_transform,
                schema=asch,
                coordinate_space=elem_coord_space,
            )

        # Mismatch in transform output axes and point cloud axes.
        bad_transform = soma.ScaleTransform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("xbad", "ybad"),
            scale_factors=[-1, 1],
        )
        with pytest.raises(ValueError):
            scene.add_new_point_cloud_dataframe(
                "ptc",
                subcollection="obsl",
                transform=bad_transform,
                schema=asch,
                coordinate_space=elem_coord_space,
            )

        # Add the point cloud dataframe.
        scene.add_new_point_cloud_dataframe(
            "ptc",
            subcollection="obsl",
            transform=transform,
            schema=asch,
            coordinate_space=elem_coord_space,
        )

        # Check the transform.
        ptc_transform = scene.get_transform_to_point_cloud_dataframe("ptc")
        assert_transform_equal(ptc_transform, transform)


@pytest.mark.parametrize(
    "coord_transform, transform_kwargs",
    [
        (soma.AffineTransform, {"matrix": [[1, 0, 1], [0, 1, 1], [0, 0, 1]]}),
        (soma.ScaleTransform, {"scale_factors": [-1, 1]}),
        (soma.UniformScaleTransform, {"scale": 2}),
        (soma.IdentityTransform, {}),
    ],
)
@pytest.mark.parametrize("set_coord_space", [True, False])
def test_scene_set_transform_to_point_cloud(
    tmp_path, coord_transform, transform_kwargs, set_coord_space
):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/", "test_scene_set_transform_to_point_cloud"
    )

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        asch = pa.schema([("x", pa.float64()), ("y", pa.float64())])
        coord_space = soma.CoordinateSpace(
            [soma.Axis(name="x_scene"), soma.Axis(name="y_scene")]
        )

        scene.add_new_point_cloud_dataframe(
            "ptc", subcollection="obsl", transform=None, schema=asch
        )

        transform = coord_transform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("x", "y"),
            **transform_kwargs,
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_point_cloud_dataframe("ptc", transform=transform)

        scene.coordinate_space = coord_space

        # No SOMAObject named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_point_cloud_dataframe("bad", transform=transform)

        # Mismatched input axes.
        transform_bad = coord_transform(
            input_axes=("x", "y"),
            output_axes=("x", "y"),
            **transform_kwargs,
        )
        with pytest.raises(ValueError):
            scene.set_transform_to_point_cloud_dataframe("ptc", transform=transform_bad)

        # Mismatched output axes.
        transform_bad = coord_transform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("x_scene", "y_scene"),
            **transform_kwargs,
        )
        with pytest.raises(ValueError):
            scene.set_transform_to_point_cloud_dataframe("ptc", transform=transform_bad)

        # Not a PointCloudDataFrame
        scene["obsl"]["col"] = soma.Collection.create(urljoin(obsl_uri, "col"))
        with pytest.raises(TypeError):
            scene.set_transform_to_point_cloud_dataframe("col", transform=transform)

        # Transform not set
        with pytest.raises(KeyError):
            scene.get_transform_to_point_cloud_dataframe("ptc")

        if set_coord_space:
            bad_coord_space = soma.CoordinateSpace.from_axis_names(("xbad", "ybad"))
            with pytest.raises(ValueError):
                scene.set_transform_to_point_cloud_dataframe(
                    "ptc", transform=transform, coordinate_space=bad_coord_space
                )

            coord_space = soma.CoordinateSpace(
                (soma.Axis(name="x", unit="nm"), soma.Axis(name="y", unit="nm"))
            )

            point_cloud = scene.set_transform_to_point_cloud_dataframe(
                "ptc", transform=transform, coordinate_space=coord_space
            )
            actual_coord_space = point_cloud.coordinate_space
            assert actual_coord_space == coord_space

        else:
            scene.set_transform_to_point_cloud_dataframe("ptc", transform=transform)

        ptc_transform = scene.get_transform_to_point_cloud_dataframe("ptc")
        assert_transform_equal(ptc_transform, transform)

        inv_transform = transform.inverse_transform()
        ptc_inv_transform = scene.get_transform_from_point_cloud_dataframe("ptc")
        assert_transform_equal(ptc_inv_transform, inv_transform)


def test_scene_multiscale_image(tmp_path):
    baseuri = urljoin(f"{tmp_path.as_uri()}/", "test_scene_multiscale_image")

    with soma.Scene.create(baseuri) as scene:
        # Create img.
        img_uri = urljoin(baseuri, "img")
        scene["img"] = soma.Collection.create(img_uri)

        # Parameters for the multiscale image.
        transform = soma.ScaleTransform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("x", "y"),
            scale_factors=[-1, 1],
        )

        # Cannot set transform before the scene coordinate space is set.
        with pytest.raises(soma.SOMAError):
            scene.add_new_multiscale_image(
                "msi",
                "img",
                transform=transform,
                type=pa.int64(),
                level_shape=[1, 2, 3],
            )

            # The scene coordinate space must be set before registering
            scene.set_transform_to_multiscale_image("msi", transform=transform)

        # Set the scene multiscale image.
        scene_coord_space = soma.CoordinateSpace(
            [soma.Axis(name="x_scene"), soma.Axis(name="y_scene")]
        )
        scene.coordinate_space = scene_coord_space

        # Mismatch in transform input axes and scene coordinate space axes.
        bad_transform = soma.ScaleTransform(
            input_axes=("xbad", "ybad"),
            output_axes=("x", "y"),
            scale_factors=[-1, 1],
        )
        with pytest.raises(ValueError):
            scene.add_new_multiscale_image(
                "msi",
                "img",
                transform=bad_transform,
                type=pa.int64(),
                level_shape=[1, 2, 3],
            )

        # Mismatch in transform output axes and multiscale image coordinate space axes.
        bad_transform = soma.ScaleTransform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("xbad", "ybad"),
            scale_factors=[-1, 1],
        )
        with pytest.raises(ValueError):
            scene.add_new_multiscale_image(
                "msi",
                "img",
                transform=bad_transform,
                type=pa.int64(),
                level_shape=[1, 2, 3],
            )

        # Add the multiscale image.
        scene.add_new_multiscale_image(
            "msi",
            "img",
            transform=transform,
            type=pa.int64(),
            level_shape=[1, 2, 3],
        )

        # Check the transform.
        msi_transform = scene.get_transform_to_multiscale_image("msi")
        assert_transform_equal(msi_transform, transform)


@pytest.mark.parametrize(
    "coord_transform, transform_kwargs",
    [
        (soma.AffineTransform, {"matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}),
        (soma.ScaleTransform, {"scale_factors": [1, 1]}),
        (soma.UniformScaleTransform, {"scale": 1}),
        (soma.IdentityTransform, {}),
    ],
)
@pytest.mark.parametrize("set_coord_space", [True, False])
def test_scene_set_transfrom_to_multiscale_image(
    tmp_path, coord_transform, transform_kwargs, set_coord_space
):
    baseuri = urljoin(
        f"{tmp_path.as_uri()}/", "test_scene_set_transform_to_multiscale_image"
    )

    with soma.Scene.create(baseuri) as scene:
        obsl_uri = urljoin(baseuri, "obsl")
        scene["obsl"] = soma.Collection.create(obsl_uri)

        img_uri = urljoin(baseuri, "img")
        scene["img"] = soma.Collection.create(img_uri)

        coord_space = soma.CoordinateSpace(
            [soma.Axis(name="x_scene"), soma.Axis(name="y_scene")]
        )

        # TODO Add transform directly to add_new_multiscale_image
        scene.add_new_multiscale_image(
            "msi",
            "img",
            transform=None,
            type=pa.int64(),
            level_shape=[3, 8, 9],
        )

        transform = coord_transform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("x", "y"),
            **transform_kwargs,
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_multiscale_image("msi", transform=transform)

        scene.coordinate_space = coord_space

        # No MultiscaleImage named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_multiscale_image("bad", transform=transform)

        # Transform not set
        with pytest.raises(KeyError):
            scene.get_transform_to_multiscale_image("msi")

        # Not a MultiscaleImage
        scene["img"]["col"] = soma.Collection.create(urljoin(img_uri, "col"))
        with pytest.raises(TypeError):
            scene.set_transform_to_multiscale_image("col", transform=transform)

        # Mismatched input axes.
        transform_bad = coord_transform(
            input_axes=("x", "y"),
            output_axes=("x", "y"),
            **transform_kwargs,
        )
        with pytest.raises(ValueError):
            scene.set_transform_to_multiscale_image("msi", transform=transform_bad)

        # Mismatched output axes.
        transform_bad = coord_transform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("x_scene", "y_scene"),
            **transform_kwargs,
        )
        with pytest.raises(ValueError):
            scene.set_transform_to_multiscale_image("msi", transform=transform_bad)

        if set_coord_space:
            bad_coord_space = soma.CoordinateSpace.from_axis_names(("xbad", "ybad"))
            with pytest.raises(ValueError):
                scene.set_transform_to_multiscale_image(
                    "msi", transform=transform, coordinate_space=bad_coord_space
                )

            coord_space = soma.CoordinateSpace(
                (soma.Axis(name="x", unit="nm"), soma.Axis(name="y", unit="nm"))
            )

            msi = scene.set_transform_to_multiscale_image(
                "msi", transform=transform, coordinate_space=coord_space
            )
            actual_coord_space = msi.coordinate_space
            assert actual_coord_space == coord_space

        else:
            msi = scene.set_transform_to_multiscale_image("msi", transform=transform)

        msi_transform = scene.get_transform_to_multiscale_image("msi")
        assert_transform_equal(msi_transform, transform)

        inv_transform = transform.inverse_transform()
        msi_transform = scene.get_transform_from_multiscale_image("msi")
        assert_transform_equal(msi_transform, inv_transform)

        # Set a level to test get transform with level.
        # -- Original size: (3, 8, 9)
        # -- This level: (3, 4, 3)
        # -- x_scale = 3 / 9 = 1 / 3
        # -- y_scale = 4 / 8 = 0.5
        scale_transform = soma.ScaleTransform(
            input_axes=("x", "y"),
            output_axes=("x", "y"),
            scale_factors=[1 / 3, 0.5],
        )
        msi.add_new_level("lowres", shape=(3, 4, 3))

        # Check the transform to the "lowres" level.
        transform_to_level = scale_transform @ transform
        msi_transform = scene.get_transform_to_multiscale_image("msi", level="lowres")
        assert_transform_equal(msi_transform, transform_to_level)

        # Check the transform from the "lowres" level.
        transform_from_level = transform_to_level.inverse_transform()
        msi_transform = scene.get_transform_from_multiscale_image("msi", level="lowres")
        assert_transform_equal(msi_transform, transform_from_level)


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
        coord_space = soma.CoordinateSpace(
            [soma.Axis(name="x_scene"), soma.Axis(name="y_scene")]
        )

        # TODO replace with Scene.add_new_geometry_dataframe when implemented
        scene["obsl"]["gdf"] = soma.GeometryDataFrame.create(gdf_uri, schema=asch)

        transform = coord_transform(
            input_axes=("x_scene", "y_scene"),
            output_axes=("x", "y"),
            **transform_kwargs,
        )

        # The scene coordinate space must be set before registering
        with pytest.raises(soma.SOMAError):
            scene.set_transform_to_geometry_dataframe("gdf", transform=transform)

        scene.coordinate_space = coord_space

        # No SOMAObject named 'bad' in Scene
        with pytest.raises(KeyError):
            scene.set_transform_to_geometry_dataframe("bad", transform=transform)

        # Not a GeometryDataFrame
        scene["obsl"]["col"] = soma.Collection.create(urljoin(obsl_uri, "col"))
        with pytest.raises(typeguard.TypeCheckError):
            scene.set_transform_to_geometry_dataframe("col", transform=transform)

        # Transform not set
        with pytest.raises(KeyError):
            scene.get_transform_to_geometry_dataframe("gdf")

        scene.set_transform_to_geometry_dataframe("gdf", transform=transform)

        gdf_transform = scene.get_transform_to_geometry_dataframe("gdf")
        assert_transform_equal(gdf_transform, transform=transform)
