import itertools
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

import tiledbsoma as soma

from .test_experiment_query import add_dataframe, add_sparse_array

shapely = pytest.importorskip("shapely")
sd = pytest.importorskip("spatialdata")


def add_multiscale_image(
    scene: soma.Scene,
    key: str,
    shapes: Tuple[Tuple[int, ...]],
):
    with scene.add_new_multiscale_image(
        key,
        "img",
        type=pa.uint8(),
        level_shape=shapes[0],
        transform=soma.IdentityTransform(("x", "y"), ("x", "y")),
    ) as image:
        coords = (slice(None), slice(None), slice(None))
        image["level0"].write(
            coords,
            pa.Tensor.from_numpy(
                np.random.randint(0, 255, size=shapes[0], dtype=np.uint8)
            ),
        )

        for index, level_shape in enumerate(shapes[1:]):
            level = image.add_new_level(f"level{index+1}", shape=level_shape)
            level.write(
                coords,
                pa.Tensor.from_numpy(
                    np.random.randint(0, 255, size=level_shape, dtype=np.uint8)
                ),
            )


def add_presence_dataframe(
    coll: soma.Collection,
    key: str,
    max_joinids: int,
    joinids: np.ndarray,
    scene_ids: List[str],
) -> None:
    df = coll.add_new_dataframe(
        key,
        schema=pa.schema(
            [
                ("soma_joinid", pa.int64()),
                ("scene_id", pa.string()),
                ("data", pa.bool_()),
            ]
        ),
        domain=((0, max_joinids - 1), ("", "")),
        index_column_names=("soma_joinid", "scene_id"),
    )
    df.write(
        pa.Table.from_pydict(
            {
                "soma_joinid": joinids,
                "scene_id": scene_ids,
                "data": len(scene_ids) * [True],
            }
        )
    )


def add_point_cloud_dataframe(
    scene: soma.Scene,
    subcoll: Union[str, List[str]],
    key: str,
    data: Dict[str, np.ndarray],
    circles: bool,
):
    with scene.add_new_point_cloud_dataframe(
        key,
        subcoll,
        transform=soma.IdentityTransform(("x", "y"), ("x", "y")),
        schema=pa.schema([("x", pa.float64()), ("y", pa.float64())]),
        domain=[[-1, 1], [-1, 1], [0, 63]],
    ) as points:

        if circles:
            points.metadata["soma_geometry"] = 2.0
            points.metadata["soma_geometry_type"] = "radius"

        points.write(pa.Table.from_pydict(data))


def add_scene(
    coll: soma.Collection,
    key: str,
    *,
    points: Dict[Tuple[Union[str, Tuple[str, ...]], str], Dict[str, np.ndarray]],
    circles: Dict[Tuple[Union[str, Tuple[str, ...]], str], Dict[str, np.ndarray]],
    images: Dict[str, Tuple[Tuple[int, ...]]],
) -> None:
    with coll.add_new_collection(key, soma.Scene, coordinate_space=("x", "y")) as scene:
        scene.add_new_collection("obsl")
        varl = scene.add_new_collection("varl")
        varl.add_new_collection("RNA")
        varl.add_new_collection("other")
        scene.add_new_collection("img")

        for (subcoll, key), data in points.items():
            add_point_cloud_dataframe(
                scene,
                subcoll if isinstance(subcoll, str) else list(subcoll),
                key,
                data,
                circles=False,
            )

        for (subcoll, key), data in circles.items():
            add_point_cloud_dataframe(
                scene,
                subcoll if isinstance(subcoll, str) else list(subcoll),
                key,
                data,
                circles=True,
            )

        for key, shapes in images.items():
            add_multiscale_image(scene, key, shapes)


@pytest.fixture(scope="module")
def soma_spatial_experiment(tmp_path_factory) -> soma.Experiment:
    """Creates a multi-scene experiment where the scenes contain disjoint observations."""
    uri = (tmp_path_factory.mktemp("base") / "exp").as_uri()

    n_obs = 64
    n_vars = 100

    with soma.Experiment.create(uri) as exp:

        # Add spatial data.
        spatial = exp.add_new_collection("spatial")
        for index in range(4):
            x, y = np.meshgrid(
                np.linspace(-1.0, 1.0, num=4), np.linspace(-1.0, 1.0, num=4)
            )
            point_df = {
                "x": x.flatten()[:8],
                "y": y.flatten()[:8],
                "soma_joinid": np.arange(index * 16, index * 16 + 8, dtype=np.int64),
            }
            circle_df = {
                "x": x.flatten()[8:],
                "y": y.flatten()[8:],
                "soma_joinid": np.arange(
                    index * 16 + 8, (index + 1) * 16, dtype=np.int64
                ),
            }
            add_scene(
                spatial,
                f"scene{index}",
                points={
                    ("obsl", "points1"): point_df,
                    (("varl", "RNA"), "points2"): point_df,
                    (("varl", "other"), "points3"): point_df,
                },
                circles={
                    ("obsl", "shapes1"): circle_df,
                    (("varl", "RNA"), "shapes2"): circle_df,
                    (("varl", "other"), "shapes3"): circle_df,
                },
                images={"tissue": ((3, 16, 8),)},
            )

        # Create obs dataframe and obs/spatial presence matrix.
        add_dataframe(exp, "obs", n_obs)
        add_presence_dataframe(
            exp,
            "obs_spatial_presence",
            n_obs,
            np.arange(n_obs),
            list(itertools.chain(*(16 * [f"scene{id}"] for id in range(4)))),
        )

        # Add measurement data.
        ms = exp.add_new_collection("ms")
        rna = ms.add_new_collection("RNA", soma.Measurement)
        add_dataframe(rna, "var", n_vars)
        rna_x = rna.add_new_collection("X", soma.Collection)
        add_sparse_array(rna_x, "data", (n_obs, n_vars))

        # Add var/spatial presence matrix that specifies all scenes have variables 20-80
        var_ids, scene_ids = np.meshgrid(
            np.arange(20, 81), np.array(["scene0", "scene1", "scene2", "scene3"])
        )
        add_presence_dataframe(
            rna,
            "var_spatial_presence",
            n_vars,
            var_ids.flatten(),
            scene_ids.flatten().tolist(),
        )

    exp = soma.Experiment.open(uri)

    # Check coordinate spaces.
    for index in range(4):
        scene = exp.spatial[f"scene{index}"]
        assert scene.coordinate_space is not None

    return exp


def check_for_scene_data(sdata, has_scenes: List[bool]):
    for index, has_scene in enumerate(has_scenes):
        if not has_scene:
            continue

        x, y = np.meshgrid(np.linspace(-1.0, 1.0, num=4), np.linspace(-1.0, 1.0, num=4))
        expected_points = pd.DataFrame.from_dict(
            {
                "x": x.flatten()[:8],
                "y": y.flatten()[:8],
                "soma_joinid": np.arange(
                    index * 16, (index + 1) * 16 - 8, dtype=np.int64
                ),
            }
        )
        expected_shapes = pd.DataFrame.from_dict(
            {
                "soma_joinid": np.arange(
                    index * 16 + 8, (index + 1) * 16, dtype=np.int64
                ),
                "radius": 8 * [2.0],
                "geometry": shapely.points(
                    list(zip(x.flatten()[8:], y.flatten()[8:]))
                ).tolist(),
            }
        )

        scene_id = f"scene{index}"
        points1 = sdata.points[f"{scene_id}_points1"]
        assert all(points1 == expected_points)
        points2 = sdata.points[f"{scene_id}_RNA_points2"]
        assert all(points2 == expected_points)
        shapes1 = sdata.shapes[f"{scene_id}_shapes1"]
        assert all(shapes1 == expected_shapes)
        shapes2 = sdata.shapes[f"{scene_id}_RNA_shapes2"]
        assert all(shapes2 == expected_shapes)
        image = sdata.images[f"{scene_id}_tissue"]
        assert image.shape == (3, 16, 8)


def test_spatial_experiment_query_none(soma_spatial_experiment):
    with soma_spatial_experiment.axis_query(
        "RNA",
        obs_query=soma.AxisQuery(value_filter="label == 'no-such-label'"),
        var_query=soma.AxisQuery(value_filter="label == 'no-such-label'"),
    ) as query:
        assert query.n_obs == 0
        assert query.n_vars == 0

        # Read to SpatialData.
        sdata = query.to_spatialdata("data")
        assert len(sdata.tables) == 1
        assert len(sdata.points) == 0
        assert len(sdata.shapes) == 0
        assert len(sdata.images) == 0

        # Check table is empty.
        ad = sdata["RNA"]
        assert ad.n_obs == 0 and ad.n_vars == 0

        # Check no region columns/metadata.
        assert "region_key" not in ad.obs
        assert "instance_key" not in ad.obs
        assert "spatialdata_attrs" not in ad.uns


def test_spatial_experiment_query_all(soma_spatial_experiment):
    with soma_spatial_experiment.axis_query("RNA") as query:
        # Check all obs/var read.
        assert query.n_obs == 64
        assert query.n_vars == 100

        # Read to SpatialData.
        sdata = query.to_spatialdata("data")

        # Check the expected scenes are included.
        scene_ids = set(str(val) for val in query.obs_scene_ids())
        expected_scene_ids = {"scene0", "scene1", "scene2", "scene3"}
        assert scene_ids == expected_scene_ids

        # Verify the correct number of assets.
        assert len(sdata.tables) == 1
        assert len(sdata.points) == 8
        assert len(sdata.shapes) == 8
        assert len(sdata.images) == 4

        check_for_scene_data(sdata, 4 * [True])

        # Check table.
        ad = sdata.tables["RNA"]

        sd_attrs = ad.uns["spatialdata_attrs"]
        assert sd_attrs["region"] == [
            "scene0_points1",
            "scene0_shapes1",
            "scene1_points1",
            "scene1_shapes1",
            "scene2_points1",
            "scene2_shapes1",
            "scene3_points1",
            "scene3_shapes1",
        ]
        assert sd_attrs["region_key"] == "region_key"
        assert sd_attrs["instance_key"] == "instance_key"

        region_df = ad.obs[["soma_joinid", "region_key", "instance_key"]]
        for scene_index in range(4):
            # Filter on points key and get the joinids and instance keys for the points.
            points_key = f"scene{scene_index}_points1"
            points_region_df = region_df[region_df["region_key"] == points_key]
            instance_keys = points_region_df["instance_key"]
            region_joinids = points_region_df["soma_joinid"].tolist()

            # Check the joinids at the points instance key match expected joinids.
            points = sdata.points[points_key].compute()
            points_joinids = points.iloc[instance_keys]["soma_joinid"].tolist()
            assert region_joinids == points_joinids

            # Filter on shape key and get the joinids and instance keys for the shapes.
            shapes_key = f"scene{scene_index}_shapes1"
            shapes_region_df = region_df[region_df["region_key"] == shapes_key]
            instance_keys = shapes_region_df["instance_key"]
            region_joinids = shapes_region_df["soma_joinid"].tolist()

            # Check the joinids at the shapes instance key match expected joinids.
            shapes = sdata.shapes[shapes_key]
            shapes_joinids = shapes.iloc[instance_keys]["soma_joinid"].to_list()
            assert region_joinids == shapes_joinids


@pytest.mark.parametrize(
    "obs_slice,has_scene",
    [
        (slice(0, 15), [True, False, False, False]),
        (slice(16, 31), [False, True, False, False]),
        (slice(32, 47), [False, False, True, False]),
        (slice(48, 63), [False, False, False, True]),
    ],
)
def test_spatial_experiment_query_obs_slice(
    soma_spatial_experiment, obs_slice, has_scene
):
    nscene = has_scene.count(True)
    with soma_spatial_experiment.axis_query(
        "RNA", obs_query=soma.AxisQuery(coords=(obs_slice,))
    ) as query:
        # Check all obs/var read.
        assert query.n_obs == obs_slice.stop - obs_slice.start + 1
        assert query.n_vars == 100

        # Check the expected scenes are included.
        scene_ids = set(str(val) for val in query.obs_scene_ids())
        expected_scene_ids = set(
            f"scene{index}" for index, val in enumerate(has_scene) if val
        )
        assert scene_ids == expected_scene_ids

        # Read to SpatialData.
        sdata = query.to_spatialdata("data")

        # Verify the correct number of assets.
        assert len(sdata.tables) == 1
        assert len(sdata.points) == 2 * nscene
        assert len(sdata.shapes) == 2 * nscene
        assert len(sdata.images) == nscene

        check_for_scene_data(sdata, has_scene)


@pytest.mark.parametrize(
    "var_slice,has_scene",
    [
        (slice(0, 10), [False, False, False, False]),
        (slice(30, 70), [True, True, True, True]),
    ],
)
def test_spatial_experiment_query_var_slice(
    soma_spatial_experiment, var_slice, has_scene
):
    nscene = has_scene.count(True)
    with soma_spatial_experiment.axis_query(
        "RNA", var_query=soma.AxisQuery(coords=(var_slice,))
    ) as query:
        assert query.n_obs == 64
        assert query.n_vars == var_slice.stop - var_slice.start + 1

        # Check the expected scenes are included.
        scene_ids = set(str(val) for val in query.var_scene_ids())
        expected_scene_ids = set(
            f"scene{index}" for index, val in enumerate(has_scene) if val
        )
        assert scene_ids == expected_scene_ids

        # Read to SpatialData.
        sdata = query.to_spatialdata("data", scene_presence_mode="var")

        # Verify the correct number of assets.
        assert len(sdata.tables) == 1
        assert len(sdata.points) == 2 * nscene
        assert len(sdata.shapes) == 2 * nscene
        assert len(sdata.images) == nscene

        check_for_scene_data(sdata, has_scene)
