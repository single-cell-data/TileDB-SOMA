# Copyright (c) 2024 TileDB, Inc,
#
# Licensed under the MIT License.

"""Experimental ingestion methods.

This module contains experimental methods to generate Spatial SOMA artifacts
start from other formats.

Do NOT merge into main.
"""

import json
import os
import pathlib
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
)

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scanpy
from anndata import AnnData
from PIL import Image

from .. import (
    Axis,
    Collection,
    CompositeTransform,
    CoordinateSystem,
    DataFrame,
    DenseNDArray,
    Experiment,
    ScaleTransform,
    Scene,
    SparseNDArray,
)
from .._constants import SOMA_JOINID
from .._tiledb_object import AnyTileDBObject
from .._types import IngestMode
from ..io import from_anndata
from ..io.ingest import (
    IngestCtx,
    IngestionParams,
    _create_or_open_collection,
    _maybe_set,
    _write_dataframe_impl,
)

if TYPE_CHECKING:
    from somacore.options import PlatformConfig

    from .._types import Path
    from ..io._registration import ExperimentAmbientLabelMapping
    from ..io.ingest import AdditionalMetadata
    from ..options import SOMATileDBContext


def from_cxg_spatial_h5ad(
    input_h5ad_path: pathlib.Path,
    experiment_uri: str,
    measurement_name: str,
    scene_name: str,
    *,
    context: Optional["SOMATileDBContext"] = None,
    platform_config: Optional["PlatformConfig"] = None,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    X_layer_name: str = "data",
    raw_X_layer_name: str = "data",
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
    X_kind: Union[Type[SparseNDArray], Type[DenseNDArray]] = SparseNDArray,
    registration_mapping: Optional["ExperimentAmbientLabelMapping"] = None,
    uns_keys: Optional[Sequence[str]] = None,
    additional_metadata: "AdditionalMetadata" = None,
):
    """
    This function reads cellxgene schema compliant H5AD file and writes
    to the SOMA data format rooted at the given `Experiment` URI.

    NOTE:
    1. This function is for testing/prototyping purposes only
    2. This function expects an h5ad that has spatial data therefore must be
    compliant with the following cellxgene schema 5.1.0:
    https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.1.0/schema.md

    Lifecycle:
        Experimental
    """

    if ingest_mode != "write":
        raise NotImplementedError(
            f'the only ingest_mode currently supported is "write"; got "{ingest_mode}"'
        )

    adata = ad.read_h5ad(input_h5ad_path)

    spatial_dict = adata.uns["spatial"]
    if not spatial_dict["is_single"]:
        raise NotImplementedError(
            "Only spatial datasets with uns['spatial']['is_single'] == True are supported"
        )

    # create a directory to store some intermediate files
    filename, _ = os.path.splitext(os.path.basename(input_h5ad_path))
    spatial_assets_dir = f"{filename}/spatial"
    os.makedirs(spatial_assets_dir, exist_ok=True)

    # get scale_factors dictionary
    library_id = (spatial_dict.keys() - {"is_single"}).pop()
    scale_factors = spatial_dict[library_id]["scalefactors"]

    # store tissue_postions.csv
    tissue_pos_df = pd.DataFrame(adata.obsm["spatial"])
    tissue_pos_df.index = adata.obs.index
    tissue_positions_file_path = f"{spatial_assets_dir}/tissue_positions.csv"
    tissue_pos_df.to_csv(
        tissue_positions_file_path,
        index_label="barcode",
        header=["pxl_col_in_fullres", "pxl_row_in_fullres"],
    )

    # store spatial images
    images_source_dict = spatial_dict[library_id]["images"]
    image_paths = {"fullres": None, "hires": None, "lowres": None}

    for key, img_array in images_source_dict.items():
        if img_array.dtype == np.float32:
            img_array = (img_array * 255).astype(np.uint8)
        img = Image.fromarray(img_array)
        img_file_path = f"{spatial_assets_dir}/{key}.png"
        img.save(img_file_path)
        image_paths[key] = img_file_path

    return _write_visium_data_to_experiment_uri(
        experiment_uri=experiment_uri,
        measurement_name=measurement_name,
        scene_name=scene_name,
        adata=adata,
        obs_id_name=obs_id_name,
        var_id_name=var_id_name,
        X_layer_name=X_layer_name,
        raw_X_layer_name=raw_X_layer_name,
        X_kind=X_kind,
        uns_keys=uns_keys,
        additional_metadata=additional_metadata,
        registration_mapping=registration_mapping,
        ingest_mode=ingest_mode,
        use_relative_uri=use_relative_uri,
        context=context,
        platform_config=platform_config,
        scale_factors=scale_factors,
        input_hires=image_paths["hires"],
        input_lowres=image_paths["lowres"],
        input_fullres=image_paths["fullres"],
        input_tissue_positions=tissue_positions_file_path,
    )


def from_visium(
    experiment_uri: str,
    input_path: "Path",
    measurement_name: str,
    scene_name: str,
    *,
    context: Optional["SOMATileDBContext"] = None,
    platform_config: Optional["PlatformConfig"] = None,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    X_layer_name: str = "data",
    raw_X_layer_name: str = "data",
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
    X_kind: Union[Type[SparseNDArray], Type[DenseNDArray]] = SparseNDArray,
    registration_mapping: Optional["ExperimentAmbientLabelMapping"] = None,
    uns_keys: Optional[Sequence[str]] = None,
    additional_metadata: "AdditionalMetadata" = None,
    use_raw_counts: bool = True,
) -> str:
    """Reads a 10x Visium dataset and writes it to an :class:`Experiment`.

    This function is for ingesting Visium data for prototyping and testing the
    proposed spatial design.

    TODO: Args list

    WARNING: This was only tested for Space Ranger version 2 output.

    Lifecycle:
        Experimental
    """

    if ingest_mode != "write":
        raise NotImplementedError(
            f'the only ingest_mode currently supported is "write"; got "{ingest_mode}"'
        )

    # Get input file locations.
    input_path = pathlib.Path(input_path)

    input_gene_expression = (
        input_path / "raw_feature_bc_matrix.h5"
        if use_raw_counts
        else input_path / "filtered_feature_bc_matrix.h5"
    )

    # TODO: Generalize - this is hard-coded for Space Ranger version 2
    input_tissue_positions = input_path / "spatial/tissue_positions.csv"

    # Get JSON scale factors.
    input_scale_factors = input_path / "spatial/scalefactors_json.json"
    with open(input_scale_factors, mode="r", encoding="utf-8") as scale_factors_json:
        scale_factors = json.load(scale_factors_json)

    # TODO: Generalize - hard-coded for Space Ranger version 2
    input_hires = input_path / "spatial/tissue_hires_image.png"
    input_lowres = input_path / "spatial/tissue_lowres_image.png"
    input_fullres = None

    adata = scanpy.read_10x_h5(input_gene_expression)

    return _write_visium_data_to_experiment_uri(
        experiment_uri=experiment_uri,
        measurement_name=measurement_name,
        scene_name=scene_name,
        adata=adata,
        obs_id_name=obs_id_name,
        var_id_name=var_id_name,
        X_layer_name=X_layer_name,
        raw_X_layer_name=raw_X_layer_name,
        X_kind=X_kind,
        uns_keys=uns_keys,
        additional_metadata=additional_metadata,
        registration_mapping=registration_mapping,
        ingest_mode=ingest_mode,
        use_relative_uri=use_relative_uri,
        context=context,
        platform_config=platform_config,
        scale_factors=scale_factors,
        input_hires=input_hires,
        input_lowres=input_lowres,
        input_fullres=input_fullres,
        input_tissue_positions=input_tissue_positions,
    )


def _write_visium_data_to_experiment_uri(
    experiment_uri: str,
    measurement_name: str,
    scene_name: str,
    *,
    adata: AnnData,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    X_layer_name: str = "data",
    raw_X_layer_name: str = "data",
    X_kind: Union[Type[SparseNDArray], Type[DenseNDArray]] = SparseNDArray,
    uns_keys: Optional[Sequence[str]] = None,
    additional_metadata: "AdditionalMetadata" = None,
    registration_mapping: Optional["ExperimentAmbientLabelMapping"] = None,
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
    context: Optional["SOMATileDBContext"] = None,
    platform_config: Optional["PlatformConfig"] = None,
    scale_factors: Dict[str, Any],
    input_hires: Optional[pathlib.Path],
    input_lowres: Optional[pathlib.Path],
    input_fullres: Optional[pathlib.Path],
    input_tissue_positions: pathlib.Path,
) -> str:
    uri = from_anndata(
        experiment_uri,
        adata,
        measurement_name,
        context=context,
        platform_config=platform_config,
        obs_id_name=obs_id_name,
        var_id_name=var_id_name,
        X_layer_name=X_layer_name,
        raw_X_layer_name=raw_X_layer_name,
        ingest_mode=ingest_mode,
        use_relative_uri=use_relative_uri,
        X_kind=X_kind,
        registration_mapping=registration_mapping,
        uns_keys=uns_keys,
        additional_metadata=additional_metadata,
    )

    ingest_ctx: IngestCtx = {
        "context": context,
        "ingestion_params": IngestionParams(ingest_mode, registration_mapping),
        "additional_metadata": additional_metadata,
    }

    pixels_per_spot_radius = 0.5 * scale_factors["spot_diameter_fullres"]
    fullres_to_coords_scale = 65 / scale_factors["spot_diameter_fullres"]

    # Create axes and transformations
    coordinate_system = CoordinateSystem(
        (
            Axis(axis_name="y", axis_type="space", axis_unit="micrometer"),
            Axis(axis_name="x", axis_type="space", axis_unit="micrometer"),
        )
    )

    spots_to_coords = CompositeTransform(
        (ScaleTransform((fullres_to_coords_scale, fullres_to_coords_scale)),)
    )

    # TODO: The `obs_df` should be in dataframe with only soma_joinid and obs_id. Not
    # currently bothering to check/enforce this.
    with Experiment.open(uri, mode="r", context=context) as experiment:
        obs_df = experiment.obs.read().concat().to_pandas()

    # Add spatial information to the experiment.
    with Experiment.open(uri, mode="w", context=context) as experiment:
        spatial_uri = f"{uri}/spatial"
        with _create_or_open_collection(
            Collection[Union[DataFrame, Scene]], spatial_uri, **ingest_ctx
        ) as spatial:
            _maybe_set(
                experiment, "spatial", spatial, use_relative_uri=use_relative_uri
            )
            scene_uri = f"{spatial_uri}/{scene_name}"
            with _create_or_open_collection(Scene, scene_uri, **ingest_ctx) as scene:
                _maybe_set(
                    spatial, scene_name, scene, use_relative_uri=use_relative_uri
                )
                scene.metadata["soma_scene_coordinates"] = coordinate_system.to_json()

                # Write image data and add to the scene.
                images_uri = f"{scene_uri}/img"
                with _write_visium_images(
                    images_uri,
                    scale_factors,
                    input_hires=input_hires,
                    input_lowres=input_lowres,
                    input_fullres=input_fullres,
                    use_relative_uri=use_relative_uri,
                    **ingest_ctx,
                ) as images:
                    images.metadata["soma_scene_coords"] = spots_to_coords.to_json()
                    _maybe_set(scene, "img", images, use_relative_uri=use_relative_uri)
                scene.metadata.update(
                    {"soma_asset_transform_images": spots_to_coords.to_json()}
                )

                obsl_uri = f"{scene_uri}/obsl"

                # Write spot data and add to the scene.
                with _write_visium_spot_dataframe(
                    obsl_uri,
                    input_tissue_positions,
                    pixels_per_spot_radius,
                    obs_df,
                    obs_id_name,
                    **ingest_ctx,
                ) as obsl:
                    _maybe_set(scene, "obsl", obsl, use_relative_uri=use_relative_uri)
                    obsl.metadata.update(
                        {"soma_asset_transform_loc": spots_to_coords.to_json()}
                    )

                varl_uri = f"{scene_uri}/varl"
                with _create_or_open_collection(
                    Collection[Collection[AnyTileDBObject]], varl_uri, **ingest_ctx
                ) as varl:
                    _maybe_set(scene, "varl", varl, use_relative_uri=use_relative_uri)
        return uri


def _write_visium_spot_dataframe(
    df_uri: str,
    input_tissue_positions: pathlib.Path,
    spot_radius: float,
    obs_df: pd.DataFrame,
    id_column_name: str,
    *,
    ingestion_params: IngestionParams,
    additional_metadata: "AdditionalMetadata" = None,
    platform_config: Optional["PlatformConfig"] = None,
    context: Optional["SOMATileDBContext"] = None,
) -> DataFrame:
    """TODO: Add _write_visium_spot_dataframe docs"""
    df = (
        pd.read_csv(input_tissue_positions)
        .rename(
            columns={
                "barcode": id_column_name,
                "pxl_col_in_fullres": "y",
                "pxl_row_in_fullres": "x",
            }
        )
        .assign(_soma_geometry=np.double(spot_radius))
    )
    df = pd.merge(obs_df, df, how="inner", on=id_column_name)
    return _write_dataframe_impl(
        df,
        df_uri,
        id_column_name,
        ingestion_params=ingestion_params,
        additional_metadata=additional_metadata,
        index_column_names=("y", "x", SOMA_JOINID),
        platform_config=platform_config,
        context=context,
    )


def _write_visium_images(
    uri: str,
    scale_factors: Dict[str, Any],
    *,
    input_hires: Optional[pathlib.Path],
    input_lowres: Optional[pathlib.Path],
    input_fullres: Optional[pathlib.Path],
    ingestion_params: IngestionParams,
    additional_metadata: "AdditionalMetadata" = None,
    platform_config: Optional["PlatformConfig"] = None,
    context: Optional["SOMATileDBContext"] = None,
    use_relative_uri: Optional[bool] = None,
) -> Collection[DenseNDArray]:
    input_images: Dict[str, Tuple[pathlib.Path, List[float]]] = {}
    if input_fullres is not None:
        input_images["fullres"] = (input_fullres, [1.0, 1.0, 1.0])
    if input_hires is not None:
        scale = 1.0 / scale_factors["tissue_hires_scalef"]
        input_images["hires"] = (input_hires, [1.0, scale, scale])
    if input_lowres is not None:
        scale = 1.0 / scale_factors["tissue_lowres_scalef"]
        input_images["lowres"] = (input_lowres, [1.0, scale, scale])
    axes_metadata = [
        {"name": "c", "type": "channel"},
        {"name": "y", "type": "space", "unit": "micrometer"},
        {"name": "x", "type": "space", "unit": "micrometer"},
    ]
    return _write_multiscale_images(
        uri,
        input_images,
        axes_metadata=axes_metadata,
        ingestion_params=ingestion_params,
        additional_metadata=additional_metadata,
        platform_config=platform_config,
        context=context,
        use_relative_uri=use_relative_uri,
    )


def _write_multiscale_images(
    uri: str,
    input_images: Dict[str, Tuple[pathlib.Path, List[float]]],
    *,
    axes_metadata: List[Dict[str, str]],
    ingestion_params: IngestionParams,
    additional_metadata: "AdditionalMetadata" = None,
    platform_config: Optional["PlatformConfig"] = None,
    context: Optional["SOMATileDBContext"] = None,
    use_relative_uri: Optional[bool] = None,
) -> Collection[DenseNDArray]:
    """TODO: Write full docs for this function

    TODO: Need to add in collection level metadata. In this case it will be

    """
    collection = _create_or_open_collection(
        Collection[DenseNDArray],
        uri,
        ingestion_params=ingestion_params,
        additional_metadata=additional_metadata,
        context=context,
    )
    dataset_metadata = {}
    for image_name, (image_path, image_scales) in input_images.items():
        dataset_metadata[f"soma_asset_transform_{image_name}"] = CompositeTransform(
            (ScaleTransform(image_scales),)
        ).to_json()
        image_uri = f"{uri}/{image_name}"

        # TODO: Need to create new imaging type with dimensions 'c', 'y', 'x'
        im = np.transpose(np.array(Image.open(image_path)), (2, 0, 1))
        image_array = DenseNDArray.create(
            image_uri,
            type=pa.from_numpy_dtype(im.dtype),
            shape=im.shape,
            platform_config=platform_config,
            context=context,
        )
        tensor = pa.Tensor.from_numpy(im)
        image_array.write(
            (slice(None), slice(None), slice(None)),
            tensor,
            platform_config=platform_config,
        )
        _maybe_set(
            collection, image_name, image_array, use_relative_uri=use_relative_uri
        )
    collection.metadata.update(dataset_metadata)
    return collection
