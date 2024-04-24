# Copyright (c) 2024 TileDB, Inc,
#
# Licensed under the MIT License.

"""Experimental ingestion methods.

This module contains experimental methods to generate Spatial SOMA artifacts
start from other formats.

Do NOT merge into main.
"""

import json
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

import numpy as np
import pandas as pd
import pyarrow as pa
import scanpy
from PIL import Image

from .. import Collection, DataFrame, DenseNDArray, Experiment, Scene, SparseNDArray
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
    from ..io.ingeset import AdditionalMetadata
    from ..options import SOMATileDBContext


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
    input_scale_factors = input_path / "spatial/scalefactors_json.json"

    # TODO: Generalize - hard-coded for Space Ranger version 2
    input_hires = input_path / "spatial/tissue_hires_image.png"
    input_lowres = input_path / "spatial/tissue_lowres_image.png"
    input_fullres = None

    # Create the
    anndata = scanpy.read_10x_h5(input_gene_expression)
    uri = from_anndata(
        experiment_uri,
        anndata,
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

    # Get JSON scale factors.
    with open(input_scale_factors, mode="r", encoding="utf-8") as scale_factors_json:
        scale_factors = json.load(scale_factors_json)

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
            with _create_or_open_collection(
                Collection[Collection[AnyTileDBObject]], scene_uri, **ingest_ctx
            ) as scene:
                _maybe_set(
                    spatial, scene_name, scene, use_relative_uri=use_relative_uri
                )

                scene_exp_uri = f"{scene_uri}/exp"
                with _create_or_open_collection(
                    Collection[AnyTileDBObject], scene_exp_uri, **ingest_ctx
                ) as scene_exp:
                    _maybe_set(
                        scene, "exp", scene_exp, use_relative_uri=use_relative_uri
                    )

                    obs_locations_uri = f"{scene_exp_uri}/obs_locations"

                    # Write spot data and add to the scene.
                    with _write_visium_spot_dataframe(
                        obs_locations_uri,
                        input_tissue_positions,
                        scale_factors,
                        obs_df,
                        obs_id_name,
                        **ingest_ctx,
                    ) as obs_locations:
                        _maybe_set(
                            scene_exp,
                            "obs_locations",
                            obs_locations,
                            use_relative_uri=use_relative_uri,
                        )

                    # Write image data and add to the scene.
                    images_uri = f"{scene_exp_uri}/images"
                    with _write_visium_images(
                        images_uri,
                        scale_factors,
                        input_hires=input_hires,
                        input_lowres=input_lowres,
                        input_fullres=input_fullres,
                        use_relative_uri=use_relative_uri,
                        **ingest_ctx,
                    ) as images:
                        _maybe_set(
                            scene_exp,
                            "images",
                            images,
                            use_relative_uri=use_relative_uri,
                        )

                scene_ms_uri = f"{scene_uri}/ms"
                with _create_or_open_collection(
                    Collection[Collection[AnyTileDBObject]], scene_ms_uri, **ingest_ctx
                ) as scene_ms:
                    _maybe_set(scene, "ms", scene_ms, use_relative_uri=use_relative_uri)
        return uri


def _write_visium_spot_dataframe(
    df_uri: str,
    input_tissue_positions: pathlib.Path,
    scale_factors: Dict[str, Any],
    obs_df: pd.DataFrame,
    id_column_name: str,
    *,
    ingestion_params: IngestionParams,
    additional_metadata: "AdditionalMetadata" = None,
    platform_config: Optional["PlatformConfig"] = None,
    context: Optional["SOMATileDBContext"] = None,
) -> DataFrame:
    """TODO: Add _write_visium_spot_dataframe docs"""
    # Create the
    spot_radius = 0.5 * scale_factors["spot_diameter_fullres"]
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
    datasets_metadata = []
    for image_name, (image_path, image_scales) in input_images.items():
        datasets_metadata.append(
            {
                "path": image_name,
                "coordinateTransforms": [{"type": "scale", "scale": image_scales}],
            }
        )
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
    metadata_blob = json.dumps(
        {
            "multiscales": [
                {
                    "version": "0.1.0-dev",
                    "name": "visium-example",
                    "datasets": datasets_metadata,
                }
            ]
        }
    )
    collection.metadata.update({"multiscales": metadata_blob})
    return collection
