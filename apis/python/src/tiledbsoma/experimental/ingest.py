# Copyright (c) 2024 TileDB, Inc
#
# Licensed under the MIT License.

"""Experimental ingestion methods.

This module contains experimental methods to generate Spatial SOMA artifacts
start from other formats.

Do NOT merge into main.
"""

import json
import os
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    Any,
    Dict,
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
import pyarrow.compute as pacomp
import scanpy
from anndata import AnnData
from PIL import Image
from somacore import Axis, CoordinateSpace, IdentityTransform

from .. import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    MultiscaleImage,
    PointCloud,
    Scene,
    SparseNDArray,
    _util,
    logging,
)
from .._arrow_types import df_to_arrow
from .._exception import (
    AlreadyExistsError,
    NotCreateableError,
    SOMAError,
)
from .._soma_object import AnySOMAObject
from .._types import IngestMode
from ..io import from_anndata
from ..io.ingest import (
    IngestCtx,
    IngestionParams,
    _create_or_open_collection,
    _maybe_set,
    _write_arrow_table,
    add_metadata,
)
from ..options._tiledb_create_write_options import (
    TileDBCreateOptions,
    TileDBWriteOptions,
)

if TYPE_CHECKING:
    from somacore.options import PlatformConfig

    from ..io._common import AdditionalMetadata
    from ..io._registration import ExperimentAmbientLabelMapping
    from ..options import SOMATileDBContext


def from_cxg_spatial_h5ad(
    input_h5ad_path: Path,
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
    write_obs_spatial_presence: bool = True,
    write_var_spatial_presence: bool = False,
) -> str:
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
    image_paths: Dict[str, Optional[str]] = {
        "fullres": None,
        "hires": None,
        "lowres": None,
    }

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
        input_tissue_positions=Path(tissue_positions_file_path),
        write_obs_spatial_presence=write_obs_spatial_presence,
        write_var_spatial_presence=write_var_spatial_presence,
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
    write_obs_spatial_presence: bool = False,
    write_var_spatial_presence: bool = False,
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
    input_path = Path(input_path)

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
        write_obs_spatial_presence=write_obs_spatial_presence,
        write_var_spatial_presence=write_var_spatial_presence,
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
    input_hires: Union[None, str, Path],
    input_lowres: Union[None, str, Path],
    input_fullres: Union[None, str, Path],
    input_tissue_positions: Path,
    write_obs_spatial_presence: bool = False,
    write_var_spatial_presence: bool = False,
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

    pixels_per_spot_diameter = scale_factors["spot_diameter_fullres"]

    # Get the size of the fullres image.
    # TODO: This is a hack.
    if input_fullres is not None:
        with Image.open(input_fullres) as im:
            ref_shape: Tuple[int, ...] = np.array(im).shape
    elif input_hires is not None:
        with Image.open(input_hires) as im:
            width, height, nchannel = np.array(im).shape  # type: ignore
        scale = scale_factors["tissue_hires_scalef"]
        ref_shape = (
            int(np.round(width / scale)),
            int(np.round(height / scale)),
            nchannel,
        )
    elif input_lowres is not None:
        with Image.open(input_lowres) as im:
            width, height, nchannel = np.array(im).shape  # type: ignore
        scale = scale_factors["tissue_lowres_scalef"]
        ref_shape = (
            int(np.round(width / scale)),
            int(np.round(height / scale)),
            nchannel,
        )

    # Create axes and transformations
    coord_space = CoordinateSpace(
        (Axis(name="x", units="pixels"), Axis(name="y", units="pixels"))
    )

    with Experiment.open(uri, mode="r", context=context) as exp:
        obs_df = (
            exp.obs.read(column_names=["soma_joinid", "obs_id"]).concat().to_pandas()
        )
        if write_obs_spatial_presence or write_var_spatial_presence:
            x_layer = exp.ms[measurement_name].X[X_layer_name].read().tables().concat()
            if write_obs_spatial_presence:
                obs_id = pacomp.unique(x_layer["soma_dim_0"])
            if write_var_spatial_presence:
                var_id = pacomp.unique(x_layer["soma_dim_1"])

    # Add spatial information to the experiment.
    with Experiment.open(uri, mode="w", context=context) as exp:
        spatial_uri = f"{uri}/spatial"
        with _create_or_open_collection(
            Collection[Scene], spatial_uri, **ingest_ctx
        ) as spatial:
            _maybe_set(exp, "spatial", spatial, use_relative_uri=use_relative_uri)
            scene_uri = f"{spatial_uri}/{scene_name}"
            with _create_or_open_collection(Scene, scene_uri, **ingest_ctx) as scene:
                _maybe_set(
                    spatial, scene_name, scene, use_relative_uri=use_relative_uri
                )
                scene.coordinate_space = coord_space

                # Write image data and add to the scene.
                img_uri = f"{scene_uri}/img"
                with _create_or_open_collection(
                    Collection[MultiscaleImage], img_uri, **ingest_ctx
                ) as img:
                    _maybe_set(scene, "img", img, use_relative_uri=use_relative_uri)
                    if any(
                        x is not None
                        for x in (input_hires, input_lowres, input_fullres)
                    ):
                        tissue_uri = f"{img_uri}/tissue"
                        with MultiscaleImage.create(
                            tissue_uri,
                            image_type="YXC",
                            type=pa.uint8(),
                            reference_level_shape=ref_shape,
                            axis_names=("x", "y", "c"),
                            context=ingest_ctx.get("context"),
                        ) as tissue:
                            add_metadata(tissue, ingest_ctx.get("additional_metadata"))
                            _maybe_set(
                                img,
                                "tissue",
                                tissue,
                                use_relative_uri=use_relative_uri,
                            )
                            _write_visium_images(
                                tissue,
                                scale_factors,
                                input_hires=input_hires,
                                input_lowres=input_lowres,
                                input_fullres=input_fullres,
                                use_relative_uri=use_relative_uri,
                                **ingest_ctx,
                            )
                            tissue.coordinate_space = coord_space
                            scene.register_multiscale_image(
                                "tissue",
                                IdentityTransform(("x", "y"), ("x", "y")),
                            )

                obsl_uri = f"{scene_uri}/obsl"
                with _create_or_open_collection(
                    Collection[AnySOMAObject], obsl_uri, **ingest_ctx
                ) as obsl:
                    _maybe_set(scene, "obsl", obsl, use_relative_uri=use_relative_uri)

                    loc_uri = f"{obsl_uri}/loc"
                    # Write spot data and add to the scene.
                    with _write_visium_spots(
                        loc_uri,
                        input_tissue_positions,
                        pixels_per_spot_diameter,
                        obs_df,
                        obs_id_name,
                        **ingest_ctx,
                    ) as loc:
                        _maybe_set(obsl, "loc", loc, use_relative_uri=use_relative_uri)
                        scene.register_point_cloud(
                            "loc", IdentityTransform(("x", "y"), ("x", "y"))
                        )

                varl_uri = f"{scene_uri}/varl"
                with _create_or_open_collection(
                    Collection[Collection[AnySOMAObject]], varl_uri, **ingest_ctx
                ) as varl:
                    _maybe_set(scene, "varl", varl, use_relative_uri=use_relative_uri)

        # Create the obs presence matrix.
        if write_obs_spatial_presence:
            obs_spatial_presence_uri = f"{uri}/obs_spatial_presence"
            obs_spatial_presence = _write_scene_presence_dataframe(
                obs_id, scene_name, obs_spatial_presence_uri, **ingest_ctx
            )
            _maybe_set(
                exp,
                "obs_spatial_presence",
                obs_spatial_presence,
                use_relative_uri=use_relative_uri,
            )
        if write_var_spatial_presence:
            var_spatial_presence_uri = (
                f"{uri}/ms/{measurement_name}/var_spatial_presence"
            )
            var_spatial_presence = _write_scene_presence_dataframe(
                var_id, scene_name, var_spatial_presence_uri, **ingest_ctx
            )
            meas = exp.ms[measurement_name]
            _maybe_set(
                meas,
                "var_spatial_presence",
                var_spatial_presence,
                use_relative_uri=use_relative_uri,
            )
    return uri


def _write_scene_presence_dataframe(
    joinids: pa.array,
    scene_id: str,
    df_uri: str,
    *,
    ingestion_params: IngestionParams,
    additional_metadata: "AdditionalMetadata" = None,
    platform_config: Optional["PlatformConfig"] = None,
    context: Optional["SOMATileDBContext"] = None,
) -> DataFrame:
    s = _util.get_start_stamp()
    arrow_table = None  # TODO: fixme

    try:
        soma_df = DataFrame.create(
            df_uri,
            schema=pa.schema(
                [
                    ("soma_joinid", pa.int64()),
                    ("scene_id", pa.string()),
                    ("data", pa.bool_()),
                ]
            ),
            index_column_names=("soma_joinid", "scene_id"),
            platform_config=platform_config,
            context=context,
        )
    except (AlreadyExistsError, NotCreateableError):
        if ingestion_params.error_if_already_exists:
            raise SOMAError(f"{df_uri} already exists")
        soma_df = DataFrame.open(df_uri, "w", context=context)

    if ingestion_params.write_schema_no_data:
        logging.log_io(
            f"Wrote schema {df_uri}",
            _util.format_elapsed(s, f"FINISH WRITING SCHEMA {df_uri}"),
        )
        add_metadata(soma_df, additional_metadata)
        return soma_df

    tiledb_create_options = TileDBCreateOptions.from_platform_config(platform_config)
    tiledb_write_options = TileDBWriteOptions.from_platform_config(platform_config)

    if joinids:
        nvalues = len(joinids)
        arrow_table = pa.Table.from_pydict(
            {
                "soma_joinid": joinids,
                "scene_id": nvalues * [scene_id],
                "data": nvalues * [True],
            }
        )
        _write_arrow_table(
            arrow_table, soma_df, tiledb_create_options, tiledb_write_options
        )

    logging.log_io(
        f"Wrote   {df_uri}",
        _util.format_elapsed(s, f"FINISH WRITING {df_uri}"),
    )
    return soma_df


def _write_visium_spots(
    df_uri: str,
    input_tissue_positions: Path,
    spot_diameter: float,
    obs_df: pd.DataFrame,
    id_column_name: str,
    *,
    ingestion_params: IngestionParams,
    additional_metadata: "AdditionalMetadata" = None,
    platform_config: Optional["PlatformConfig"] = None,
    context: Optional["SOMATileDBContext"] = None,
) -> PointCloud:
    """TODO: Add _write_visium_spot_dataframe docs"""
    df = (
        pd.read_csv(input_tissue_positions)
        .rename(
            columns={
                "barcode": id_column_name,
                "pxl_row_in_fullres": "y",
                "pxl_col_in_fullres": "x",
            }
        )
        .assign(spot_diameter_fullres=np.double(spot_diameter))
    )
    df = pd.merge(obs_df, df, how="inner", on=id_column_name)
    df.drop(id_column_name, axis=1, inplace=True)

    arrow_table = df_to_arrow(df)

    soma_point_cloud = PointCloud.create(
        df_uri,
        schema=arrow_table.schema,
        platform_config=platform_config,
        context=context,
    )
    # TODO: Consider moving the following to properties in the PointCloud class
    if additional_metadata is None:
        additional_metadata = {
            "soma_geometry_type": "radius",
            "soma_geometry": 0.5 * np.double(spot_diameter),  # type: ignore
        }
    else:
        additional_metadata["soma_geometry_type"] = "radius"
        additional_metadata["soma_geometry"] = 0.5 * np.double(spot_diameter)  # type: ignore
    if ingestion_params.write_schema_no_data:
        add_metadata(soma_point_cloud, additional_metadata)
        return soma_point_cloud

    tiledb_create_options = TileDBCreateOptions.from_platform_config(platform_config)
    tiledb_write_options = TileDBWriteOptions.from_platform_config(platform_config)

    if arrow_table:
        _write_arrow_table(
            arrow_table, soma_point_cloud, tiledb_create_options, tiledb_write_options
        )

    add_metadata(soma_point_cloud, additional_metadata)
    return soma_point_cloud


def _write_visium_images(
    image_pyramid: MultiscaleImage,
    scale_factors: Dict[str, Any],
    *,
    input_hires: Union[None, str, Path],
    input_lowres: Union[None, str, Path],
    input_fullres: Union[None, str, Path],
    ingestion_params: IngestionParams,
    additional_metadata: "AdditionalMetadata" = None,
    platform_config: Optional["PlatformConfig"] = None,
    context: Optional["SOMATileDBContext"] = None,
    use_relative_uri: Optional[bool] = None,
) -> None:
    # Write metadata from scale_factors file
    for key, value in scale_factors.items():
        image_pyramid.metadata[key] = value

    # Add the different levels of zoom to the image pyramid.
    for name, image_path in (
        ("fullres", input_fullres),
        ("hires", input_hires),
        ("lowres", input_lowres),
    ):
        if image_path is None:
            continue
        with Image.open(image_path) as im:
            im_data = pa.Tensor.from_numpy(np.array(im))
        im_array = image_pyramid.add_new_level(
            name, type=pa.uint8(), shape=im_data.shape
        )
        im_array.write(
            (slice(None), slice(None), slice(None)),
            im_data,
            platform_config=platform_config,
        )
