# Copyright (c) 2024 TileDB, Inc
#
# Licensed under the MIT License.

"""Experimental ingestion methods.

This module contains experimental methods to generate Spatial SOMA artifacts
start from other formats.

Do NOT merge into main.
"""

import json
import warnings
from pathlib import Path
from typing import (
    TYPE_CHECKING,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
    Union,
)

import attrs
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pacomp
import scanpy
from typing_extensions import Self

try:
    from PIL import Image
except ImportError as err:
    warnings.warn("Experimental spatial ingestor requires the `pillow` package.")
    raise err


from somacore import Axis, CoordinateSpace, IdentityTransform, ScaleTransform

from .. import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    ImageProperties,
    MultiscaleImage,
    PointCloudDataFrame,
    Scene,
    SparseNDArray,
    _util,
    logging,
)
from .._arrow_types import df_to_arrow
from .._constants import SPATIAL_DISCLAIMER
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
from ._util import _read_visium_software_version

if TYPE_CHECKING:
    from somacore.options import PlatformConfig

    from ..io._common import AdditionalMetadata
    from ..io._registration import ExperimentAmbientLabelMapping
    from ..options import SOMATileDBContext


def path_validator(instance, attribute, value: Path) -> None:  # type: ignore[no-untyped-def]
    if not value.exists():
        raise OSError(f"Path {value} does not exist")


def optional_path_converter(value: Optional[Union[str, Path]]) -> Optional[Path]:
    return None if value is None else Path(value)


def optional_path_validator(instance, attribute, x: Optional[Path]) -> None:  # type: ignore[no-untyped-def]
    if x is not None and not x.exists():
        raise OSError(f"Path {x} does not exist")


@attrs.define(kw_only=True)
class VisiumPaths:

    @classmethod
    def from_base_folder(
        cls,
        base_path: Union[str, Path],
        *,
        gene_expression: Optional[Union[str, Path]] = None,
        scale_factors: Optional[Union[str, Path]] = None,
        tissue_positions: Optional[Union[str, Path]] = None,
        fullres_image: Optional[Union[str, Path]] = None,
        hires_image: Optional[Union[str, Path]] = None,
        lowres_image: Optional[Union[str, Path]] = None,
        use_raw_counts: bool = False,
        version: Optional[Union[int, Tuple[int, int, int]]] = None,
    ) -> Self:
        """Create ingestion files from Space Ranger output directory.

        This method attempts to find the required Visium assets from an output
        directory from Space Ranger. The path for all files can be directly
        specified instead.

        The full resolution image is an input file to Space Ranger. In order to
        include it, the ``fullres_image`` argument must be specified.
        """
        base_path = Path(base_path)

        if gene_expression is None:
            gene_expression_suffix = (
                "raw_feature_bc_matrix.h5"
                if use_raw_counts
                else "filtered_feature_bc_matrix.h5"
            )

            possible_paths = list(base_path.glob(f"*{gene_expression_suffix}"))
            if len(possible_paths) == 0:
                raise OSError(
                    f"No expression matrix ending in {gene_expression_suffix} found "
                    f"in {base_path}. If the file has been renamed, it can be directly "
                    f"specified with the `gene_expression` argument."
                )
            if len(possible_paths) > 1:
                raise OSError(
                    f"Multiple files ending in {gene_expression_suffix} "
                    f"found in {base_path}. The desired file must be specified with "
                    f"the `gene_expression` argument."
                )
            gene_expression = possible_paths[0]

        return cls.from_spatial_folder(
            base_path / "spatial",
            gene_expression=gene_expression,
            scale_factors=scale_factors,
            tissue_positions=tissue_positions,
            fullres_image=fullres_image,
            hires_image=hires_image,
            lowres_image=lowres_image,
            version=version,
        )

    @classmethod
    def from_spatial_folder(
        cls,
        spatial_dir: Union[str, Path],
        *,
        gene_expression: Union[str, Path],
        scale_factors: Optional[Union[str, Path]] = None,
        tissue_positions: Optional[Union[str, Path]] = None,
        fullres_image: Optional[Union[str, Path]] = None,
        hires_image: Optional[Union[str, Path]] = None,
        lowres_image: Optional[Union[str, Path]] = None,
        version: Optional[Union[int, Tuple[int, int, int]]] = None,
    ) -> Self:
        """Create ingestion files from Space Ranger spatial output directory
        and the gene expression file.

        This method attempts to find the required Visium assets from the spatial output
        directory from Space Ranger. The path for all files can be directly specified
        instead.

        The full resolution image is an input file to Space Ranger. In order to
        include it, the `fullres_image` argument must be specified.
        """
        spatial_dir = Path(spatial_dir)

        # Attempt to read the Space Ranger version if it is not already set.
        if version is None:
            try:
                version = _read_visium_software_version(gene_expression)
            except (KeyError, ValueError):
                warnings.warn(
                    "Unable to determine Space Ranger vesion from gene expression file."
                )
        major_version = version[0] if isinstance(version, tuple) else version

        if tissue_positions is None:
            if major_version == 1:
                possible_file_names = [
                    "tissue_positions_list.csv",
                    "tissue_positions.csv",
                ]
            else:
                possible_file_names = [
                    "tissue_positions.csv",
                    "tissue_positions_list.csv",
                ]
            for possible in possible_file_names:
                tissue_positions = spatial_dir / possible
                if tissue_positions.exists():
                    break
            else:
                raise OSError(
                    f"No tissue position file found in {spatial_dir}. Tried files: "
                    f"{possible_file_names}. If the file has been renamed it can be "
                    f"directly specified using argument `tissue_positions`."
                )

        if scale_factors is None:
            scale_factors = spatial_dir / "scalefactors_json.json"

        if hires_image is None:
            hires_image = spatial_dir / "tissue_hires_image.png"
        if lowres_image is None:
            lowres_image = spatial_dir / "tissue_lowres_image.png"

        return cls(
            gene_expression=gene_expression,
            scale_factors=scale_factors,
            tissue_positions=tissue_positions,
            fullres_image=fullres_image,
            hires_image=hires_image,
            lowres_image=lowres_image,
            version=version,
        )

    gene_expression: Path = attrs.field(converter=Path, validator=path_validator)
    scale_factors: Path = attrs.field(converter=Path, validator=path_validator)
    tissue_positions: Path = attrs.field(converter=Path, validator=path_validator)
    fullres_image: Optional[Path] = attrs.field(
        converter=optional_path_converter, validator=optional_path_validator
    )
    hires_image: Optional[Path] = attrs.field(
        converter=optional_path_converter, validator=optional_path_validator
    )

    lowres_image: Optional[Path] = attrs.field(
        converter=optional_path_converter, validator=optional_path_validator
    )
    version: Optional[Union[int, Tuple[int, int, int]]] = attrs.field(default=None)

    @version.validator
    def _validate_version(  # type: ignore[no-untyped-def]
        self, attribute, value: Optional[Union[int, Tuple[int, int, int]]]
    ) -> None:
        major_version = value[0] if isinstance(value, tuple) else value
        if major_version is not None and major_version != 2:
            warnings.warn(
                f"Support for Space Ranger version {value} has not been tests."
            )

    @property
    def has_image(self) -> bool:
        return (
            self.fullres_image is not None
            or self.hires_image is not None
            or self.lowres_image is not None
        )


def from_visium(
    experiment_uri: str,
    input_path: Union[Path, VisiumPaths],
    measurement_name: str,
    scene_name: str,
    *,
    context: Optional["SOMATileDBContext"] = None,
    platform_config: Optional["PlatformConfig"] = None,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    X_layer_name: str = "data",
    raw_X_layer_name: str = "data",
    image_name: str = "tissue",
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
    X_kind: Union[Type[SparseNDArray], Type[DenseNDArray]] = SparseNDArray,
    registration_mapping: Optional["ExperimentAmbientLabelMapping"] = None,
    uns_keys: Optional[Sequence[str]] = None,
    additional_metadata: "AdditionalMetadata" = None,
    use_raw_counts: bool = False,
    write_obs_spatial_presence: bool = False,
    write_var_spatial_presence: bool = False,
) -> str:
    """Reads a 10x Visium dataset and writes it to an :class:`Experiment`.

    This function is for ingesting Visium data for prototyping and testing the
    proposed spatial design.

    WARNING: This was only tested for Space Ranger version 2 output.

    Args:
        experiment_uri: The experiment to create or update.

        input_path: A path to the base directory storing SpaceRanger output or
            a ``VisiumPaths`` object.

        measurement_name: The name of the measurement to store data in.

        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.

        platform_config: Platform-specific options used to create this array, provided in the form
          ``{\"tiledb\": {\"create\": {\"sparse_nd_array_dim_zstd_level\": 7}}}``.

        obs_id_name/var_id_name: Which AnnData ``obs`` and ``var`` columns, respectively, to use
          for append mode.

          Values of this column will be used to decide which obs/var rows in appended
          inputs are distinct from the ones already stored, for the assignment of ``soma_joinid``.  If
          this column exists in the input data, as a named index or a non-index column name, it will
          be used. If this column doesn't exist in the input data, and if the index is nameless or
          named ``index``, that index will be given this name when written to the SOMA experiment's
          ``obs`` / ``var``.

          NOTE: it is not necessary for this column to be the index-column
          name in the input AnnData objects ``obs``/``var``.

        X_layer_name: SOMA array name for the AnnData's ``X`` matrix.

        raw_X_layer_name: SOMA array name for the AnnData's ``raw/X`` matrix.

        image_name: SOMA multiscale image name for the multiscale image of the
            Space Ranger output images.

        ingest_mode: The ingestion type to perform:

            - ``write``: Writes all data, creating new layers if the SOMA already exists.
            - ``resume``: Adds data to an existing SOMA, skipping writing data
              that was previously written. Useful for continuing after a partial
              or interrupted ingestion operation.
            - ``schema_only``: Creates groups and the array schema, without
              writing any data to the array. Useful to prepare for appending
              multiple H5AD files to a single SOMA.

        X_kind: Which type of matrix is used to store dense X data from the
          H5AD file: ``DenseNDArray`` or ``SparseNDArray``.

        registration_mapping: Does not need to be supplied when ingesting a single
          H5AD/AnnData object into a single :class:`Experiment`. When multiple inputs
          are to be ingested into a single experiment, there are two steps. First:

          .. code-block:: python

              import tiledbsoma.io
              rd = tiledbsoma.io.register_h5ads(
                  experiment_uri,
                  h5ad_file_names,
                  measurement_name="RNA",
                  obs_field_name="obs_id",
                  var_field_name="var_id",
                  context=context,
              )

          Once that's been done, the data ingests per se may be done in any order,
          or in parallel, via for each ``h5ad_file_name``:

          .. code-block:: python

              tiledbsoma.io.from_h5ad(
                  experiment_uri,
                  h5ad_file_name,
                  measurement_name="RNA",
                  ingest_mode="write",
                  registration_mapping=rd,
              )

        uns_keys: Only ingest the specified top-level ``uns`` keys.
          The default is to ingest them all. Use ``uns_keys=[]``
          to not ingest any ``uns`` keys.

        additional_metadata: Optional metadata to add to the ``Experiment`` and all descendents.
          This is a coarse-grained mechanism for setting key-value pairs on all SOMA objects in an
          ``Experiment`` hierarchy. Metadata for particular objects is more commonly set like:

          .. code-block:: python

              with soma.open(uri, 'w') as exp:
                  exp.metadata.update({"aaa": "BBB"})
                  exp.obs.metadata.update({"ccc": 123})

        use_raw_counts: If ``True`` ingest the raw gene expression data, otherwise
            use the filtered gene expression data. Only used if ``input_path`` is
            not a ``VisiumPaths`` object. Defaults to ``False``.

        write_obs_spatial_presence: If ``True`` create and write data to the ``obs``
            presence matrix. Defaults to ``False``.

        write_var_spatial_presence: If ``True`` create and write data to the ``var``
            presence matrix. Defaults to ``False``.

    Returns:
        The URI of the newly created experiment.

    Lifecycle:
        Experimental
    """

    if ingest_mode != "write":
        raise NotImplementedError(
            f'the only ingest_mode currently supported is "write"; got "{ingest_mode}"'
        )

    # Disclaimer about the experimental nature of the generated experiment.
    warnings.warn(SPATIAL_DISCLAIMER, stacklevel=2)

    # Get input file locations.
    input_paths = (
        input_path
        if isinstance(input_path, VisiumPaths)
        else VisiumPaths.from_base_folder(input_path, use_raw_counts=use_raw_counts)
    )

    # Get JSON scale factors.
    with open(
        input_paths.scale_factors, mode="r", encoding="utf-8"
    ) as scale_factors_json:
        scale_factors = json.load(scale_factors_json)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = scanpy.read_10x_h5(input_paths.gene_expression)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
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

    # Set the ingestion parameters.
    ingest_ctx: IngestCtx = {
        "context": context,
        "ingestion_params": IngestionParams(ingest_mode, registration_mapping),
        "additional_metadata": additional_metadata,
    }

    # Get the spot diameters from teh scale factors file.
    pixels_per_spot_diameter = scale_factors["spot_diameter_fullres"]

    # Create a list of image paths.
    # -- Each item contains: level name, image path, and scale factors to fullres.
    image_paths: List[Tuple[str, Path, Optional[float]]] = []
    if input_paths.fullres_image is not None:
        image_paths.append(("fullres", Path(input_paths.fullres_image), None))
    if input_paths.hires_image is not None:
        scale = scale_factors["tissue_hires_scalef"]
        image_paths.append(("hires", Path(input_paths.hires_image), scale))
    if input_paths.lowres_image is not None:
        scale = scale_factors["tissue_lowres_scalef"]
        image_paths.append(("lowres", Path(input_paths.lowres_image), scale))

    # Create axes and transformations
    # mypy false positive https://github.com/python/mypy/issues/5313
    coord_space = CoordinateSpace(
        (Axis(name="x", unit="pixels"), Axis(name="y", unit="pixels"))  # type: ignore[arg-type]
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
    with Experiment.open(experiment_uri, mode="w", context=context) as exp:
        spatial_uri = _util.uri_joinpath(experiment_uri, "spatial")
        with _create_or_open_collection(
            Collection[Scene], spatial_uri, **ingest_ctx
        ) as spatial:
            _maybe_set(exp, "spatial", spatial, use_relative_uri=use_relative_uri)
            scene_uri = _util.uri_joinpath(spatial_uri, scene_name)
            with _create_or_open_scene(scene_uri, **ingest_ctx) as scene:
                _maybe_set(
                    spatial, scene_name, scene, use_relative_uri=use_relative_uri
                )
                scene.coordinate_space = coord_space

                img_uri = _util.uri_joinpath(scene_uri, "img")
                with _create_or_open_collection(
                    Collection[MultiscaleImage], img_uri, **ingest_ctx
                ) as img:
                    _maybe_set(scene, "img", img, use_relative_uri=use_relative_uri)

                    # Write image data and add to the scene.
                    if image_paths:
                        tissue_uri = _util.uri_joinpath(img_uri, image_name)
                        with _create_visium_tissue_images(
                            tissue_uri,
                            image_paths,
                            use_relative_uri=use_relative_uri,
                            **ingest_ctx,
                        ) as tissue_image:
                            _maybe_set(
                                img,
                                image_name,
                                tissue_image,
                                use_relative_uri=use_relative_uri,
                            )

                            # Add scale factors as extra metadata.
                            for key, value in scale_factors.items():
                                tissue_image.metadata[key] = value

                            scale = image_paths[0][2]
                            if scale is None:
                                scene.set_transform_to_multiscale_image(
                                    image_name,
                                    IdentityTransform(("x", "y"), ("x", "y")),
                                )
                            else:
                                level_props: ImageProperties = (
                                    tissue_image.level_properties(0)
                                )
                                updated_scales = (
                                    level_props.width
                                    / np.round(level_props.width / scale),
                                    level_props.height
                                    / np.round(level_props.height / scale),
                                )
                                scene.set_transform_to_multiscale_image(
                                    image_name,
                                    ScaleTransform(
                                        ("x", "y"), ("x", "y"), updated_scales
                                    ),
                                )
                            tissue_image.coordinate_space = CoordinateSpace(
                                (Axis(name="x", unit="pixels"), Axis(name="y", unit="pixels"))  # type: ignore[arg-type]
                            )

                obsl_uri = _util.uri_joinpath(scene_uri, "obsl")
                with _create_or_open_collection(
                    Collection[AnySOMAObject], obsl_uri, **ingest_ctx
                ) as obsl:
                    _maybe_set(scene, "obsl", obsl, use_relative_uri=use_relative_uri)

                    # Write spot data and add to the scene.
                    loc_uri = _util.uri_joinpath(obsl_uri, "loc")
                    with _write_visium_spots(
                        loc_uri,
                        input_paths.tissue_positions,
                        pixels_per_spot_diameter,
                        obs_df,
                        obs_id_name,
                        **ingest_ctx,
                    ) as loc:
                        _maybe_set(obsl, "loc", loc, use_relative_uri=use_relative_uri)
                        scene.set_transform_to_point_cloud_dataframe(
                            "loc", IdentityTransform(("x", "y"), ("x", "y"))
                        )
                        loc.coordinate_space = coord_space

                varl_uri = _util.uri_joinpath(scene_uri, "varl")
                with _create_or_open_collection(
                    Collection[Collection[AnySOMAObject]], varl_uri, **ingest_ctx
                ) as varl:
                    _maybe_set(scene, "varl", varl, use_relative_uri=use_relative_uri)

        # Create the obs presence matrix.
        if write_obs_spatial_presence:
            obs_spatial_presence_uri = _util.uri_joinpath(uri, "obs_spatial_presence")
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
            var_spatial_presence_uri = _util.uri_joinpath(
                _util.uri_joinpath(_util.uri_joinpath(uri, "ms"), measurement_name),
                "var_spatial_presence",
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
) -> PointCloudDataFrame:
    """Creates, opens, and writes data to a ``PointCloudDataFrame`` with the spot
    locations and metadata. Returns the open dataframe for writing.
    """
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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        soma_point_cloud = PointCloudDataFrame.create(
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


def _create_or_open_scene(
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: Optional["SOMATileDBContext"],
    additional_metadata: "AdditionalMetadata" = None,
) -> Scene:
    """Creates or opens a ``Scene`` and returns it open for writing."""
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            scene = Scene.create(uri, context=context)
    except (AlreadyExistsError, NotCreateableError):
        # It already exists. Are we resuming?
        if ingestion_params.error_if_already_exists:
            raise SOMAError(f"{uri} already exists")
        scene = Scene.open(uri, "w", context=context)

    add_metadata(scene, additional_metadata)
    return scene


def _create_visium_tissue_images(
    uri: str,
    image_paths: List[Tuple[str, Path, Optional[float]]],
    *,
    additional_metadata: "AdditionalMetadata" = None,
    platform_config: Optional["PlatformConfig"] = None,
    context: Optional["SOMATileDBContext"] = None,
    ingestion_params: IngestionParams,
    use_relative_uri: Optional[bool] = None,
) -> MultiscaleImage:
    """Creates, opens, and writes a ``MultiscaleImage`` with the provide
    visium resolutions levels and returns the open image for writing.
    """

    # Open the first image to get the base size.
    with Image.open(image_paths[0][1]) as im:
        im_data_numpy = np.array(im)
        ref_shape: Tuple[int, ...] = im_data_numpy.shape
        im_data = pa.Tensor.from_numpy(im_data_numpy)

    # Create the multiscale image.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        image_pyramid = MultiscaleImage.create(
            uri,
            type=pa.uint8(),
            reference_level_shape=ref_shape,
            axis_names=("y", "x", "c"),
            axis_types=("height", "width", "channel"),
            context=context,
        )

    # Add additional metadata.
    add_metadata(image_pyramid, additional_metadata)

    # Add and write the first level.
    im_array = image_pyramid.add_new_level(image_paths[0][0], shape=ref_shape)
    im_array.write(
        (slice(None), slice(None), slice(None)),
        im_data,
        platform_config=platform_config,
    )
    im_array.close()

    # Add the remaining levels.
    for name, image_path, _ in image_paths[1:]:
        with Image.open(image_path) as im:
            im_data = pa.Tensor.from_numpy(np.array(im))
        im_array = image_pyramid.add_new_level(name, shape=im_data.shape)
        im_array.write(
            (slice(None), slice(None), slice(None)),
            im_data,
            platform_config=platform_config,
        )
        im_array.close()

    return image_pyramid
