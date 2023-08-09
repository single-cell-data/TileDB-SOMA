# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Ingestion methods.

This module contains methods to generate SOMA artifacts starting from
other formats. Currently only ``.h5ad`` (`AnnData <https://anndata.readthedocs.io/>`_) is supported.
"""

import math
import time
from typing import (
    Any,
    List,
    Mapping,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
    overload,
)

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb
from anndata._core.sparse_dataset import SparseDataset
from somacore.options import PlatformConfig

from .. import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    Measurement,
    SparseNDArray,
    _factory,
    _util,
    eta,
    logging,
)
from .._arrow_types import df_to_arrow, tiledb_type_from_arrow_type
from .._collection import AnyTileDBCollection
from .._common_nd_array import NDArray
from .._constants import SOMA_JOINID
from .._exception import DoesNotExistError, SOMAError
from .._funcs import typeguard_ignore
from .._tdb_handles import RawHandle
from .._tiledb_array import TileDBArray
from .._tiledb_object import AnyTileDBObject, TileDBObject
from .._types import INGEST_MODES, IngestMode, NPNDArray, Path
from ..options import SOMATileDBContext
from ..options._soma_tiledb_context import _validate_soma_tiledb_context
from ..options._tiledb_create_options import TileDBCreateOptions
from . import conversions
from .registration import signatures

SparseMatrix = Union[sp.csr_matrix, sp.csc_matrix, SparseDataset]
DenseMatrix = Union[NPNDArray, h5py.Dataset]
Matrix = Union[DenseMatrix, SparseMatrix]
_NDArr = TypeVar("_NDArr", bound=NDArray)
_TDBO = TypeVar("_TDBO", bound=TileDBObject[RawHandle])

# ----------------------------------------------------------------
class IngestionParams:
    """Maps from user-level ingest modes to a set of implementation-level boolean flags."""

    write_schema_no_data: bool
    error_if_already_exists: bool
    skip_existing_nonempty_domain: bool

    def __init__(self, ingest_mode: str) -> None:
        if ingest_mode == "write":
            self.write_schema_no_data = False
            self.error_if_already_exists = True
            self.skip_existing_nonempty_domain = False

        elif ingest_mode == "schema_only":
            self.write_schema_no_data = True
            self.error_if_already_exists = True
            self.skip_existing_nonempty_domain = False

        elif ingest_mode == "resume":
            self.write_schema_no_data = False
            self.error_if_already_exists = False
            self.skip_existing_nonempty_domain = True

        else:
            raise SOMAError(
                f'expected ingest_mode to be one of {INGEST_MODES}; got "{ingest_mode}"'
            )


# ----------------------------------------------------------------
def from_h5ad(
    experiment_uri: str,
    input_path: Path,
    measurement_name: str,
    *,
    context: Optional[SOMATileDBContext] = None,
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
    X_kind: Union[Type[SparseNDArray], Type[DenseNDArray]] = SparseNDArray,
) -> str:
    """Reads an ``.h5ad`` file and writes it to an :class:`Experiment`.

    Measurement data is stored in a :class:`Measurement` in the experiment's
    ``ms`` field, with the key provided by ``measurement_name``. Data elements
    are available at the standard fields (``var``, ``X``, etc.). Unstructured
    data from ``uns`` is partially supported (structured arrays and non-numeric
    NDArrays are skipped), and is available at the measurement's ``uns`` key
    (i.e., at ``your_experiment.ms[measurement_name]["uns"]``).

    Args:
        experiment_uri: The experiment to create or update.

        input_path: A path to an input H5AD file.

        measurement_name: The name of the measurement to store data in.

        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.

        platform_config: Platform-specific options used to create this array, provided in the form
        ``{"tiledb": {"create": {"sparse_nd_array_dim_zstd_level": 7}}}`` nested keys.

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

    Returns:
        The URI of the newly created experiment.

    Lifecycle:
        Experimental.
    """
    if ingest_mode not in INGEST_MODES:
        raise SOMAError(
            f'expected ingest_mode to be one of {INGEST_MODES}; got "{ingest_mode}"'
        )

    if isinstance(input_path, ad.AnnData):
        raise TypeError("input path is an AnnData object -- did you want from_anndata?")

    context = _validate_soma_tiledb_context(context)

    s = _util.get_start_stamp()
    logging.log_io(None, f"START  Experiment.from_h5ad {input_path}")

    logging.log_io(None, f"START  READING {input_path}")

    anndata = ad.read_h5ad(input_path, backed="r")

    logging.log_io(None, _util.format_elapsed(s, f"FINISH READING {input_path}"))

    uri = from_anndata(
        experiment_uri,
        anndata,
        measurement_name,
        context=context,
        platform_config=platform_config,
        ingest_mode=ingest_mode,
        use_relative_uri=use_relative_uri,
        X_kind=X_kind,
    )

    logging.log_io(
        None, _util.format_elapsed(s, f"FINISH Experiment.from_h5ad {input_path} {uri}")
    )
    return uri


# ----------------------------------------------------------------
def from_anndata(
    experiment_uri: str,
    anndata: ad.AnnData,
    measurement_name: str,
    *,
    context: Optional[SOMATileDBContext] = None,
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
    X_kind: Union[Type[SparseNDArray], Type[DenseNDArray]] = SparseNDArray,
) -> str:
    """Writes an `AnnData <https://anndata.readthedocs.io/>`_ object to an :class:`Experiment`.

    Measurement data is stored in a :class:`Measurement` in the experiment's
    ``ms`` field, with the key provided by ``measurement_name``. Data elements
    are available at the standard fields (``var``, ``X``, etc.). Unstructured
    data from ``uns`` is partially supported (structured arrays and non-numeric
    NDArrays are skipped), and is available at the measurement's ``uns`` key
    (i.e., at ``your_experiment.ms[measurement_name]["uns"]``).

    Args:
        experiment_uri: The experiment to create or update.

        input_path: A path to an input H5AD file.

        measurement_name: The name of the measurement to store data in.

        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.

        platform_config:
            Platform-specific options used to create this array, provided in the form
            ``{"tiledb": {"create": {"sparse_nd_array_dim_zstd_level": 7}}}``.

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

    Returns:
        The URI of the newly created experiment.

    Lifecycle:
        Experimental.
    """
    if ingest_mode not in INGEST_MODES:
        raise SOMAError(
            f'expected ingest_mode to be one of {INGEST_MODES}; got "{ingest_mode}"'
        )

    # Map the user-level ingest mode to a set of implementation-level boolean flags
    ingestion_params = IngestionParams(ingest_mode)

    if not isinstance(anndata, ad.AnnData):
        raise TypeError(
            "Second argument is not an AnnData object -- did you want from_h5ad?"
        )

    context = _validate_soma_tiledb_context(context)

    # Without _at least_ an index, there is nothing to indicate the dimension indices.
    if anndata.obs.index.empty or anndata.var.index.empty:
        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

    s = _util.get_start_stamp()
    logging.log_io(None, "START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()

    logging.log_io(None, _util.format_elapsed(s, "FINISH DECATEGORICALIZING"))

    s = _util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {experiment_uri}")

    # Must be done first, to create the parent directory.
    experiment = _create_or_open_coll(
        Experiment, experiment_uri, ingestion_params=ingestion_params, context=context
    )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # OBS
    df_uri = _util.uri_joinpath(experiment_uri, "obs")
    with _write_dataframe(
        df_uri,
        conversions.decategoricalize_obs_or_var(anndata.obs),
        id_column_name="obs_id",
        platform_config=platform_config,
        context=context,
        ingestion_params=ingestion_params,
    ) as obs:
        _maybe_set(experiment, "obs", obs, use_relative_uri=use_relative_uri)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS

    # For local disk and S3, this is the same as experiment.ms.uri.
    # For ingest_mode="resume" on TileDB Cloud, experiment_uri will be of the form
    # tiledb://namespace/s3://bucket/path/to/exp whereas experiment.ms.uri will
    # be of the form tiledb://namespace/uuid. Only for the former is it suitable
    # to append "/ms" so that is what we do here.
    experiment_ms_uri = f"{experiment_uri}/ms"

    with _create_or_open_coll(
        Collection[Measurement],
        experiment_ms_uri,
        ingestion_params=ingestion_params,
        context=context,
    ) as ms:
        _maybe_set(experiment, "ms", ms, use_relative_uri=use_relative_uri)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # MS/meas
        measurement_uri = _util.uri_joinpath(experiment_ms_uri, measurement_name)
        with _create_or_open_coll(
            Measurement,
            measurement_uri,
            ingestion_params=ingestion_params,
            context=context,
        ) as measurement:
            _maybe_set(
                ms, measurement_name, measurement, use_relative_uri=use_relative_uri
            )

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # ms/meas/uns
            _maybe_ingest_uns(
                measurement,
                anndata.uns,
                platform_config=platform_config,
                context=context,
                ingestion_params=ingestion_params,
                use_relative_uri=use_relative_uri,
            )

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # MS/meas/VAR
            with _write_dataframe(
                _util.uri_joinpath(measurement_uri, "var"),
                conversions.decategoricalize_obs_or_var(anndata.var),
                id_column_name="var_id",
                platform_config=platform_config,
                context=context,
                ingestion_params=ingestion_params,
            ) as var:
                _maybe_set(measurement, "var", var, use_relative_uri=use_relative_uri)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # MS/meas/X/DATA

            measurement_X_uri = _util.uri_joinpath(measurement_uri, "X")
            with _create_or_open_coll(
                Collection,
                measurement_X_uri,
                ingestion_params=ingestion_params,
                context=context,
            ) as x:
                _maybe_set(measurement, "X", x, use_relative_uri=use_relative_uri)

                # Since we did ``anndata = ad.read_h5ad(path_to_h5ad, "r")`` with the "r":
                # * If we do ``anndata.X[:]`` we're loading all of a CSR/CSC/etc into memory.
                # * If we do ``anndata.X`` we're getting a pageable object which can be loaded
                #   chunkwise into memory.
                # Using the latter allows us to ingest larger .h5ad files without OOMing.

                with _create_from_matrix(
                    X_kind,
                    _util.uri_joinpath(measurement_X_uri, "data"),
                    anndata.X,
                    ingestion_params=ingestion_params,
                    platform_config=platform_config,
                    context=context,
                ) as data:
                    _maybe_set(x, "data", data, use_relative_uri=use_relative_uri)

                for layer_name, layer in anndata.layers.items():
                    with _create_from_matrix(
                        X_kind,
                        _util.uri_joinpath(measurement_X_uri, layer_name),
                        layer,
                        ingestion_params=ingestion_params,
                        platform_config=platform_config,
                        context=context,
                    ) as layer_data:
                        _maybe_set(
                            x, layer_name, layer_data, use_relative_uri=use_relative_uri
                        )

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # MS/meas/OBSM,VARM,OBSP,VARP
                if len(anndata.obsm.keys()) > 0:  # do not create an empty collection
                    obsm_uri = _util.uri_joinpath(measurement_uri, "obsm")
                    with _create_or_open_coll(
                        Collection,
                        obsm_uri,
                        ingestion_params=ingestion_params,
                        context=context,
                    ) as obsm:
                        _maybe_set(
                            measurement, "obsm", obsm, use_relative_uri=use_relative_uri
                        )
                        for key in anndata.obsm.keys():
                            with _create_from_matrix(
                                # TODO (https://github.com/single-cell-data/TileDB-SOMA/issues/1245):
                                # consider a use-dense flag at the tiledbsoma.io API
                                # DenseNDArray,
                                SparseNDArray,
                                _util.uri_joinpath(obsm.uri, key),
                                conversions.to_tiledb_supported_array_type(
                                    key, anndata.obsm[key]
                                ),
                                ingestion_params=ingestion_params,
                                platform_config=platform_config,
                                context=context,
                            ) as arr:
                                _maybe_set(
                                    obsm, key, arr, use_relative_uri=use_relative_uri
                                )
                            arr.close()
                    measurement.obsm.close()

                if len(anndata.varm.keys()) > 0:  # do not create an empty collection
                    _util.uri_joinpath(measurement_uri, "varm")
                    with _create_or_open_coll(
                        Collection,
                        _util.uri_joinpath(measurement.uri, "varm"),
                        ingestion_params=ingestion_params,
                        context=context,
                    ) as varm:
                        _maybe_set(
                            measurement, "varm", varm, use_relative_uri=use_relative_uri
                        )
                        for key in anndata.varm.keys():
                            with _create_from_matrix(
                                # TODO (https://github.com/single-cell-data/TileDB-SOMA/issues/1245):
                                # consider a use-dense flag at the tiledbsoma.io API
                                # DenseNDArray,
                                SparseNDArray,
                                _util.uri_joinpath(varm.uri, key),
                                conversions.to_tiledb_supported_array_type(
                                    key, anndata.varm[key]
                                ),
                                ingestion_params=ingestion_params,
                                platform_config=platform_config,
                                context=context,
                            ) as darr:
                                _maybe_set(
                                    varm,
                                    key,
                                    darr,
                                    use_relative_uri=use_relative_uri,
                                )

                if len(anndata.obsp.keys()) > 0:  # do not create an empty collection
                    _util.uri_joinpath(measurement_uri, "obsp")
                    with _create_or_open_coll(
                        Collection,
                        _util.uri_joinpath(measurement.uri, "obsp"),
                        ingestion_params=ingestion_params,
                        context=context,
                    ) as obsp:
                        _maybe_set(
                            measurement, "obsp", obsp, use_relative_uri=use_relative_uri
                        )
                        for key in anndata.obsp.keys():
                            with _create_from_matrix(
                                SparseNDArray,
                                _util.uri_joinpath(obsp.uri, key),
                                conversions.to_tiledb_supported_array_type(
                                    key, anndata.obsp[key]
                                ),
                                ingestion_params=ingestion_params,
                                platform_config=platform_config,
                                context=context,
                            ) as sarr:
                                _maybe_set(
                                    obsp,
                                    key,
                                    sarr,
                                    use_relative_uri=use_relative_uri,
                                )

                if len(anndata.varp.keys()) > 0:  # do not create an empty collection
                    _util.uri_joinpath(measurement_uri, "obsp")
                    with _create_or_open_coll(
                        Collection,
                        _util.uri_joinpath(measurement.uri, "varp"),
                        ingestion_params=ingestion_params,
                        context=context,
                    ) as varp:
                        _maybe_set(
                            measurement, "varp", varp, use_relative_uri=use_relative_uri
                        )
                        for key in anndata.varp.keys():
                            with _create_from_matrix(
                                SparseNDArray,
                                _util.uri_joinpath(varp.uri, key),
                                conversions.to_tiledb_supported_array_type(
                                    key, anndata.varp[key]
                                ),
                                ingestion_params=ingestion_params,
                                platform_config=platform_config,
                                context=context,
                            ) as sarr:
                                _maybe_set(
                                    varp,
                                    key,
                                    sarr,
                                    use_relative_uri=use_relative_uri,
                                )

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # MS/RAW
                if anndata.raw is not None:
                    raw_uri = _util.uri_joinpath(experiment_ms_uri, "raw")
                    with _create_or_open_coll(
                        Measurement,
                        raw_uri,
                        ingestion_params=ingestion_params,
                        context=context,
                    ) as raw_measurement:
                        _maybe_set(
                            ms,
                            "raw",
                            raw_measurement,
                            use_relative_uri=use_relative_uri,
                        )

                        with _write_dataframe(
                            _util.uri_joinpath(raw_uri, "var"),
                            conversions.decategoricalize_obs_or_var(anndata.raw.var),
                            id_column_name="var_id",
                            ingestion_params=ingestion_params,
                            platform_config=platform_config,
                            context=context,
                        ) as var:
                            _maybe_set(
                                raw_measurement,
                                "var",
                                var,
                                use_relative_uri=use_relative_uri,
                            )

                        raw_X_uri = _util.uri_joinpath(raw_uri, "X")
                        with _create_or_open_coll(
                            Collection,
                            raw_X_uri,
                            ingestion_params=ingestion_params,
                            context=context,
                        ) as rm_x:
                            _maybe_set(
                                raw_measurement,
                                "X",
                                rm_x,
                                use_relative_uri=use_relative_uri,
                            )

                            with _create_from_matrix(
                                SparseNDArray,
                                _util.uri_joinpath(raw_X_uri, "data"),
                                anndata.raw.X,
                                ingestion_params=ingestion_params,
                                platform_config=platform_config,
                                context=context,
                            ) as rm_x_data:
                                _maybe_set(
                                    rm_x,
                                    "data",
                                    rm_x_data,
                                    use_relative_uri=use_relative_uri,
                                )

    experiment.close()

    logging.log_io(
        f"Wrote   {experiment.uri}",
        _util.format_elapsed(s, f"FINISH WRITING {experiment.uri}"),
    )
    return experiment.uri


def _maybe_set(
    coll: AnyTileDBCollection,
    key: str,
    value: AnyTileDBObject,
    *,
    use_relative_uri: Optional[bool],
) -> None:
    if coll.closed or coll.mode != "w":
        raise SOMAError(f"Collection must be open for write: {coll.uri}")
    try:
        coll.set(key, value, use_relative_uri=use_relative_uri)
    except SOMAError:
        # This is already a member of the collection.
        pass


@overload
def _create_or_open_coll(
    cls: Type[Experiment],
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: Optional[SOMATileDBContext],
) -> Experiment:
    ...


@overload
def _create_or_open_coll(
    cls: Type[Measurement],
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: Optional[SOMATileDBContext],
) -> Measurement:
    ...


@overload
def _create_or_open_coll(
    cls: Type[Collection[_TDBO]],
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: Optional[SOMATileDBContext],
) -> Collection[_TDBO]:
    ...


@typeguard_ignore
def _create_or_open_coll(
    cls: Type[Any],
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: Optional[SOMATileDBContext],
) -> Any:
    try:
        thing = cls.open(uri, "w", context=context)
    except DoesNotExistError:
        pass  # This is always OK; make a new one.
    else:
        # It already exists. Are we resuming?
        if ingestion_params.error_if_already_exists:
            raise SOMAError(f"{uri} already exists")
        return thing

    return cls.create(uri, context=context)


def _write_dataframe(
    df_uri: str,
    df: pd.DataFrame,
    id_column_name: Optional[str],
    *,
    ingestion_params: IngestionParams,
    platform_config: Optional[PlatformConfig] = None,
    context: Optional[SOMATileDBContext] = None,
) -> DataFrame:
    df[SOMA_JOINID] = np.arange(len(df), dtype=np.int64)

    df.reset_index(inplace=True)
    if id_column_name is not None:
        df.rename(columns={"index": id_column_name}, inplace=True)
    df.set_index(SOMA_JOINID, inplace=True)

    return _write_dataframe_impl(
        df,
        df_uri,
        ingestion_params=ingestion_params,
        platform_config=platform_config,
        context=context,
    )


def _write_dataframe_impl(
    df: pd.DataFrame,
    df_uri: str,
    *,
    ingestion_params: IngestionParams,
    platform_config: Optional[PlatformConfig] = None,
    context: Optional[SOMATileDBContext] = None,
) -> DataFrame:
    s = _util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {df_uri}")

    arrow_table = df_to_arrow(df)

    try:
        soma_df = _factory.open(df_uri, "w", soma_type=DataFrame, context=context)
    except DoesNotExistError:
        soma_df = DataFrame.create(
            df_uri,
            schema=arrow_table.schema,
            platform_config=platform_config,
            context=context,
        )
    else:
        if ingestion_params.skip_existing_nonempty_domain:
            storage_ned = _read_nonempty_domain(soma_df)
            dim_range = ((int(df.index.min()), int(df.index.max())),)
            if _chunk_is_contained_in(dim_range, storage_ned):
                logging.log_io(
                    f"Skipped {soma_df.uri}",
                    _util.format_elapsed(s, f"SKIPPED {soma_df.uri}"),
                )
                return soma_df
        else:
            raise SOMAError(f"{soma_df.uri} already exists")

    if ingestion_params.write_schema_no_data:
        logging.log_io(
            f"Wrote schema {soma_df.uri}",
            _util.format_elapsed(s, f"FINISH WRITING SCHEMA {soma_df.uri}"),
        )
        return soma_df

    soma_df.write(arrow_table)
    logging.log_io(
        f"Wrote   {soma_df.uri}",
        _util.format_elapsed(s, f"FINISH WRITING {soma_df.uri}"),
    )
    return soma_df


@typeguard_ignore
def create_from_matrix(
    cls: Type[_NDArr],
    uri: str,
    matrix: Union[Matrix, h5py.Dataset],
    platform_config: Optional[PlatformConfig] = None,
    ingest_mode: IngestMode = "write",
    context: Optional[SOMATileDBContext] = None,
) -> _NDArr:
    """
    Create and populate the ``soma_matrix`` from the contents of ``matrix``.

    Lifecycle:
        Experimental.
    """
    return _create_from_matrix(
        cls,
        uri,
        matrix,
        ingestion_params=IngestionParams(ingest_mode),
        platform_config=platform_config,
        context=context,
    )


@typeguard_ignore
def _create_from_matrix(
    cls: Type[_NDArr],
    uri: str,
    matrix: Union[Matrix, h5py.Dataset],
    *,
    ingestion_params: IngestionParams,
    platform_config: Optional[PlatformConfig] = None,
    context: Optional[SOMATileDBContext] = None,
) -> _NDArr:
    """
    Internal helper for user-facing ``create_from_matrix``.
    """
    # SparseDataset has no ndim but it has a shape
    if len(matrix.shape) != 2:
        raise ValueError(f"expected matrix.shape == 2; got {matrix.shape}")

    s = _util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {uri}")

    try:
        soma_ndarray = cls.open(
            uri, "w", platform_config=platform_config, context=context
        )
    except DoesNotExistError:
        # A SparseNDArray must be appendable in soma.io.
        shape = [None for _ in matrix.shape] if cls.is_sparse else matrix.shape
        soma_ndarray = cls.create(
            uri,
            type=pa.from_numpy_dtype(matrix.dtype),
            shape=shape,
            platform_config=platform_config,
            context=context,
        )
    else:
        if ingestion_params.error_if_already_exists:
            raise SOMAError(f"{soma_ndarray.uri} already exists")

    if ingestion_params.write_schema_no_data:
        logging.log_io(
            f"Wrote schema {soma_ndarray.uri}",
            _util.format_elapsed(s, f"FINISH WRITING SCHEMA {soma_ndarray.uri}"),
        )
        return soma_ndarray

    logging.log_io(
        f"Writing {soma_ndarray.uri}",
        _util.format_elapsed(s, f"START  WRITING {soma_ndarray.uri}"),
    )

    if isinstance(soma_ndarray, DenseNDArray):
        _write_matrix_to_denseNDArray(
            soma_ndarray,
            matrix,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(
                platform_config
            ),
            context=context,
            ingestion_params=ingestion_params,
        )
    elif isinstance(soma_ndarray, SparseNDArray):  # SOMASparseNDArray
        _write_matrix_to_sparseNDArray(
            soma_ndarray,
            matrix,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(
                platform_config
            ),
            context=context,
            ingestion_params=ingestion_params,
        )
    else:
        raise TypeError(f"unknown array type {type(soma_ndarray)}")

    logging.log_io(
        f"Wrote   {soma_ndarray.uri}",
        _util.format_elapsed(s, f"FINISH WRITING {soma_ndarray.uri}"),
    )
    return soma_ndarray


def update_obs(
    exp: Experiment,
    new_data: pd.DataFrame,
    *,
    default_index_name: str = "obs_id",
    platform_config: Optional[PlatformConfig] = None,
) -> None:
    """
    TO DO: WRITE ME

    Lifecycle:
        Experimental.
    """
    _update_dataframe(
        exp.obs,
        new_data,
        default_index_name=default_index_name,
        platform_config=platform_config,
    )


def update_var(
    exp: Experiment,
    measurement_name: str,
    new_data: pd.DataFrame,
    *,
    default_index_name: str = "var_id",
    platform_config: Optional[PlatformConfig] = None,
    # XXX filters tiledb_create_options: TileDBCreateOptions,
) -> None:
    """
    TO DO: WRITE ME

    Lifecycle:
        Experimental.
    """
    # TODO: check measurement exists
    _update_dataframe(
        exp.ms[measurement_name].var,
        new_data,
        default_index_name=default_index_name,
        platform_config=platform_config,
    )


def _update_dataframe(
    sdf: DataFrame,
    new_data: pd.DataFrame,
    *,
    default_index_name: str,
    platform_config: Optional[PlatformConfig],
) -> None:
    """
    TO DO: WRITE ME

    Lifecycle:
        Experimental.
    """
    if sdf.closed or sdf.mode != "w":
        raise SOMAError(f"DataFrame must be open for write: {sdf.uri}")
    old_sig = signatures._string_dict_from_arrow_schema(sdf.schema)
    new_sig = signatures._string_dict_from_pandas_dataframe(
        new_data, default_index_name
    )

    old_keys = set(old_sig.keys())
    new_keys = set(new_sig.keys())

    drop_keys = old_keys.difference(new_keys)
    add_keys = new_keys.difference(old_keys)
    common_keys = old_keys.intersection(new_keys)

    tiledb_create_options = TileDBCreateOptions.from_platform_config(platform_config)

    msgs = []
    for key in common_keys:
        old_type = old_sig[key]
        new_type = new_sig[key]
        if old_type != new_type:
            msgs.append(f"{key} type {old_type} != {new_type}")
    if msgs:
        msg = ", ".join(msgs)
        raise ValueError(f"unsupported type updates: {msg}")

    se = tiledb.ArraySchemaEvolution(sdf.context.tiledb_ctx)
    for drop_key in drop_keys:
        se.drop_attribute(drop_key)

    arrow_table = df_to_arrow(new_data)
    arrow_schema = arrow_table.schema.remove_metadata()

    for add_key in add_keys:
        # Don't directly use the new dataframe's dtypes. Go through the
        # to-Arrow-schema logic, and back, as this recapitulates the original
        # schema-creation logic.
        atype = arrow_schema.field(add_key).type
        dtype = tiledb_type_from_arrow_type(atype)
        filters = tiledb_create_options.attr_filters_tiledb(add_key, ["ZstdFilter"])
        se.add_attribute(
            tiledb.Attr(
                name=add_key,
                dtype=dtype,
                filters=filters,
            )
        )

    se.array_evolve(uri=sdf.uri)


def add_X_layer(
    exp: Experiment,
    measurement_name: str,
    X_layer_name: str,
    # E.g. a scipy.csr_matrix from scanpy analysis:
    X_layer_data: Union[Matrix, h5py.Dataset],
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
) -> None:
    """This is useful for adding X data, for example from
    `Scanpy <https://scanpy.readthedocs.io/>`_'s ``scanpy.pp.normalize_total``,
    ``scanpy.pp.log1p``, etc.

    Use ``ingest_mode="resume"`` to not error out if the schema already exists.

    Lifecycle:
        Experimental.
    """
    if exp.closed or exp.mode != "w":
        raise SOMAError(f"Experiment must be open for write: {exp.uri}")
    add_matrix_to_collection(
        exp,
        measurement_name,
        "X",
        X_layer_name,
        X_layer_data,
        use_relative_uri=use_relative_uri,
    )


def add_matrix_to_collection(
    exp: Experiment,
    measurement_name: str,
    collection_name: str,
    matrix_name: str,
    # E.g. a scipy.csr_matrix from scanpy analysis:
    matrix_data: Union[Matrix, h5py.Dataset],
    ingest_mode: IngestMode = "write",
    use_relative_uri: Optional[bool] = None,
    context: Optional[SOMATileDBContext] = None,
) -> None:
    """This is useful for adding X/obsp/varm/etc data, for example from
    `Scanpy <https://scanpy.readthedocs.io/>`_'s ``scanpy.pp.normalize_total``,
    ``scanpy.pp.log1p``, etc.

    Use ``ingest_mode="resume"`` to not error out if the schema already exists.

    Lifecycle:
        Experimental.
    """

    ingestion_params = IngestionParams(ingest_mode)

    # For local disk and S3, creation and storage URIs are identical.  For
    # cloud, creation URIs look like tiledb://namespace/s3://bucket/path/to/obj
    # whereas storage URIs (for the same object) look like
    # tiledb://namespace/uuid.  When the caller passes a creation URI (which
    # they must) via exp.uri, we need to follow that.
    extend_creation_uri = exp.uri.startswith("tiledb://")

    with exp.ms[measurement_name] as meas:
        if extend_creation_uri:
            coll_uri = f"{exp.uri}/ms/{measurement_name}/{collection_name}"
        else:
            coll_uri = f"{meas.uri}/{collection_name}"

        if collection_name in meas:
            coll = cast(Collection[RawHandle], meas[collection_name])
        else:
            coll = _create_or_open_coll(
                Collection,
                coll_uri,
                ingestion_params=ingestion_params,
                context=context,
            )
            _maybe_set(meas, collection_name, coll, use_relative_uri=use_relative_uri)
        with coll:
            matrix_uri = f"{coll_uri}/{matrix_name}"

            with _create_from_matrix(
                SparseNDArray,
                matrix_uri,
                matrix_data,
                ingestion_params=ingestion_params,
                context=context,
            ) as sparse_nd_array:
                _maybe_set(
                    coll,
                    matrix_name,
                    sparse_nd_array,
                    use_relative_uri=use_relative_uri,
                )


def _write_matrix_to_denseNDArray(
    soma_ndarray: DenseNDArray,
    matrix: Union[Matrix, h5py.Dataset],
    tiledb_create_options: TileDBCreateOptions,
    context: Optional[SOMATileDBContext],
    ingestion_params: IngestionParams,
) -> None:
    """Write a matrix to an empty DenseNDArray"""

    # There is a chunk-by-chunk already-done check for resume mode, below.
    # This full-matrix-level check here might seem redundant, but in fact it's important:
    # * By checking input bounds against storage NED here, we can see if the entire matrix
    #   was already ingested and avoid even loading chunks;
    # * By checking chunkwise we can catch the case where a matrix was already *partly*
    #   ingested.
    # * Of course, this also helps us catch already-completed writes in the non-chunked case.
    # TODO: make sure we're not using an old timestamp for this
    storage_ned = None
    if ingestion_params.skip_existing_nonempty_domain:
        # This lets us check for already-ingested chunks, when in resume-ingest mode.
        storage_ned = _read_nonempty_domain(soma_ndarray)
        matrix_bounds = [
            (0, int(n - 1)) for n in matrix.shape
        ]  # Cast for lint in case np.int64
        logging.log_io(
            None,
            f"Input bounds {tuple(matrix_bounds)} storage non-empty domain {storage_ned}",
        )
        if _chunk_is_contained_in(matrix_bounds, storage_ned):
            logging.log_io(
                f"Skipped {soma_ndarray.uri}", f"SKIPPED WRITING {soma_ndarray.uri}"
            )
            return

    # Write all at once?
    if not tiledb_create_options.write_X_chunked:
        if not isinstance(matrix, np.ndarray):
            matrix = matrix.toarray()
        soma_ndarray.write((slice(None),), pa.Tensor.from_numpy(matrix))
        return

    # OR, write in chunks
    eta_tracker = eta.Tracker()
    nrow, ncol = matrix.shape
    i = 0
    # Number of rows to chunk by. Dense writes, so this is a constant.
    chunk_size = int(math.ceil(tiledb_create_options.goal_chunk_nnz / ncol))
    while i < nrow:
        t1 = time.time()
        i2 = i + chunk_size

        # Print doubly-inclusive lo..hi like 0..17 and 18..31.
        chunk_percent = min(100, 100 * (i2 - 1) / nrow)
        logging.log_io(
            None,
            "START  chunk rows %d..%d of %d (%.3f%%)"
            % (i, i2 - 1, nrow, chunk_percent),
        )

        chunk = matrix[i:i2, :]

        if ingestion_params.skip_existing_nonempty_domain and storage_ned is not None:
            chunk_bounds = matrix_bounds
            chunk_bounds[0] = (
                int(i),
                int(i2 - 1),
            )  # Cast for lint in case np.int64
            if _chunk_is_contained_in_axis(chunk_bounds, storage_ned, 0):
                # Print doubly inclusive lo..hi like 0..17 and 18..31.
                logging.log_io(
                    "... %7.3f%% done" % chunk_percent,
                    "SKIP   chunk rows %d..%d of %d (%.3f%%)"
                    % (i, i2 - 1, nrow, chunk_percent),
                )
                i = i2
                continue

        if isinstance(chunk, np.ndarray):
            tensor = pa.Tensor.from_numpy(chunk)
        else:
            tensor = pa.Tensor.from_numpy(chunk.toarray())
        soma_ndarray.write((slice(i, i2), slice(None)), tensor)

        t2 = time.time()
        chunk_seconds = t2 - t1
        eta_seconds = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

        if chunk_percent < 100:
            logging.log_io(
                "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                "FINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                % (chunk_seconds, chunk_percent, eta_seconds),
            )

        i = i2

    return


def _read_nonempty_domain(arr: TileDBArray) -> Any:
    try:
        return arr._handle.reader.nonempty_domain()
    except SOMAError:
        # This means that we're open in write-only mode.
        # Reopen the array in read mode.
        pass

    cls = type(arr)
    with cls.open(arr.uri, "r", platform_config=None, context=arr.context) as readarr:
        return readarr._handle.reader.nonempty_domain()


def _find_sparse_chunk_size(
    matrix: SparseMatrix, start_index: int, axis: int, goal_chunk_nnz: int
) -> int:
    """Given a sparse matrix and a start index, return a step size, on the stride axis,
    which will achieve the cumulative nnz desired.

    Args:
        matrix:
            The input scipy.sparse matrix.
        start_index:
            The index at which to start a chunk.
        axis:
            The stride axis, across which to find a chunk.
        goal_chunk_nnz:
            Desired number of non-zero array entries for the chunk.
    """
    chunk_size = 1
    sum_nnz = 0
    coords: List[Union[slice, int]] = [slice(None), slice(None)]

    # Empirically we find:
    # * If the input matrix is sp.csr_matrix or sp.csc_matrix then getting all these nnz values is
    #   quick.
    # * If the input matrix is anndata._core.sparse_dataset.SparseDataset -- which happens with
    #   out-of-core anndata reads -- then getting all these nnz values is prohibitively expensive.
    # * It turns out that getting a sample is quite sufficient. We do this regardless of whether
    #   the matrix is anndata._core.sparse_dataset.SparseDataset or not.
    # * The max_rows is manually defined after running experiments with 60GB .h5ad files.
    count = 0
    max_rows = 100

    for index in range(start_index, matrix.shape[axis]):
        count += 1
        coords[axis] = index
        sum_nnz += matrix[tuple(coords)].nnz
        if sum_nnz > goal_chunk_nnz:
            break
        if count > max_rows:
            break
        chunk_size += 1

    if sum_nnz == 0:  # completely empty sparse array (corner case)
        return 1

    if sum_nnz > goal_chunk_nnz:
        return chunk_size

    # Solve the equation:
    #
    # sum_nnz              count
    # -------          =  -------
    # goal_chunk_nnz       result
    chunk_size = int(count * goal_chunk_nnz / sum_nnz)
    if chunk_size < 1:
        chunk_size = 1
    return chunk_size


def _write_matrix_to_sparseNDArray(
    soma_ndarray: SparseNDArray,
    matrix: Matrix,
    tiledb_create_options: TileDBCreateOptions,
    context: Optional[SOMATileDBContext],
    ingestion_params: IngestionParams,
) -> None:
    """Write a matrix to an empty DenseNDArray"""

    def _coo_to_table(mat_coo: sp.coo_matrix, axis: int = 0, base: int = 0) -> pa.Table:
        pydict = {
            "soma_data": mat_coo.data,
            "soma_dim_0": mat_coo.row + base if base > 0 and axis == 0 else mat_coo.row,
            "soma_dim_1": mat_coo.col + base if base > 0 and axis == 1 else mat_coo.col,
        }
        return pa.Table.from_pydict(pydict)

    # There is a chunk-by-chunk already-done check for resume mode, below.
    # This full-matrix-level check here might seem redundant, but in fact it's important:
    # * By checking input bounds against storage NED here, we can see if the entire matrix
    #   was already ingested and avoid even loading chunks;
    # * By checking chunkwise we can catch the case where a matrix was already *partly*
    #   ingested.
    # * Of course, this also helps us catch already-completed writes in the non-chunked case.
    # TODO: make sure we're not using an old timestamp for this
    storage_ned = None
    if ingestion_params.skip_existing_nonempty_domain:
        # This lets us check for already-ingested chunks, when in resume-ingest mode.
        # THIS IS A HACK AND ONLY WORKS BECAUSE WE ARE DOING THIS BEFORE ALL WRITES.
        storage_ned = _read_nonempty_domain(soma_ndarray)
        matrix_bounds = [
            (0, int(n - 1)) for n in matrix.shape
        ]  # Cast for lint in case np.int64
        logging.log_io(
            None,
            f"Input bounds {tuple(matrix_bounds)} storage non-empty domain {storage_ned}",
        )
        if _chunk_is_contained_in(matrix_bounds, storage_ned):
            logging.log_io(
                f"Skipped {soma_ndarray.uri}", f"SKIPPED WRITING {soma_ndarray.uri}"
            )
            return

    # Write all at once?
    if not tiledb_create_options.write_X_chunked:
        soma_ndarray.write(_coo_to_table(sp.coo_matrix(matrix)))
        return

    # Or, write in chunks, striding across the most efficient slice axis

    stride_axis = 0
    if sp.isspmatrix_csc(matrix):
        # E.g. if we used anndata.X[:]
        stride_axis = 1
    if isinstance(matrix, SparseDataset) and matrix.format_str == "csc":
        # E.g. if we used anndata.X without the [:]
        stride_axis = 1

    dim_max_size = matrix.shape[stride_axis]

    eta_tracker = eta.Tracker()
    goal_chunk_nnz = tiledb_create_options.goal_chunk_nnz

    coords = [slice(None), slice(None)]
    i = 0
    while i < dim_max_size:
        t1 = time.time()

        # Chunk size on the stride axis
        if isinstance(matrix, (np.ndarray, h5py.Dataset)):
            # These are dense, being ingested as sparse.
            # Example:
            # * goal_chunk_nnz = 100M
            # * An obsm element has shape (32458, 2)
            # * stride_axis = 0 (ingest row-wise)
            # * We want ceiling of 100M / 2 to obtain chunk_size = 50M rows
            # * This attains goal_chunk_nnz = 100_000_000 since we have 50M rows
            #   with 2 elements each
            # * Result: we divide by the shape, slotted by the non-stride axis
            non_stride_axis = 1 - stride_axis
            chunk_size = int(math.ceil(goal_chunk_nnz / matrix.shape[non_stride_axis]))
        else:
            chunk_size = _find_sparse_chunk_size(  # type: ignore [unreachable]
                matrix, i, stride_axis, goal_chunk_nnz
            )

        i2 = i + chunk_size

        coords[stride_axis] = slice(i, i2)
        chunk_coo = sp.coo_matrix(matrix[tuple(coords)])

        chunk_percent = min(100, 100 * (i2 - 1) / dim_max_size)

        if ingestion_params.skip_existing_nonempty_domain and storage_ned is not None:
            chunk_bounds = matrix_bounds
            chunk_bounds[stride_axis] = (
                int(i),
                int(i2 - 1),
            )  # Cast for lint in case np.int64
            if _chunk_is_contained_in_axis(chunk_bounds, storage_ned, stride_axis):
                # Print doubly inclusive lo..hi like 0..17 and 18..31.
                logging.log_io(
                    "... %7.3f%% done" % chunk_percent,
                    "SKIP   chunk rows %d..%d of %d (%.3f%%), nnz=%d"
                    % (i, i2 - 1, dim_max_size, chunk_percent, chunk_coo.nnz),
                )
                i = i2
                continue

        # Print doubly inclusive lo..hi like 0..17 and 18..31.
        logging.log_io(
            None,
            "START  chunk rows %d..%d of %d (%.3f%%), nnz=%d"
            % (i, i2 - 1, dim_max_size, chunk_percent, chunk_coo.nnz),
        )

        soma_ndarray.write(_coo_to_table(chunk_coo, stride_axis, i))

        t2 = time.time()
        chunk_seconds = t2 - t1
        eta_seconds = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

        if chunk_percent < 100:
            logging.log_io(
                "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                "FINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                % (chunk_seconds, chunk_percent, eta_seconds),
            )

        i = i2


def _chunk_is_contained_in(
    chunk_bounds: Sequence[Tuple[int, int]],
    storage_nonempty_domain: Optional[Sequence[Tuple[Optional[int], Optional[int]]]],
) -> bool:
    """
    Determines if a dim range is included within the array's non-empty domain.  Ranges are inclusive
    on both endpoints.  This is a helper for resume-ingest mode.

    We say "bounds" not "MBR" with the "M" for minimum: a sparse matrix might not _have_ any
    elements for some initial/final rows or columns. Suppose an input array has shape 100 x 200, so
    bounds ``((0, 99), (0, 199))`` -- and also suppose there are no matrix elements for column 1.
    Also suppose the matrix has already been written to TileDB-SOMA storage. The TileDB non-empty
    domain _is_ tight -- it'd say ``((0, 99), (3, 197))`` for example.  When we come back for a
    resume-mode ingest, we'd see the input bounds aren't contained within the storage non-empty
    domain, and erroneously declare that the data need to be rewritten.

    This is why we take the stride axis as an argument. In resume mode, it's our contract with the
    user that they declare they are retrying the exact same input file -- and we do our best to
    fulfill their ask by checking the dimension being strided on.
    """
    if storage_nonempty_domain is None:
        return False

    if len(chunk_bounds) != len(storage_nonempty_domain):
        raise SOMAError(
            f"internal error: ingest data ndim {len(chunk_bounds)} != storage ndim {len(storage_nonempty_domain)}"
        )
    for i in range(len(chunk_bounds)):
        if not _chunk_is_contained_in_axis(chunk_bounds, storage_nonempty_domain, i):
            return False
    return True


def _chunk_is_contained_in_axis(
    chunk_bounds: Sequence[Tuple[int, int]],
    storage_nonempty_domain: Sequence[Tuple[Optional[int], Optional[int]]],
    stride_axis: int,
) -> bool:
    """Helper function for ``_chunk_is_contained_in``."""
    storage_lo, storage_hi = storage_nonempty_domain[stride_axis]
    if storage_lo is None or storage_hi is None:
        # E.g. an array has had its schema created but no data written yet
        return False

    chunk_lo, chunk_hi = chunk_bounds[stride_axis]
    if chunk_lo < storage_lo or chunk_lo > storage_hi:
        return False
    if chunk_hi < storage_lo or chunk_hi > storage_hi:
        return False

    return True


def _maybe_ingest_uns(
    m: Measurement,
    uns: Mapping[str, object],
    *,
    platform_config: Optional[PlatformConfig],
    context: Optional[SOMATileDBContext],
    ingestion_params: IngestionParams,
    use_relative_uri: Optional[bool],
) -> None:
    # Don't try to ingest an empty uns.
    if not uns:
        return
    _ingest_uns_dict(
        m,
        "uns",
        uns,
        platform_config=platform_config,
        context=context,
        ingestion_params=ingestion_params,
        use_relative_uri=use_relative_uri,
    )


def _ingest_uns_dict(
    parent: AnyTileDBCollection,
    parent_key: str,
    dct: Mapping[str, object],
    *,
    platform_config: Optional[PlatformConfig],
    context: Optional[SOMATileDBContext],
    ingestion_params: IngestionParams,
    use_relative_uri: Optional[bool],
) -> None:

    with _create_or_open_coll(
        Collection,
        _util.uri_joinpath(parent.uri, parent_key),
        ingestion_params=ingestion_params,
        context=context,
    ) as coll:
        _maybe_set(parent, parent_key, coll, use_relative_uri=use_relative_uri)
        coll.metadata["soma_tiledbsoma_type"] = "uns"
        for key, value in dct.items():
            _ingest_uns_node(
                coll,
                key,
                value,
                platform_config=platform_config,
                context=context,
                ingestion_params=ingestion_params,
                use_relative_uri=use_relative_uri,
            )

    msg = f"Wrote   {coll.uri} (uns collection)"
    logging.log_io(msg, msg)


def _ingest_uns_node(
    coll: Any,
    key: Any,
    value: Any,
    *,
    platform_config: Optional[PlatformConfig],
    context: Optional[SOMATileDBContext],
    ingestion_params: IngestionParams,
    use_relative_uri: Optional[bool],
) -> None:

    if isinstance(value, np.generic):
        # This is some kind of numpy scalar value. Metadata entries
        # only accept native Python types, so unwrap it.
        value = value.item()

    if isinstance(value, (int, float, str)):
        # Primitives get set on the metadata.
        coll.metadata[key] = value
        return

    if isinstance(value, Mapping):
        # Mappings are represented as sub-dictionaries.
        _ingest_uns_dict(
            coll,
            key,
            value,
            platform_config=platform_config,
            context=context,
            ingestion_params=ingestion_params,
            use_relative_uri=use_relative_uri,
        )
        return

    if isinstance(value, pd.DataFrame):
        with _write_dataframe(
            _util.uri_joinpath(coll.uri, key),
            value,
            None,
            platform_config=platform_config,
            context=context,
            ingestion_params=ingestion_params,
        ) as df:
            _maybe_set(coll, key, df, use_relative_uri=use_relative_uri)
        return

    if isinstance(value, list) or "numpy" in str(type(value)):
        value = np.asarray(value)
    if isinstance(value, np.ndarray):
        if value.dtype.names is not None:
            msg = f"Skipped {coll.uri}[{key!r}]" " (uns): unsupported structured array"
            # This is a structured array, which we do not support.
            logging.log_io(msg, msg)
            return

        if value.dtype.char in ("U", "O"):
            # In the wild it's quite common to see arrays of strings in uns data.
            # Frequent example: uns["louvain_colors"].
            _ingest_uns_string_array(
                coll,
                key,
                value,
                platform_config,
                context=context,
                use_relative_uri=use_relative_uri,
                ingestion_params=ingestion_params,
            )
        else:
            _ingest_uns_ndarray(
                coll,
                key,
                value,
                platform_config,
                context=context,
                use_relative_uri=use_relative_uri,
                ingestion_params=ingestion_params,
            )
        return

    msg = (
        f"Skipped {coll.uri}[{key!r}]" f" (uns object): unrecognized type {type(value)}"
    )
    logging.log_io(msg, msg)


def _ingest_uns_string_array(
    coll: AnyTileDBCollection,
    key: str,
    value: NPNDArray,
    platform_config: Optional[PlatformConfig],
    context: Optional[SOMATileDBContext],
    *,
    use_relative_uri: Optional[bool],
    ingestion_params: IngestionParams,
) -> None:
    """
    Ingest an uns string array. In the SOMA data model, we have NDArrays _of number only_ ...
    so we need to make this a SOMADataFrame.

    Ideally we don't want to an index column "soma_joinid" -- "index", maybe.
    However, ``SOMADataFrame`` _requires_ that soma_joinid be present, either
    as an index column, or as a data column. The former is less confusing.
    """
    if len(value.shape) != 1:
        msg = (
            f"Skipped {coll.uri}[{key!r}]"
            f" (uns object): string-array is not one-dimensional"
        )
        logging.log_io(msg, msg)
        return

    n = len(value)
    df_uri = _util.uri_joinpath(coll.uri, key)
    df = pd.DataFrame(
        data={
            "soma_joinid": np.arange(n, dtype=np.int64),
            "values": [str(e) if e else "" for e in value],
        }
    )
    df.set_index("soma_joinid", inplace=True)

    with _write_dataframe_impl(
        df,
        df_uri,
        ingestion_params=ingestion_params,
        platform_config=platform_config,
        context=context,
    ) as soma_df:
        _maybe_set(coll, key, soma_df, use_relative_uri=use_relative_uri)


def _ingest_uns_ndarray(
    coll: AnyTileDBCollection,
    key: str,
    value: NPNDArray,
    platform_config: Optional[PlatformConfig],
    context: Optional[SOMATileDBContext],
    *,
    use_relative_uri: Optional[bool],
    ingestion_params: IngestionParams,
) -> None:
    arr_uri = _util.uri_joinpath(coll.uri, key)

    if any(e <= 0 for e in value.shape):
        msg = f"Skipped {arr_uri} (uns ndarray): zero in shape {value.shape}"
        logging.log_io(msg, msg)
        return

    try:
        pa_dtype = pa.from_numpy_dtype(value.dtype)
    except pa.ArrowNotImplementedError:
        msg = (
            f"Skipped {arr_uri} (uns ndarray):"
            f" unsupported dtype {value.dtype!r} ({value.dtype})"
        )
        logging.log_io(msg, msg)
        return
    try:
        soma_arr = _factory.open(arr_uri, "w", soma_type=DenseNDArray, context=context)
    except DoesNotExistError:
        soma_arr = DenseNDArray.create(
            arr_uri,
            type=pa_dtype,
            shape=value.shape,
            platform_config=platform_config,
            context=context,
        )

    # If resume mode: don't re-write existing data. This is the user's explicit request
    # that we not re-write things that have already been written.
    if ingestion_params.skip_existing_nonempty_domain:
        storage_ned = _read_nonempty_domain(soma_arr)
        dim_range = ((0, value.shape[0] - 1),)
        if _chunk_is_contained_in(dim_range, storage_ned):
            logging.log_io(
                f"Skipped {soma_arr.uri}",
                f"Skipped {soma_arr.uri}",
            )
            return

    with soma_arr:
        _maybe_set(coll, key, soma_arr, use_relative_uri=use_relative_uri)
        soma_arr.write(
            (),
            pa.Tensor.from_numpy(value),
            platform_config=platform_config,
        )
    msg = f"Wrote   {soma_arr.uri} (uns ndarray)"
    logging.log_io(msg, msg)


# ----------------------------------------------------------------
def to_h5ad(
    experiment: Experiment,
    h5ad_path: Path,
    measurement_name: str,
    *,
    X_layer_name: str = "data",
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
) -> None:
    """Converts the experiment group to `AnnData <https://anndata.readthedocs.io/>`_
    format and writes it to the specified ``.h5ad`` file.

    Lifecycle:
        Experimental.
    """
    s = _util.get_start_stamp()
    logging.log_io(None, f"START  Experiment.to_h5ad -> {h5ad_path}")

    anndata = to_anndata(
        experiment,
        measurement_name=measurement_name,
        obs_id_name=obs_id_name,
        var_id_name=var_id_name,
        X_layer_name=X_layer_name,
    )

    s2 = _util.get_start_stamp()
    logging.log_io(None, f"START  write {h5ad_path}")

    anndata.write_h5ad(h5ad_path)

    logging.log_io(None, _util.format_elapsed(s2, f"FINISH write {h5ad_path}"))

    logging.log_io(
        None, _util.format_elapsed(s, f"FINISH Experiment.to_h5ad -> {h5ad_path}")
    )


# ----------------------------------------------------------------
def to_anndata(
    experiment: Experiment,
    measurement_name: str,
    *,
    X_layer_name: str = "data",
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
) -> ad.AnnData:
    """Converts the experiment group to `AnnData <https://anndata.readthedocs.io/>`_
    format. Choice of matrix formats is following what we often see in input
    ``.h5ad`` files:

    * ``X`` as ``scipy.sparse.csr_matrix``
    * ``obs``,``var`` as ``pandas.dataframe``
    * ``obsm``,``varm`` arrays as ``numpy.ndarray``
    * ``obsp``,``varp`` arrays as ``scipy.sparse.csr_matrix``

    Lifecycle:
        Experimental.
    """

    s = _util.get_start_stamp()
    logging.log_io(None, "START  Experiment.to_anndata")

    if measurement_name not in experiment.ms.keys():
        raise ValueError(
            f"requested measurement name {measurement_name} not found in input: {experiment.ms.keys()}"
        )
    measurement = experiment.ms[measurement_name]

    obs_df = experiment.obs.read().concat().to_pandas()
    obs_df.drop([SOMA_JOINID], axis=1, inplace=True)
    if obs_id_name not in obs_df.keys():
        raise ValueError(
            f"requested obs IDs column name {obs_id_name} not found in input: {obs_df.keys()}"
        )
    obs_df.set_index(obs_id_name, inplace=True)

    var_df = measurement.var.read().concat().to_pandas()
    var_df.drop([SOMA_JOINID], axis=1, inplace=True)
    if var_id_name not in var_df.keys():
        raise ValueError(
            f"requested var IDs column name {var_id_name} not found in input: {var_df.keys()}"
        )
    var_df.set_index(var_id_name, inplace=True)

    nobs = len(obs_df.index)
    nvar = len(var_df.index)

    if X_layer_name not in measurement.X:
        raise ValueError(
            f"X_layer_name {X_layer_name} not found in data: {measurement.X.keys()}"
        )
    X_data = measurement.X[X_layer_name]
    X_csr = None
    X_dtype = None  # some datasets have no X
    if isinstance(X_data, DenseNDArray):
        X_ndarray = X_data.read((slice(None), slice(None))).to_numpy()
        X_dtype = X_ndarray.dtype
    elif isinstance(X_data, SparseNDArray):
        X_mat = X_data.read().tables().concat().to_pandas()  # TODO: CSR/CSC options ...
        X_csr = conversions.csr_from_tiledb_df(X_mat, nobs, nvar)
        X_dtype = X_csr.dtype
    else:
        raise TypeError(f"Unexpected NDArray type {type(X_data)}")

    obsm = {}
    if "obsm" in measurement:
        for key in measurement.obsm.keys():
            shape = measurement.obsm[key].shape
            if len(shape) != 2:
                raise ValueError(f"expected shape == 2; got {shape}")
            if isinstance(measurement.obsm[key], DenseNDArray):
                obj = cast(DenseNDArray, measurement.obsm[key])
                matrix = obj.read().to_numpy()
                # The spelling ``sp.csr_array`` is more idiomatic but doesn't exist until Python 3.8
                obsm[key] = matrix
            else:
                # obsp is nobs x nobs.
                # obsm is nobs x some number -- number of PCA components, etc.
                matrix = measurement.obsm[key].read().tables().concat().to_pandas()
                nobs_times_width, coo_column_count = matrix.shape
                if coo_column_count != 3:
                    raise SOMAError(
                        f"internal error: expect COO width of 3; got {coo_column_count}"
                    )
                if nobs_times_width % nobs != 0:
                    raise SOMAError(
                        f"internal error: encountered non-rectangular obsm[{key}]: {nobs} does not divide {nobs_times_width}"
                    )
                obsm[key] = conversions.csr_from_tiledb_df(
                    matrix, nobs, nobs_times_width // nobs
                ).toarray()

    varm = {}
    if "varm" in measurement:
        for key in measurement.varm.keys():
            shape = measurement.varm[key].shape
            if len(shape) != 2:
                raise ValueError(f"expected shape == 2; got {shape}")
            if isinstance(measurement.varm[key], DenseNDArray):
                obj = cast(DenseNDArray, measurement.varm[key])
                matrix = obj.read().to_numpy()
                # The spelling ``sp.csr_array`` is more idiomatic but doesn't exist until Python 3.8
                varm[key] = matrix
            else:
                # varp is nvar x nvar.
                # varm is nvar x some number -- number of PCs, etc.
                matrix = measurement.varm[key].read().tables().concat().to_pandas()
                nvar_times_width, coo_column_count = matrix.shape
                if coo_column_count != 3:
                    raise SOMAError(
                        f"internal error: expect COO width of 3; got {coo_column_count}"
                    )
                if nvar_times_width % nvar != 0:
                    raise SOMAError(
                        f"internal error: encountered non-rectangular varm[{key}]: {nvar} does not divide {nvar_times_width}"
                    )
                varm[key] = conversions.csr_from_tiledb_df(
                    matrix, nvar, nvar_times_width // nvar
                ).toarray()

    obsp = {}
    if "obsp" in measurement:
        for key in measurement.obsp.keys():
            matrix = measurement.obsp[key].read().tables().concat().to_pandas()
            obsp[key] = conversions.csr_from_tiledb_df(matrix, nobs, nobs)

    varp = {}
    if "varp" in measurement:
        for key in measurement.varp.keys():
            matrix = measurement.varp[key].read().tables().concat().to_pandas()
            varp[key] = conversions.csr_from_tiledb_df(matrix, nvar, nvar)

    anndata = ad.AnnData(
        X=X_csr if X_csr is not None else X_ndarray,
        obs=obs_df,
        var=var_df,
        obsm=obsm,
        varm=varm,
        obsp=obsp,
        varp=varp,
        dtype=X_dtype,
    )

    logging.log_io(None, _util.format_elapsed(s, "FINISH Experiment.to_anndata"))

    return anndata
