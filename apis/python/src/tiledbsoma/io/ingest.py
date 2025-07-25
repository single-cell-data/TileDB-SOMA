# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Ingestion methods.

This module contains methods to generate SOMA artifacts starting from
other formats. Currently only ``.h5ad`` (`AnnData <https://anndata.readthedocs.io/>`_) is supported.
"""

from __future__ import annotations

import contextlib
import json
import math
import multiprocessing
import os
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from functools import partial
from itertools import repeat
from typing import (
    Any,
    Iterable,
    Literal,
    Mapping,
    Sequence,
    TypedDict,
    TypeVar,
    cast,
    no_type_check,
    overload,
)

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
from more_itertools import batched
from typing_extensions import deprecated

# As of anndata 0.11 we get a warning importing anndata.experimental.
# But anndata.abc doesn't exist in anndata 0.10. And anndata 0.11 doesn't
# exist for Python 3.9. And we have not yet dropped support for Python 3.9.
try:
    from anndata.abc import CSCDataset
except (AttributeError, ModuleNotFoundError):
    from anndata.experimental import CSCDataset

from somacore.options import PlatformConfig
from typing_extensions import get_args

from tiledbsoma import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    Measurement,
    PointCloudDataFrame,
    SparseNDArray,
    _factory,
    _util,
    eta,
    logging,
)
from tiledbsoma._collection import AnyTileDBCollection, CollectionBase
from tiledbsoma._common_nd_array import NDArray
from tiledbsoma._constants import SOMA_DATAFRAME_ORIGINAL_INDEX_NAME_JSON, SOMA_JOINID
from tiledbsoma._exception import AlreadyExistsError, DoesNotExistError, NotCreateableError, SOMAError
from tiledbsoma._soma_array import SOMAArray
from tiledbsoma._soma_object import AnySOMAObject, SOMAObject
from tiledbsoma._tdb_handles import RawHandle
from tiledbsoma._types import _INGEST_MODES, INGEST_MODES, IngestMode, NPNDArray, Path, _IngestMode
from tiledbsoma.options import SOMATileDBContext
from tiledbsoma.options._soma_tiledb_context import _validate_soma_tiledb_context
from tiledbsoma.options._tiledb_create_write_options import TileDBCreateOptions, TileDBWriteOptions

from . import conversions
from ._common import (
    _TILEDBSOMA_TYPE,
    _UNS_OUTGEST_COLUMN_NAME_1D,
    _UNS_OUTGEST_HINT_1D,
    _UNS_OUTGEST_HINT_2D,
    _UNS_OUTGEST_HINT_KEY,
    AdditionalMetadata,
    Matrix,
    SparseMatrix,
    UnsMapping,
    UnsNode,
)
from ._registration import (
    AxisIDMapping,
    ExperimentAmbientLabelMapping,
    ExperimentIDMapping,
)
from ._util import get_arrow_str_format, read_h5ad

_NDArr = TypeVar("_NDArr", bound=NDArray)
_TDBO = TypeVar("_TDBO", bound=SOMAObject[RawHandle])


def add_metadata(obj: SOMAObject[Any], additional_metadata: AdditionalMetadata) -> None:
    if additional_metadata:
        obj.verify_open_for_writing()
        obj.metadata.update(additional_metadata)


# ----------------------------------------------------------------
class IngestionParams:
    """Maps from user-level ingest modes to a set of implementation-level boolean flags."""

    write_schema_no_data: bool
    error_if_already_exists: bool
    skip_existing_nonempty_domain: bool
    appending: bool

    def __init__(
        self,
        ingest_mode: _IngestMode,
        label_mapping: ExperimentAmbientLabelMapping | None,
    ) -> None:
        if ingest_mode == "schema_only":
            self.write_schema_no_data = True
            self.error_if_already_exists = False
            self.skip_existing_nonempty_domain = False
            self.appending = False

        elif ingest_mode == "write":
            if label_mapping is None:
                self.write_schema_no_data = False
                self.error_if_already_exists = True
                self.skip_existing_nonempty_domain = False
                self.appending = False
            else:
                # append mode, but, the user supplying non-null registration information suffices
                # for us to understand "append"
                self.write_schema_no_data = False
                self.error_if_already_exists = False
                self.skip_existing_nonempty_domain = False
                self.appending = True

        elif ingest_mode == "resume":
            if label_mapping is None:
                self.write_schema_no_data = False
                self.error_if_already_exists = False
                self.skip_existing_nonempty_domain = True
                self.appending = False
            else:
                # resume-append mode, but, the user supplying non-null registration information
                # suffices for us to understand "resume-append"
                self.write_schema_no_data = False
                self.error_if_already_exists = False
                self.skip_existing_nonempty_domain = True
                self.appending = True

        elif ingest_mode == "update":
            self.write_schema_no_data = False
            self.error_if_already_exists = False
            self.skip_existing_nonempty_domain = False
            self.appending = False

        else:
            raise SOMAError(f'expected ingest_mode to be one of {_INGEST_MODES}; got "{ingest_mode}"')


# The tiledbsoma.io._registration package is private. These are the two sole user-facing API
# entrypoints for append-mode soma_joinid registration.
def register_h5ads(
    experiment_uri: str | None,
    h5ad_file_names: Sequence[str] | str,
    *,
    measurement_name: str,
    obs_field_name: str,
    var_field_name: str,
    append_obsm_varm: bool = False,
    context: SOMATileDBContext | None = None,
    use_multiprocessing: bool = False,
    allow_duplicate_obs_ids: bool = False,
) -> ExperimentAmbientLabelMapping:
    """Extends registration data from the baseline, already-written SOMA
    experiment to include multiple H5AD input files. See ``from_h5ad`` and
    ``from_anndata`` on-line help.

    The registration process will raise an error if any `obs` IDs (from `obs_field_name`)
    are duplicated across the combination of all inputs and the target SOMA Experiment.
    You can set `allow_duplicate_obs_ids=True` to bypass this check if you are adding a
    new Measurement to existing observations.

    If enabled via the ``use_multiprocessing`` parameter, this function will use multiprocessing
    to register each H5AD in parallel. In cases with many files, this can produce a performance
    benefit. Regardless of ``use_multiprocessing``, H5ADs will be registered concurrently -- you
    can control the concurrency using the ``soma.compute_concurrency_level`` configuration
    parameter in the ``context`` argument.
    """
    if isinstance(h5ad_file_names, str):
        h5ad_file_names = [h5ad_file_names]

    context = _validate_soma_tiledb_context(context)
    concurrency_level = _concurrency_level(context)

    logging.log_io(None, f"Loading per-axis metadata for {len(h5ad_file_names)} files.")
    executor_context: contextlib.AbstractContextManager[ProcessPoolExecutor | ThreadPoolExecutor]
    if use_multiprocessing:
        if multiprocessing.get_start_method() == "fork":
            warnings.warn(
                "Multiprocessing `fork` start method is inherently unsafe -- use `spawn`. See `multiprocessing.set_start_method()`",
                stacklevel=2,
            )
        executor_context = ProcessPoolExecutor(max_workers=concurrency_level)
    else:
        executor_context = contextlib.nullcontext(enter_result=context.threadpool)

    with executor_context as executor:
        axes_metadata = list(
            executor.map(
                ExperimentAmbientLabelMapping._load_axes_metadata_from_h5ads,
                batched(
                    h5ad_file_names,
                    math.ceil(len(h5ad_file_names) / concurrency_level),
                ),
                repeat(obs_field_name),
                repeat(var_field_name),
                repeat(
                    partial(
                        ExperimentAmbientLabelMapping._validate_anndata,
                        append_obsm_varm,
                    ),
                ),
            ),
        )
    logging.log_io(None, "Loaded per-axis metadata")

    return ExperimentAmbientLabelMapping._register_common(
        experiment_uri,
        axes_metadata,
        measurement_name=measurement_name,
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
        context=context,
        allow_duplicate_obs_ids=allow_duplicate_obs_ids,
    )


def register_anndatas(
    experiment_uri: str | None,
    adatas: Iterable[ad.AnnData] | ad.AnnData,
    *,
    measurement_name: str,
    obs_field_name: str,
    var_field_name: str,
    append_obsm_varm: bool = False,
    context: SOMATileDBContext | None = None,
    allow_duplicate_obs_ids: bool = False,
) -> ExperimentAmbientLabelMapping:
    """Extends registration data from the baseline, already-written SOMA
    experiment to include multiple H5AD input files. See ``from_h5ad`` and
    ``from_anndata`` on-line help.
    """
    if isinstance(adatas, ad.AnnData):
        adatas = [adatas]

    context = _validate_soma_tiledb_context(context)

    axes_metadata = [
        ExperimentAmbientLabelMapping._load_axes_metadata_from_anndatas(
            adatas,
            obs_field_name,
            var_field_name,
            partial(ExperimentAmbientLabelMapping._validate_anndata, append_obsm_varm),
        ),
    ]

    return ExperimentAmbientLabelMapping._register_common(
        experiment_uri,
        axes_metadata,
        measurement_name=measurement_name,
        obs_field_name=obs_field_name,
        var_field_name=var_field_name,
        context=context,
        allow_duplicate_obs_ids=allow_duplicate_obs_ids,
    )


def from_h5ad(
    experiment_uri: str,
    input_path: Path,
    measurement_name: str,
    *,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None = None,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    X_layer_name: str = "data",
    raw_X_layer_name: str = "data",
    ingest_mode: IngestMode = "write",
    use_relative_uri: bool | None = None,
    X_kind: type[SparseNDArray] | type[DenseNDArray] = SparseNDArray,
    registration_mapping: ExperimentAmbientLabelMapping | None = None,
    uns_keys: Sequence[str] | None = None,
    additional_metadata: AdditionalMetadata = None,
) -> str:
    r"""Reads an ``.h5ad`` file and writes it to an :class:`Experiment`.

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
          name in the input AnnData object's ``obs``/``var``.

        X_layer_name: SOMA array name for the AnnData's ``X`` matrix.

        raw_X_layer_name: SOMA array name for the AnnData's ``raw/X`` matrix.

        ingest_mode: The ingestion type to perform:

            - ``write``: Writes all data, creating new layers if the SOMA already exists.
            - ``resume``: (deprecated) Adds data to an existing SOMA, skipping writing data
              that was previously written. Useful for continuing after a partial
              or interrupted ingestion operation.
            - ``schema_only``: Creates groups and the array schema, without
              writing any data to the array. Useful to prepare for appending
              multiple H5AD files to a single SOMA.

          The 'resume' ingest_mode is deprecated and will be removed in a future version. The
          current implementation has a known issue that can can cause multi-dataset appends to
          not resume correctly.

          The recommended and safest approach for recovering from a failed ingestion is to delete
          the partially written SOMA Experiment and restart the ingestion process from the original
          input files or a known-good backup.

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

    Returns:
        The URI of the newly created experiment.

    Lifecycle:
        Maturing.
    """
    if ingest_mode not in INGEST_MODES:
        raise SOMAError(f'expected ingest_mode to be one of {INGEST_MODES}; got "{ingest_mode}"')
    _check_for_deprecated_modes(ingest_mode)

    if isinstance(input_path, ad.AnnData):
        raise TypeError("input path is an AnnData object -- did you want from_anndata?")

    context = _validate_soma_tiledb_context(context)

    s = _util.get_start_stamp()
    logging.log_io(None, f"START  Experiment.from_h5ad {input_path}")

    logging.log_io(None, f"START  READING {input_path}")

    with read_h5ad(input_path, mode="r", ctx=context) as anndata:
        logging.log_io(None, _util.format_elapsed(s, f"FINISH READING {input_path}"))

        uri = _from_anndata(
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

    logging.log_io(None, _util.format_elapsed(s, f"FINISH Experiment.from_h5ad {input_path} {uri}"))
    return uri


class IngestCtx(TypedDict):
    """Convenience type-alias for kwargs passed to ingest functions."""

    context: SOMATileDBContext | None
    ingestion_params: IngestionParams
    additional_metadata: AdditionalMetadata


class IngestPlatformCtx(IngestCtx):
    """Convenience type-alias for kwargs passed to ingest functions.

    Extends :class:`IngestCtx`, adds ``platform_config``.
    """

    platform_config: PlatformConfig | None


# ----------------------------------------------------------------
def from_anndata(
    experiment_uri: str,
    anndata: ad.AnnData,
    measurement_name: str,
    *,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None = None,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    X_layer_name: str = "data",
    raw_X_layer_name: str = "data",
    ingest_mode: IngestMode = "write",
    use_relative_uri: bool | None = None,
    X_kind: type[SparseNDArray] | type[DenseNDArray] = SparseNDArray,
    registration_mapping: ExperimentAmbientLabelMapping | None = None,
    uns_keys: Sequence[str] | None = None,
    additional_metadata: AdditionalMetadata = None,
) -> str:
    """Writes an `AnnData <https://anndata.readthedocs.io/>`_ object to an :class:`Experiment`.

    Usage is the same as ``from_h5ad`` except that you can use this function when the AnnData object
    is already loaded into memory.

    Lifecycle:
        Maturing.
    """
    if ingest_mode not in INGEST_MODES:
        raise SOMAError(f'expected ingest_mode to be one of {INGEST_MODES}; got "{ingest_mode}"')
    _check_for_deprecated_modes(ingest_mode)

    return _from_anndata(
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


def _from_anndata(
    experiment_uri: str,
    anndata: ad.AnnData,
    measurement_name: str,
    *,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None = None,
    obs_id_name: str = "obs_id",
    var_id_name: str = "var_id",
    X_layer_name: str = "data",
    raw_X_layer_name: str = "data",
    ingest_mode: IngestMode = "write",
    use_relative_uri: bool | None = None,
    X_kind: type[SparseNDArray] | type[DenseNDArray] = SparseNDArray,
    registration_mapping: ExperimentAmbientLabelMapping | None = None,
    uns_keys: Sequence[str] | None = None,
    additional_metadata: AdditionalMetadata = None,
) -> str:
    """Private helper function."""
    # Map the user-level ingest mode to a set of implementation-level boolean flags
    ingestion_params = IngestionParams(ingest_mode, registration_mapping)

    if ingestion_params.appending and X_kind == DenseNDArray:
        raise ValueError("dense X is not supported for append mode")

    if not isinstance(anndata, ad.AnnData):
        raise TypeError("Second argument is not an AnnData object -- did you want from_h5ad?")

    for ad_key in ["obsm", "obsp", "varm", "varp"]:
        for key, val in getattr(anndata, ad_key).items():
            if not isinstance(val, get_args(Matrix)):
                raise TypeError(
                    f"{ad_key} value at {key} is not of type {[cl.__name__ for cl in get_args(Matrix)]}: {type(val)}",
                )

    # For single ingest (no append):
    #
    # * obs, var, X, etc. array indices map 1-1 to soma_joinid, soma_dim_0, soma_dim_1, etc.
    #
    # * There is no need to look at a particular barcode column etc.
    #
    # For append mode:
    #
    # * We require mappings to be pre-computed and passed in
    #
    # * The mappings are from obs-label/var-label to soma_joinid, soma_dim_0, soma_dim_1
    #   for all the H5AD/AnnData objects being ingested
    #
    # * Here we select out the renumberings for the obs, var, X, etc. array indices
    if registration_mapping is None:
        joinid_maps = ExperimentIDMapping.from_anndata(anndata, measurement_name=measurement_name)
        filter_existing_obs_joinids: bool = True
    else:
        if not registration_mapping.prepared and Experiment.exists(experiment_uri):
            raise SOMAError(
                "Experiment must be prepared prior to ingestion. Please call ``registration_map.prepare_experiment`` method.",
            )
        joinid_maps = registration_mapping.id_mappings_for_anndata(anndata, measurement_name=measurement_name)
        filter_existing_obs_joinids = registration_mapping.obs_axis.allow_duplicate_ids

    context = _validate_soma_tiledb_context(context)

    # Without _at least_ one index, there is nothing to indicate the dimension indices.
    if anndata.obs.index.empty or anndata.var.index.empty:
        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

    s = _util.get_start_stamp()
    logging.log_io(None, "START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()

    logging.log_io(None, _util.format_elapsed(s, "FINISH DECATEGORICALIZING"))

    s = _util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {experiment_uri}")

    ingest_ctx: IngestCtx = {
        "context": context,
        "ingestion_params": ingestion_params,
        "additional_metadata": additional_metadata,
    }
    ingest_platform_ctx: IngestPlatformCtx = dict(**ingest_ctx, platform_config=platform_config)

    # Must be done first, to create the parent directory.
    experiment = _create_or_open_collection(Experiment, experiment_uri, **ingest_ctx)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # OBS
    df_uri = _util.uri_joinpath(experiment_uri, "obs")
    with _write_dataframe(
        df_uri,
        conversions.obs_or_var_to_tiledb_supported_array_type(anndata.obs),
        id_column_name=obs_id_name,
        axis_mapping=joinid_maps.obs_axis,
        **ingest_platform_ctx,
        filter_existing_joinids=filter_existing_obs_joinids,
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

    with _create_or_open_collection(
        Collection[Measurement],
        experiment_ms_uri,
        **ingest_ctx,
    ) as ms:
        _maybe_set(experiment, "ms", ms, use_relative_uri=use_relative_uri)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # MS/meas
        measurement_uri = _util.uri_joinpath(experiment_ms_uri, _util.sanitize_key(measurement_name))
        with _create_or_open_collection(
            Measurement,
            measurement_uri,
            **ingest_ctx,
        ) as measurement:
            _maybe_set(ms, measurement_name, measurement, use_relative_uri=use_relative_uri)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # ms/meas/uns
            _maybe_ingest_uns(
                measurement,
                anndata.uns,
                use_relative_uri=use_relative_uri,
                uns_keys=uns_keys,
                **ingest_platform_ctx,
            )

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # MS/meas/VAR
            with _write_dataframe(
                _util.uri_joinpath(measurement_uri, "var"),
                conversions.obs_or_var_to_tiledb_supported_array_type(anndata.var),
                id_column_name=var_id_name,
                # Layer existence is pre-checked in the registration phase
                axis_mapping=joinid_maps.var_axes[measurement_name],
                **ingest_platform_ctx,
            ) as var:
                _maybe_set(measurement, "var", var, use_relative_uri=use_relative_uri)

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            # MS/meas/X/DATA
            measurement_X_uri = _util.uri_joinpath(measurement_uri, "X")
            with _create_or_open_collection(
                Collection,
                measurement_X_uri,
                **ingest_ctx,
            ) as x:
                _maybe_set(measurement, "X", x, use_relative_uri=use_relative_uri)

                # Since we did ``anndata = ad.read_h5ad(path_to_h5ad, "r")`` with the "r":
                # * If we do ``anndata.X[:]`` we're loading all of a CSR/CSC/etc into memory.
                # * If we do ``anndata.X`` we're getting a pageable object which can be loaded
                #   chunkwise into memory.
                # Using the latter allows us to ingest larger .h5ad files without OOMing.

                # Some AnnData objects have no X at all. This might be a missing attribute
                # anndata.X; it might be anndata.X present, but with value None.
                try:
                    has_X = anndata.X is not None
                except (NameError, KeyError):
                    # We need to check both -- different exception types occur depending
                    # on whether the anndata object is read in backing mode or not.
                    has_X = False

                if has_X:
                    with _create_from_matrix(
                        X_kind,
                        _util.uri_joinpath(measurement_X_uri, _util.sanitize_key(X_layer_name)),
                        anndata.X,
                        axis_0_mapping=joinid_maps.obs_axis,
                        axis_1_mapping=joinid_maps.var_axes[measurement_name],
                        **ingest_platform_ctx,
                    ) as data:
                        _maybe_set(x, X_layer_name, data, use_relative_uri=use_relative_uri)

                for layer_name, layer in anndata.layers.items():
                    with _create_from_matrix(
                        X_kind,
                        _util.uri_joinpath(measurement_X_uri, _util.sanitize_key(layer_name)),
                        layer,
                        axis_0_mapping=joinid_maps.obs_axis,
                        axis_1_mapping=joinid_maps.var_axes[measurement_name],
                        **ingest_platform_ctx,
                    ) as layer_data:
                        _maybe_set(x, layer_name, layer_data, use_relative_uri=use_relative_uri)

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # MS/meas/OBSM,VARM,OBSP,VARP
                def _ingest_obs_var_m_p(
                    ad_key: Literal["obsm", "varm", "obsp", "varp"],
                    axis_0_mapping: AxisIDMapping,
                    axis_1_mapping: AxisIDMapping | None = None,
                ) -> None:
                    ad_val = getattr(anndata, ad_key)
                    if len(ad_val.keys()) > 0:  # do not create an empty collection
                        ad_val_uri = _util.uri_joinpath(measurement_uri, _util.sanitize_key(ad_key))
                        with _create_or_open_collection(
                            Collection,
                            ad_val_uri,
                            **ingest_ctx,
                        ) as coll:
                            _maybe_set(
                                measurement,
                                ad_key,
                                coll,
                                use_relative_uri=use_relative_uri,
                            )
                            for key in ad_val:
                                val = ad_val[key]
                                num_cols = val.shape[1]
                                axis_1_mapping_ = axis_1_mapping if axis_1_mapping else AxisIDMapping.identity(num_cols)
                                with _create_from_matrix(
                                    # TODO (https://github.com/single-cell-data/TileDB-SOMA/issues/1245):
                                    # consider a use-dense flag at the tiledbsoma.io API
                                    # DenseNDArray,
                                    SparseNDArray,
                                    _util.uri_joinpath(ad_val_uri, _util.sanitize_key(key)),
                                    conversions.to_tiledb_supported_array_type(key, val),
                                    axis_0_mapping=axis_0_mapping,
                                    axis_1_mapping=axis_1_mapping_,
                                    **ingest_platform_ctx,
                                ) as arr:
                                    _maybe_set(
                                        coll,
                                        key,
                                        arr,
                                        use_relative_uri=use_relative_uri,
                                    )

                _ingest_obs_var_m_p("obsm", joinid_maps.obs_axis)
                _ingest_obs_var_m_p("varm", joinid_maps.var_axes[measurement_name])
                _ingest_obs_var_m_p("obsp", joinid_maps.obs_axis, joinid_maps.obs_axis)
                _ingest_obs_var_m_p(
                    "varp",
                    joinid_maps.var_axes[measurement_name],
                    joinid_maps.var_axes[measurement_name],
                )

                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # MS/RAW
                if anndata.raw is not None:
                    raw_uri = _util.uri_joinpath(experiment_ms_uri, "raw")
                    with _create_or_open_collection(
                        Measurement,
                        raw_uri,
                        **ingest_ctx,
                    ) as raw_measurement:
                        _maybe_set(
                            ms,
                            "raw",
                            raw_measurement,
                            use_relative_uri=use_relative_uri,
                        )

                        with _write_dataframe(
                            _util.uri_joinpath(raw_uri, "var"),
                            conversions.obs_or_var_to_tiledb_supported_array_type(anndata.raw.var),
                            id_column_name=var_id_name,
                            axis_mapping=joinid_maps.var_axes["raw"],
                            **ingest_platform_ctx,
                        ) as var:
                            _maybe_set(
                                raw_measurement,
                                "var",
                                var,
                                use_relative_uri=use_relative_uri,
                            )

                        raw_X_uri = _util.uri_joinpath(raw_uri, "X")
                        with _create_or_open_collection(
                            Collection,
                            raw_X_uri,
                            **ingest_ctx,
                        ) as rm_x:
                            _maybe_set(
                                raw_measurement,
                                "X",
                                rm_x,
                                use_relative_uri=use_relative_uri,
                            )

                            with _create_from_matrix(
                                SparseNDArray,
                                _util.uri_joinpath(raw_X_uri, _util.sanitize_key(raw_X_layer_name)),
                                anndata.raw.X,
                                axis_0_mapping=joinid_maps.obs_axis,
                                axis_1_mapping=joinid_maps.var_axes["raw"],
                                **ingest_platform_ctx,
                            ) as rm_x_data:
                                _maybe_set(
                                    rm_x,
                                    raw_X_layer_name,
                                    rm_x_data,
                                    use_relative_uri=use_relative_uri,
                                )

    experiment.close()

    logging.log_io(
        f"Wrote   {experiment.uri}",
        _util.format_elapsed(s, f"FINISH WRITING {experiment.uri}"),
    )
    return experiment.uri


@deprecated(
    """This function is deprecated and will be removed in a future version of this package.

It is recommended to use tiledbsoma.io.from_anndata (with a registration map from tiledbsoma.io.register_anndatas or tiledbsoma.io.register_h5ads) for appending new, complete AnnData objects to an Experiment.""",
)
def append_obs(
    exp: Experiment,
    new_obs: pd.DataFrame,
    *,
    obs_id_name: str = "obs_id",
    registration_mapping: ExperimentAmbientLabelMapping,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None = None,
) -> str:
    """Writes new rows to an existing ``obs`` dataframe (this is distinct from ``update_obs``
    which mutates the entirety of the ``obs`` dataframe, e.g. to add/remove columns).

    This function is deprecated and will be removed in a future version of this package.

    It is recommended to use ``tiledbsoma.io.from_anndata`` (with a registration map from
    ``tiledbsoma.io.register_anndatas`` or ``tiledbsoma.io.register_h5ads``) for appending new,
    complete AnnData objects to an :class:`Experiment`.

    Example::

        rd = tiledbsoma.io.register_anndatas(
            exp_uri,
            [new_anndata],
            measurement_name="RNA",
            obs_field_name="obs_id",
            var_field_name="var_id",
        )

        with tiledbsoma.Experiment.open(exp_uri, "w") as exp:
            tiledbsoma.io.append_obs(
                exp, new_anndata.obs, registration_mapping=rd,
            )

    Lifecycle:
        Deprecated.
    """
    exp.verify_open_for_writing()

    # Map the user-level ingest mode to a set of implementation-level boolean flags.
    # See comments in from_anndata.
    context = _validate_soma_tiledb_context(context)
    ingestion_params = IngestionParams("write", registration_mapping)
    joinid_map = registration_mapping.obs_axis.id_mapping_from_dataframe(new_obs)

    s = _util.get_start_stamp()
    logging.log_io_same(f"Start  writing obs for {exp.obs.uri}")

    with _write_dataframe(
        exp.obs.uri,
        conversions.obs_or_var_to_tiledb_supported_array_type(new_obs),
        id_column_name=obs_id_name,
        platform_config=platform_config,
        context=context,
        ingestion_params=ingestion_params,
        axis_mapping=joinid_map,
        must_exist=True,
    ):
        logging.log_io_same(_util.format_elapsed(s, f"Finish writing obs for {exp.obs.uri}"))
    return exp.obs.uri


@deprecated(
    """This function is deprecated and will be removed in a future version of this package.

It is recommended to use tiledbsoma.io.from_anndata (with a registration map from tiledbsoma.io.register_anndatas or tiledbsoma.io.register_h5ads) for appending new, complete AnnData objects to an Experiment.""",
)
def append_var(
    exp: Experiment,
    new_var: pd.DataFrame,
    measurement_name: str,
    *,
    var_id_name: str = "var_id",
    registration_mapping: ExperimentAmbientLabelMapping,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None = None,
) -> str:
    """Writes new rows to an existing ``var`` dataframe (this is distinct from ``update_var``
    which mutates the entirety of the ``var`` dataframe, e.g. to add/remove columns).

    This function is deprecated and will be removed in a future version of this package.

    It is recommended to use ``tiledbsoma.io.from_anndata`` (with a registration map from
    ``tiledbsoma.io.register_anndatas`` or ``tiledbsoma.io.register_h5ads``) for appending new,
    complete AnnData objects to an :class:`Experiment`.

    Example::

        rd = tiledbsoma.io.register_anndatas(
            exp_uri,
            [new_anndata],
            measurement_name="RNA",
            obs_field_name="obs_id",
            var_field_name="var_id",
        )

        with tiledbsoma.Experiment.open(exp_uri, "w") as exp:
            tiledbsoma.io.append_var(
                exp, a2.var, measurement_name="RNA", registration_mapping=rd,
            )

    Lifecycle:
        Deprecated.
    """
    exp.verify_open_for_writing()
    if measurement_name not in exp.ms:
        raise SOMAError(f"Experiment {exp.uri} has no measurement named {measurement_name}")
    sdf = exp.ms[measurement_name].var

    # Map the user-level ingest mode to a set of implementation-level boolean flags.
    # See comments in from_anndata.
    context = _validate_soma_tiledb_context(context)
    ingestion_params = IngestionParams("write", registration_mapping)
    joinid_map = registration_mapping.var_axes[measurement_name].id_mapping_from_dataframe(new_var)

    s = _util.get_start_stamp()
    logging.log_io_same(f"Start  writing var for {sdf.uri}")

    with _write_dataframe(
        sdf.uri,
        conversions.obs_or_var_to_tiledb_supported_array_type(new_var),
        id_column_name=var_id_name,
        platform_config=platform_config,
        context=context,
        ingestion_params=ingestion_params,
        axis_mapping=joinid_map,
        must_exist=True,
    ):
        logging.log_io_same(_util.format_elapsed(s, f"Finish writing var for {sdf.uri}"))
    return sdf.uri


@deprecated(
    """This function is deprecated and will be removed in a future version of this package.

It is recommended to use tiledbsoma.io.from_anndata (with a registration map from tiledbsoma.io.register_anndatas or tiledbsoma.io.register_h5ads) for appending new, complete AnnData objects to an Experiment.""",
)
def append_X(
    exp: Experiment,
    new_X: Matrix | h5py.Dataset,
    measurement_name: str,
    X_layer_name: str,
    obs_ids: Sequence[str],
    var_ids: Sequence[str],
    *,
    registration_mapping: ExperimentAmbientLabelMapping,
    X_kind: type[SparseNDArray] | type[DenseNDArray] = SparseNDArray,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None = None,
) -> str:
    """Appends new data to an existing ``X`` matrix. Nominally to be used in conjunction
    with ``update_obs`` and ``update_var``, as an itemized alternative to doing
    ``from_anndata`` with a registration mapping supplied.

    This function is deprecated and will be removed in a future version of this package.

    It is recommended to use ``tiledbsoma.io.from_anndata`` (with a registration map from
    ``tiledbsoma.io.register_anndatas`` or ``tiledbsoma.io.register_h5ads``) for appending new,
    complete AnnData objects to an :class:`Experiment`.

    Example::

        rd = tiledbsoma.io.register_anndatas(
            exp_uri,
            [new_anndata],
            measurement_name="RNA",
            obs_field_name="obs_id",
            var_field_name="var_id",
        )


        with tiledbsoma.Experiment.open(exp_uri) as exp:
            tiledbsoma.io.append_X(
                exp,
                new_X=adata.X,
                measurement_name=measurement_name,
                X_layer_name=X_layer_name,
                obs_ids=list(new_anndata.obs.index),
                var_ids=list(new_anndata.var.index),
                registration_mapping=rd,
            )

    Lifecycle:
        Deprecated.
    """
    exp.verify_open_for_writing()
    if measurement_name not in exp.ms:
        raise SOMAError(f"Experiment {exp.uri} has no measurement named {measurement_name}")
    if X_layer_name not in exp.ms[measurement_name].X:
        raise SOMAError(f"Experiment {exp.uri} has no X layer named {X_layer_name} in measurement {measurement_name}")
    X = exp.ms[measurement_name].X[X_layer_name]

    # Map the user-level ingest mode to a set of implementation-level boolean flags.
    # See comments in from_anndata.
    ingestion_params = IngestionParams("write", registration_mapping)
    context = _validate_soma_tiledb_context(context)

    s = _util.get_start_stamp()
    logging.log_io_same(f"Start  writing var for {X.uri}")

    axis_0_mapping = registration_mapping.obs_axis.id_mapping_from_values(obs_ids)
    axis_1_mapping = registration_mapping.var_axes[measurement_name].id_mapping_from_values(var_ids)

    with _create_from_matrix(
        X_kind,
        X.uri,
        new_X,
        ingestion_params=ingestion_params,
        platform_config=platform_config,
        context=context,
        axis_0_mapping=axis_0_mapping,
        axis_1_mapping=axis_1_mapping,
        must_exist=True,
    ):
        logging.log_io_same(_util.format_elapsed(s, f"Finish writing X for {X.uri}"))
    return X.uri


def _maybe_set(
    coll: AnyTileDBCollection,
    key: str,
    value: AnySOMAObject,
    *,
    use_relative_uri: bool | None,
) -> None:
    coll.verify_open_for_writing()
    try:
        coll.set(key, value, use_relative_uri=use_relative_uri)
    except SOMAError:
        # This is already a member of the collection.
        pass


@overload
def _create_or_open_collection(
    cls: type[Experiment],
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: SOMATileDBContext | None,
    additional_metadata: AdditionalMetadata = None,
) -> Experiment: ...


@overload
def _create_or_open_collection(
    cls: type[Measurement],
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: SOMATileDBContext | None,
    additional_metadata: AdditionalMetadata = None,
) -> Measurement: ...


@overload
def _create_or_open_collection(
    cls: type[Collection[_TDBO]],
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: SOMATileDBContext | None,
    additional_metadata: AdditionalMetadata = None,
) -> Collection[_TDBO]: ...


@no_type_check
def _create_or_open_collection(
    cls: type[CollectionBase[_TDBO]],
    uri: str,
    *,
    ingestion_params: IngestionParams,
    context: SOMATileDBContext | None,
    additional_metadata: AdditionalMetadata = None,
) -> CollectionBase[_TDBO]:
    try:
        coll = cls.create(uri, context=context)
    except (AlreadyExistsError, NotCreateableError) as e:
        # It already exists. Are we resuming?
        if ingestion_params.error_if_already_exists:
            raise SOMAError(f"{uri} already exists") from e
        coll = cls.open(uri, "w", context=context)

    add_metadata(coll, additional_metadata)
    return coll


# Spellings compatible with 1.2.7:
@overload
def _create_or_open_coll(
    cls: type[Experiment],
    uri: str,
    *,
    ingest_mode: IngestMode,
    context: SOMATileDBContext | None,
) -> Experiment: ...


@overload
def _create_or_open_coll(
    cls: type[Measurement],
    uri: str,
    *,
    ingest_mode: IngestMode,
    context: SOMATileDBContext | None,
) -> Measurement: ...


@overload
def _create_or_open_coll(
    cls: type[Collection[_TDBO]],
    uri: str,
    *,
    ingest_mode: IngestMode,
    context: SOMATileDBContext | None,
) -> Collection[_TDBO]: ...


def _create_or_open_coll(
    cls: type[Any],
    uri: str,
    *,
    ingest_mode: IngestMode,
    context: SOMATileDBContext | None,
) -> Any:
    return _create_or_open_collection(
        cls,
        uri,
        ingestion_params=IngestionParams(ingest_mode=ingest_mode, label_mapping=None),
        context=context,
    )


def _extract_new_values_for_append_aux(
    previous_soma_dataframe: DataFrame,
    arrow_table: pa.Table,
    filter_existing_joinids: bool,
) -> pa.Table:
    """Helper function for _extract_new_values_for_append.

    This does two things:

    * Retains only the 'new' rows compared to existing storage.  Example is
      append-mode updates to the var dataframe for which it's likely that
      most/all gene IDs have already been seen.

    * Categorical of type sometype vs plain sometype:

      o If we're appending a non-categorical column to existing categorical
        storage, convert the about-to-be-written data to categorical of that
        type, to match.

      o If we're appending a categorical column to existing non-categorical
        storage, convert the about-to-be-written data to non-categorical, to
        match.

    Context: https://github.com/single-cell-data/TileDB-SOMA/issues/3353.
    Namely, we find that AnnData's to_h5ad/from_h5ad can categoricalize (without
    the user's knowledge or intention) string columns.  For example, even
    cell_id/barcode, for which there may be millions of distinct values, with no
    gain to be had from dictionary encoding, will be converted to categorical.
    We find that converting these high-cardinality enums to non-enumerated is a
    significant performance win for subsequent accesses. When we do an initial
    ingest from AnnData to TileDB-SOMA, we decategoricalize if the cardinality
    exceeds some threshold.

    All well and good -- except for one more complication which is append mode.
    Namely, if the new column has high enough cardinality that we would
    downgrade to non-categorical, but the existing storage has categorical, we
    must write the new data as categorical.  Likewise, if the new column has low
    enough cardinality that we would keep it as categorical, but the existing
    storage has non-categorical, we must write the new data as non-categorical.
    """
    if filter_existing_joinids:
        # Retain only the new rows.
        previous_sjids_table = previous_soma_dataframe.read(column_names=["soma_joinid"]).concat()
        # use numpy.isin over pyarrow.compute.is_in, as it is MUCH faster
        mask = pa.array(
            np.isin(
                arrow_table[SOMA_JOINID].to_numpy(),
                previous_sjids_table[SOMA_JOINID].to_numpy(),
                invert=True,
            ),
        )

        arrow_table = arrow_table.filter(mask)

    # Check if any new data.
    if any(column.num_chunks == 0 for column in arrow_table.columns):
        return arrow_table

    # This is a redundant, failsafe check. The append-mode registrar already
    # ensure schema homogeneity before we get here.
    old_schema = previous_soma_dataframe.schema
    new_schema = arrow_table.schema

    # Note: we may be doing an add-column that doesn't exist in tiledbsoma
    # storage but is present in the new AnnData. We don't need to change
    # anything in that case. Regardless, we can't assume that the old
    # and new schema have the same column names.

    def is_cat(field: pa.Field) -> bool:
        return cast("bool", pa.types.is_dictionary(field.type))

    # Make a quick check of the old and new schemas to see if any columns need
    # changing between non-categorical and categorical.  We're about to
    # duplicate the new data -- and we must, since pyarrow.Table is immutable --
    # so let's only do that if we need to.
    any_to_change = False
    for name in new_schema.names:
        if name not in old_schema.names:
            continue
        old_field = old_schema.field(name)
        new_field = new_schema.field(name)
        if not is_cat(old_field) and is_cat(new_field):
            any_to_change = True
            break
        if is_cat(old_field) and not is_cat(new_field):
            any_to_change = True
            break

    if any_to_change:
        fields_dict = {}
        for name in new_schema.names:
            if name not in old_schema.names:
                continue
            column = arrow_table.column(name)
            old_field = old_schema.field(name)
            new_field = new_schema.field(name)

            # Note from https://enpiar.com/arrow-site/docs/python/data.html#tables
            # we have the assertion that pa.Table columns are necessarily of type
            # pa.ChunkedArray.
            assert isinstance(column, pa.ChunkedArray)

            if not is_cat(old_field) and is_cat(new_field):
                # Convert from categorical to non-categorical.  Note that if
                # this is a pa.ChunkedArray, there is no .dictionary_decode()
                # for it, but each chunk does have a .dictionary_decode().
                column = pa.chunked_array([chunk.dictionary_decode() for chunk in column.chunks])

            elif is_cat(old_field) and not is_cat(new_field):
                # Convert from non-categorical to categorical.  Note:
                # libtiledbsoma already merges the enum mappings, e.g if the
                # storage has red, yellow, & green, but our new data has some
                # yellow, green, and orange.
                column = pa.chunked_array([chunk.dictionary_encode() for chunk in column.chunks])

            fields_dict[name] = column
        arrow_table = pa.Table.from_pydict(fields_dict)

    return arrow_table


def _extract_new_values_for_append(
    df_uri: str,
    arrow_table: pa.Table,
    filter_existing_joinids: bool,
    context: SOMATileDBContext | None = None,
) -> pa.Table:
    """For append mode: mostly we just go ahead and write the data, except var.

    Nominally:

    * Cell IDs (obs IDs) are distinct per input
      o This means "obs grows down"
    * Gene IDs (var IDs) are almost all the same per input
      o This means "var gets stacked"
    * X is obs x var
      o This means "X grows down"
    * obsm is obs x M, and varm is var x N
      o This means means they "grow down"
    * obsp is obs x obs, and varp is var x var
      o This means they _would_ grow block-diagonally -- but
        (for biological reasons) we disallow append of obsp and varp

    Reviewing the above, we can see that "just write the data" is almost always the right way to
    append -- except for var. There, if there are 100 H5ADs appending into a single SOMA experiment,
    we don't want 100 TileDB fragment-layers of mostly the same data -- we just want to write the
    _new_ genes not yet seen before (if any).
    """
    try:
        with _factory.open(df_uri, "r", soma_type=DataFrame, context=context) as previous_soma_dataframe:
            return _extract_new_values_for_append_aux(previous_soma_dataframe, arrow_table, filter_existing_joinids)

    except DoesNotExistError:
        return arrow_table


def _write_arrow_table(
    arrow_table: pa.Table,
    handle: DataFrame | SparseNDArray | PointCloudDataFrame,
    tiledb_create_options: TileDBCreateOptions,
    tiledb_write_options: TileDBWriteOptions,
) -> None:
    """Handles num-bytes capacity for remote object stores."""
    cap = tiledb_create_options.remote_cap_nbytes
    if arrow_table.nbytes > cap:
        n = len(arrow_table)
        if n < 2:
            raise SOMAError(f"single table row nbytes {arrow_table.nbytes} exceeds cap nbytes {cap}")
        m = n // 2
        _write_arrow_table(arrow_table[:m], handle, tiledb_create_options, tiledb_write_options)
        _write_arrow_table(arrow_table[m:], handle, tiledb_create_options, tiledb_write_options)
    else:
        logging.log_io(
            None,
            f"Write Arrow table num_rows={len(arrow_table)} num_bytes={arrow_table.nbytes} cap={cap}",
        )
        handle.write(arrow_table, platform_config=tiledb_write_options)


def _write_dataframe(
    df_uri: str,
    df: pd.DataFrame,
    id_column_name: str | None,
    *,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
    platform_config: PlatformConfig | None = None,
    context: SOMATileDBContext | None = None,
    axis_mapping: AxisIDMapping,
    must_exist: bool = False,
    filter_existing_joinids: bool = True,
) -> DataFrame:
    """Convert and save a pd.DataFrame as a SOMA DataFrame.

    This function adds a ``soma_joinid`` index to the provided ``pd.DataFrame``, as that is
    required by libtiledbsoma; the original ``pd.DataFrame`` index is "reset" to a column, with its
    original name stored as metadata on the SOMA ``DataFrame``.

    There are several subtleties related to naming the reset index column, and the ``id_column_name``
     parameter; see ``_prepare_df_for_ingest`` / https://github.com/single-cell-data/TileDB-SOMA/issues/2829.

    NOTE: this function mutates the input dataframe, for parsimony of memory usage. Callers should
    copy any user-provided ``pd.DataFrame``s (adata obs, var, uns, etc.) before passing them here.
    """
    original_index_metadata = conversions._prepare_df_for_ingest(df, id_column_name)
    df[SOMA_JOINID] = np.asarray(axis_mapping.data, dtype=np.int64)
    df.set_index(SOMA_JOINID, inplace=True)

    return _write_dataframe_impl(
        df,
        df_uri,
        id_column_name,
        shape=axis_mapping.get_shape(),
        ingestion_params=ingestion_params,
        additional_metadata=additional_metadata,
        original_index_metadata=original_index_metadata,
        platform_config=platform_config,
        context=context,
        must_exist=must_exist,
        filter_existing_joinids=filter_existing_joinids,
    )


def _write_dataframe_impl(
    df: pd.DataFrame,
    df_uri: str,
    id_column_name: str | None,
    *,
    shape: int,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
    original_index_metadata: str | None = None,
    platform_config: PlatformConfig | None = None,
    context: SOMATileDBContext | None = None,
    must_exist: bool = False,
    filter_existing_joinids: bool = True,
) -> DataFrame:
    """Save a Pandas DataFrame as a SOMA DataFrame.

    Expects the required ``soma_joinid`` index to have already been added to the ``pd.DataFrame``.
    """
    s = _util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {df_uri}")

    arrow_table = conversions.df_to_arrow_table(df)

    # Don't many-layer the almost-always-repeated var dataframe.
    # And for obs, if there are duplicate values for whatever reason, don't write them
    # when appending.
    if ingestion_params.appending:
        if id_column_name is None:
            # Nominally, nil id_column_name only happens for uns append and we do not append uns,
            # which is a concern for our caller. This is a second-level check.
            raise ValueError("internal coding error: id_column_name unspecified")
        arrow_table = _extract_new_values_for_append(df_uri, arrow_table, filter_existing_joinids, context)

    def check_for_containment(
        df: pd.DataFrame,
        soma_df: DataFrame,
        ingestion_params: IngestionParams,
    ) -> bool:
        """For resume mode, check if the non-empty domain has already been written."""
        if ingestion_params.skip_existing_nonempty_domain:
            storage_ned = _read_nonempty_domain(soma_df)
            dim_range = ((int(df.index.min()), int(df.index.max())),)
            if _chunk_is_contained_in(dim_range, storage_ned):
                logging.log_io(
                    f"Skipped {df_uri}",
                    _util.format_elapsed(s, f"SKIPPED {df_uri}"),
                )
                return True
        return False

    if must_exist:
        # For update_obs, update_var, append_obs, and append_var, it's
        # an error situation if the dataframe doesn't already exist.
        soma_df = DataFrame.open(df_uri, "w", context=context)
        check_for_containment(df, soma_df, ingestion_params)

    else:
        # We could (and used to) do:
        #   if exists:
        #     open
        #   else:
        #     create
        # However, for remote object stores, that's two round-trip requests
        # to the server, whether the dataframe exists or not. Instead we
        # try create, doing the open if the create threw already-exists.
        # When the dataframe doesn't exist, this is just one round-trip request,
        # and when it does, it's two (as before).
        #
        # Note that for append/update, the dataframe must exist; but for
        # resume mode, the dataframe may or may not exist. (The point of
        # resume mode is to continue from a previous ingest which ended
        # prematurely.)
        try:
            soma_df = DataFrame.create(
                df_uri,
                schema=arrow_table.schema,
                # Note: tiledbsoma.io creates dataframes with soma_joinid being the one
                # and only index column.
                domain=[[0, shape - 1]],
                platform_config=platform_config,
                context=context,
            )
            # Save the original index name for outgest. We use JSON for elegant indication of index name
            # being None (in Python anyway).
            soma_df.metadata[SOMA_DATAFRAME_ORIGINAL_INDEX_NAME_JSON] = json.dumps(original_index_metadata)
        except (AlreadyExistsError, NotCreateableError) as e:
            if ingestion_params.error_if_already_exists:
                raise SOMAError(f"{df_uri} already exists") from e
            soma_df = DataFrame.open(df_uri, "w", context=context)
            check_for_containment(df, soma_df, ingestion_params)

    if ingestion_params.write_schema_no_data:
        logging.log_io(
            f"Wrote schema {df_uri}",
            _util.format_elapsed(s, f"FINISH WRITING SCHEMA {df_uri}"),
        )
        add_metadata(soma_df, additional_metadata)
        return soma_df

    tiledb_create_options = TileDBCreateOptions.from_platform_config(platform_config)
    tiledb_write_options = TileDBWriteOptions.from_platform_config(platform_config)

    if arrow_table and not any(column.num_chunks == 0 for column in arrow_table.columns):
        _write_arrow_table(arrow_table, soma_df, tiledb_create_options, tiledb_write_options)

    add_metadata(soma_df, additional_metadata)

    logging.log_io(
        f"Wrote   {df_uri}",
        _util.format_elapsed(s, f"FINISH WRITING {df_uri}"),
    )
    return soma_df


@deprecated(
    """This function is deprecated and will be removed in a future version of this package.

To add a new matrix as a layer within an existing SOMA Experiment (e.g., to X, obsm, varm), please use the more specific functions tiledbsoma.io.add_X_layer or tiledbsoma.io.add_matrix_to_collection. If you need to create a standalone SOMA NDArray outside of a pre-defined Experiment structure, please use the direct SOMA API constructors, such as tiledbsoma.SparseNDArray.create.""",
)
def create_from_matrix(
    cls: type[_NDArr],
    uri: str,
    matrix: Matrix | h5py.Dataset,
    platform_config: PlatformConfig | None = None,
    ingest_mode: IngestMode = "write",
    context: SOMATileDBContext | None = None,
) -> _NDArr:
    """Create and populate the ``soma_matrix`` from the contents of ``matrix``.

    This function is deprecated and will be removed in a future version of this package.

    To add a new matrix as a layer within an existing SOMA :class:`Experiment` (e.g., to X, obsm, varm), please use the more specific functions ``tiledbsoma.io.add_X_layer`` or ``tiledbsoma.io.add_matrix_to_collection``. If you need to create a standalone SOMA NDArray outside of a pre-defined :class:`Experiment` structure, please use the direct SOMA API constructors, such as ``tiledbsoma.SparseNDArray.create``.

    Lifecycle:
        Deprecated.
    """
    _check_for_deprecated_modes(ingest_mode)
    return _create_from_matrix(
        cls,
        uri,
        matrix,
        ingestion_params=IngestionParams(ingest_mode, None),
        platform_config=platform_config,
        context=context,
        axis_0_mapping=AxisIDMapping.identity(matrix.shape[0]),
        axis_1_mapping=AxisIDMapping.identity(matrix.shape[1]),
    )


def _create_from_matrix(
    cls: type[_NDArr],
    uri: str,
    matrix: Matrix | h5py.Dataset,
    *,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
    platform_config: PlatformConfig | None = None,
    context: SOMATileDBContext | None = None,
    axis_0_mapping: AxisIDMapping,
    axis_1_mapping: AxisIDMapping,
    must_exist: bool = False,
) -> _NDArr:
    """Internal helper for user-facing ``create_from_matrix``."""
    # Older SparseDataset has no ndim but it has a shape
    if len(matrix.shape) != 2:
        raise ValueError(f"expected matrix.shape == 2; got {matrix.shape}")

    s = _util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {uri}")

    if must_exist:
        if ingestion_params.error_if_already_exists:
            raise SOMAError(f"{uri} already exists")
        soma_ndarray = cls.open(uri, "w", platform_config=platform_config, context=context)
    else:
        try:
            shape: Sequence[int | None] = ()
            # A SparseNDArray must be appendable in soma.io.

            # Instead of
            #   shape = tuple(int(e) for e in matrix.shape)
            # we consult the registration mapping. This is important
            # in the case when multiple H5ADs/AnnDatas are being
            # ingested to an experiment which doesn't pre-exist.
            shape = (axis_0_mapping.get_shape(), axis_1_mapping.get_shape())

            soma_ndarray = cls.create(
                uri,
                type=pa.from_numpy_dtype(matrix.dtype),
                shape=shape,
                platform_config=platform_config,
                context=context,
            )
        except (AlreadyExistsError, NotCreateableError) as e:
            if ingestion_params.error_if_already_exists:
                raise SOMAError(f"{uri} already exists") from e
            soma_ndarray = cls.open(uri, "w", platform_config=platform_config, context=context)

    if ingestion_params.write_schema_no_data:
        logging.log_io(
            f"Wrote schema {uri}",
            _util.format_elapsed(s, f"FINISH WRITING SCHEMA {uri}"),
        )
        return soma_ndarray

    logging.log_io(
        f"Writing {uri}",
        _util.format_elapsed(s, f"START  WRITING {uri}"),
    )

    if isinstance(soma_ndarray, DenseNDArray):
        _write_matrix_to_denseNDArray(
            soma_ndarray,
            matrix,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(platform_config),
            tiledb_write_options=TileDBWriteOptions.from_platform_config(platform_config),
            ingestion_params=ingestion_params,
            additional_metadata=additional_metadata,
        )
    elif isinstance(soma_ndarray, SparseNDArray):  # SOMASparseNDArray
        _write_matrix_to_sparseNDArray(
            soma_ndarray,
            matrix,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(platform_config),
            tiledb_write_options=TileDBWriteOptions.from_platform_config(platform_config),
            ingestion_params=ingestion_params,
            additional_metadata=additional_metadata,
            axis_0_mapping=axis_0_mapping,
            axis_1_mapping=axis_1_mapping,
        )
    else:
        raise TypeError(f"unknown array type {type(soma_ndarray)}")

    logging.log_io(
        f"Wrote   {uri}",
        _util.format_elapsed(s, f"FINISH WRITING {uri}"),
    )
    return soma_ndarray


def update_obs(
    exp: Experiment,
    new_data: pd.DataFrame,
    *,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None = None,
    default_index_name: str = "obs_id",
) -> None:
    """Replaces the entire ``obs`` DataFrame with the contents of a new pandas DataFrame.

    This function is designed to perform a full replacement of the ``obs`` DataFrame. It assumes the input ``new_data``
    DataFrame represents the desired final state for the ``obs`` SOMA DataFrame. The operation implicitly relies on row
    order for alignment.

    **Details:**
    * **Schema Changes:** Columns present in ``new_data`` but not in the existing ``obs`` will be added. Columns present
        in the existing ``obs`` but absent from ``new_data`` will be dropped.
    * **Row Alignment:** The function requires that the input ```new_data`` DataFrame has the exact same number of rows
        and order as the existing ``obs`` DataFrame. It does *not* perform a join based on user-defined cell IDs or
        other index columns.

    To avoid data misalignment, the following workflow is recommended:

    1. Read the *entire* existing  ``obs`` DataFrame into memory.
    2. Perform all desired modifications (updating values, adding/dropping columns)
       on this DataFrame, preserving the original row order.
    3. Pass the fully modified DataFrame to ``update_obs``'s ``new_data`` argument.

    Args:
        exp: The :class:`SOMAExperiment` whose ``obs`` is to be updated. Must be opened for write.
        new_data: A pandas DataFrame containing the final desired data for the `obs` SOMA DataFrame.
        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.

        platform_config: Platform-specific options used to update this array, provided in the form
            ``{"tiledb": {"create": {"dataframe_dim_zstd_level": 7}}}``

        default_index_name: What to call the ``new_data`` index column if it is nameless in Pandas,
            or has name ``"index"``.

    Returns:
        None

    Lifecycle:
        Maturing.
    """
    _update_dataframe(
        exp.obs,
        new_data,
        "update_obs",
        context=context,
        platform_config=platform_config,
        default_index_name=default_index_name,
    )


def update_var(
    exp: Experiment,
    new_data: pd.DataFrame,
    measurement_name: str,
    *,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None = None,
    default_index_name: str = "var_id",
) -> None:
    """Replaces the entire ``var`` DataFrame with the contents of a new pandas DataFrame.

    This function is designed to perform a full replacement of the ``var`` DataFrame. It assumes the input ``new_data``
    DataFrame represents the desired final state for the ``var`` SOMA DataFrame. The operation implicitly relies on row
    order for alignment.

    **Details:**
    * **Schema Changes:** Columns present in ``new_data`` but not in the existing ``var`` will be added. Columns present
        in the existing ``var`` but absent from ``new_data`` will be dropped.
    * **Row Alignment:** The function requires that the input ```new_data`` DataFrame has the exact same number of rows
        and order as the existing ``var`` DataFrame. It does *not* perform a join based on user-defined cell IDs or
        other index columns.

    To avoid data misalignment, the following workflow is recommended:

    1. Read the *entire* existing  ``var`` DataFrame into memory.
    2. Perform all desired modifications (updating values, adding/dropping columns)
       on this DataFrame, preserving the original row order.
    3. Pass the fully modified DataFrame to ``update_var``'s ``new_data`` argument.

    Args:
        exp: The :class:`SOMAExperiment` whose ``var`` is to be updated. Must be opened for write.
        new_data: A pandas DataFrame containing the final desired data for the `var` SOMA DataFrame.
        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.

        platform_config: Platform-specific options used to update this array, provided in the form
            ``{"tiledb": {"create": {"dataframe_dim_zstd_level": 7}}}``

        default_index_name: What to call the ``new_data`` index column if it is nameless in Pandas,
            or has name ``"index"``.

    Returns:
        None

    Lifecycle:
        Maturing.
    """
    if measurement_name not in exp.ms:
        raise ValueError(f"cannot find measurement name {measurement_name} within experiment at {exp.uri}")

    _update_dataframe(
        exp.ms[measurement_name].var,
        new_data,
        "update_var",
        context=context,
        platform_config=platform_config,
        default_index_name=default_index_name,
    )


def _update_dataframe(
    sdf: DataFrame,
    new_data: pd.DataFrame,
    caller_name: str,
    *,
    context: SOMATileDBContext | None = None,
    platform_config: PlatformConfig | None,
    default_index_name: str,
) -> None:
    """See ``update_obs`` and ``update_var``. This is common helper code shared by both."""
    sdf.verify_open_for_writing()
    old_sig = conversions._string_dict_from_arrow_schema(sdf.schema)
    new_schema = conversions.df_to_arrow_schema(new_data, default_index_name)
    new_sig = conversions._string_dict_from_arrow_schema(new_schema)

    with DataFrame.open(sdf.uri, mode="r", context=context, platform_config=platform_config) as sdf_r:
        # Until we someday support deletes, this is the correct check on the existing,
        # contiguous soma join IDs compared to the new contiguous ones about to be created.
        old_jids = sorted(e.as_py() for e in sdf_r.read(column_names=["soma_joinid"]).concat()["soma_joinid"])
        num_old_data = len(old_jids)
        num_new_data = len(new_data)
        if num_old_data != num_new_data:
            raise ValueError(
                f"{caller_name}: old and new data must have the same row count; got {num_old_data} != {num_new_data}",
            )
        new_jids = list(range(num_new_data))
        if old_jids != new_jids:
            max_jid_diffs_to_display = 10
            jid_diffs = [(old_jid, new_jid) for old_jid, new_jid in zip(old_jids, new_jids) if old_jid != new_jid]
            jid_diff_strs = [f"{old_jid} != {new_jid}" for old_jid, new_jid in jid_diffs[:max_jid_diffs_to_display]] + (
                ["…"] if len(jid_diffs) > max_jid_diffs_to_display else []
            )
            raise ValueError(
                f"{caller_name}: old data soma_joinid must be [0,{num_old_data}), found {len(jid_diffs)} diffs: {', '.join(jid_diff_strs)}",
            )

        old_keys = set(old_sig.keys())
        new_keys = set(new_sig.keys())
        drop_keys = old_keys.difference(new_keys)
        add_keys = new_keys.difference(old_keys)
        common_keys = old_keys.intersection(new_keys)

        msgs = []
        for key in common_keys:
            old_type = old_sig[key]
            new_type = new_sig[key]

            if old_type != new_type:
                msgs.append(f"{key} type {old_type} != {new_type}")
        if msgs:
            msg = ", ".join(msgs)
            raise ValueError(f"unsupported type updates: {msg}")

        # Further operations are in-place for parsimony of memory usage:
        new_data = new_data.copy()
        arrow_table = conversions.df_to_arrow_table(new_data)
        arrow_schema = arrow_table.schema.remove_metadata()

        add_attrs = {}
        add_enmrs = {}
        for add_key in add_keys:
            # Don't directly use the new dataframe's dtypes. Go through the
            # to-Arrow-schema logic, and back, as this recapitulates the original
            # schema-creation logic.
            atype = arrow_schema.field(add_key).type
            if pa.types.is_dictionary(arrow_table.schema.field(add_key).type):
                add_attrs[add_key] = get_arrow_str_format(atype.index_type)

                enmr_format = get_arrow_str_format(atype.value_type)
                add_enmrs[add_key] = (enmr_format, atype.ordered)
            else:
                add_attrs[add_key] = get_arrow_str_format(atype)

        sdf_r._handle._handle._update_dataframe_schema(list(drop_keys), add_attrs, add_enmrs)

    _write_dataframe(
        df_uri=sdf.uri,
        df=new_data,
        id_column_name=default_index_name,
        ingestion_params=IngestionParams("update", label_mapping=None),
        context=context,
        platform_config=platform_config,
        axis_mapping=AxisIDMapping.identity(new_data.shape[0]),
        must_exist=True,
    )


def update_matrix(
    soma_ndarray: SparseNDArray | DenseNDArray,
    new_data: Matrix | h5py.Dataset,
    *,
    context: SOMATileDBContext | None = None,  # noqa: ARG001
    platform_config: PlatformConfig | None = None,
) -> None:
    """Given a ``SparseNDArray`` or ``DenseNDArray`` already opened for write,
    writes the new data. It is the caller's responsibility to ensure that the
    intended shape of written contents of the array match those of the existing
    data. The intended use-case is to replace updated numerical values.

    Example::

        with tiledbsoma.Experiment.open(uri, "w") as exp:
            tiledbsoma.io.update_matrix(
                exp.ms["RNA"].X["data"],
                adata.X,
            )

    Args:
        soma_ndarray: a ``SparseNDArray`` or ``DenseNDArray`` already opened for write.

        new_data: If the ``soma_ndarray`` is sparse, a Scipy CSR/CSC matrix or
            AnnData ``CSCDataset`` / ``CSRDataset``. If the ``soma_ndarray`` is dense,
            a NumPy NDArray.

        context: Optional :class:`SOMATileDBContext` containing storage parameters, etc.

        platform_config: Platform-specific options used to update this array, provided
            in the form ``{"tiledb": {"create": {"dataframe_dim_zstd_level": 7}}}``

    Returns:
        None

    Lifecycle:
        Maturing.
    """
    # More developer-level information on why we do not -- and cannot -- check
    # shape/bounding box:
    #
    # * The TileDB-SOMA "shape" can be huge x huge, with "room for growth" --
    #   this does not track the user-level "shape" and is not intended to.
    # * The TileDB-SOMA bounding box is, by contrast, intended to track the
    #   user-level "shape" but it is not thread-safe and may be incorrect.
    #   Please see
    #   https://github.com/single-cell-data/TileDB-SOMA/issues/1969
    #   https://github.com/single-cell-data/TileDB-SOMA/issues/1971

    s = _util.get_start_stamp()
    logging.log_io(
        f"Writing {soma_ndarray.uri}",
        f"START  UPDATING {soma_ndarray.uri}",
    )

    ingestion_params = IngestionParams("write", None)

    if isinstance(soma_ndarray, DenseNDArray):
        _write_matrix_to_denseNDArray(
            soma_ndarray,
            new_data,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(platform_config),
            tiledb_write_options=TileDBWriteOptions.from_platform_config(platform_config),
            ingestion_params=ingestion_params,
            additional_metadata=None,
        )
    elif isinstance(soma_ndarray, SparseNDArray):  # SOMASparseNDArray
        _write_matrix_to_sparseNDArray(
            soma_ndarray,
            new_data,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(platform_config),
            tiledb_write_options=TileDBWriteOptions.from_platform_config(platform_config),
            ingestion_params=ingestion_params,
            additional_metadata=None,
            axis_0_mapping=AxisIDMapping.identity(new_data.shape[0]),
            axis_1_mapping=AxisIDMapping.identity(new_data.shape[1]),
        )
    else:
        raise TypeError(f"unknown array type {type(soma_ndarray)}")

    logging.log_io(
        f"Wrote   {soma_ndarray.uri}",
        _util.format_elapsed(s, f"FINISH UPDATING {soma_ndarray.uri}"),
    )


def add_X_layer(
    exp: Experiment,
    measurement_name: str,
    X_layer_name: str,
    # E.g. a scipy.csr_matrix from scanpy analysis:
    X_layer_data: Matrix | h5py.Dataset,
    ingest_mode: IngestMode = "write",
    use_relative_uri: bool | None = None,
    context: SOMATileDBContext | None = None,
) -> None:
    """This is useful for adding X data, for example from
    `Scanpy <https://scanpy.readthedocs.io/>`_'s ``scanpy.pp.normalize_total``,
    ``scanpy.pp.log1p``, etc.

    Lifecycle:
        Maturing.
    """
    exp.verify_open_for_writing()
    _check_for_deprecated_modes(ingest_mode)
    add_matrix_to_collection(
        exp,
        measurement_name,
        collection_name="X",
        matrix_name=X_layer_name,
        matrix_data=X_layer_data,
        ingest_mode=ingest_mode,
        use_relative_uri=use_relative_uri,
        context=context,
    )


def add_matrix_to_collection(
    exp: Experiment,
    measurement_name: str,
    collection_name: str,
    matrix_name: str,
    # E.g. a scipy.csr_matrix from scanpy analysis:
    matrix_data: Matrix | h5py.Dataset,
    ingest_mode: IngestMode = "write",
    use_relative_uri: bool | None = None,
    context: SOMATileDBContext | None = None,
) -> None:
    """This is useful for adding X/obsp/varm/etc data, for example from
    Scanpy's ``scanpy.pp.normalize_total``, ``scanpy.pp.log1p``, etc.

    Lifecycle:
        Maturing.
    """
    ingestion_params = IngestionParams(ingest_mode, None)
    _check_for_deprecated_modes(ingest_mode)

    # For local disk and S3, creation and storage URIs are identical.  For
    # cloud, creation URIs look like tiledb://namespace/s3://bucket/path/to/obj
    # whereas storage URIs (for the same object) look like
    # tiledb://namespace/uuid.  When the caller passes a creation URI (which
    # they must) via exp.uri, we need to follow that.
    extend_creation_uri = exp.uri.startswith("tiledb://")

    with exp.ms[measurement_name] as meas:
        if extend_creation_uri:
            coll_uri = f"{exp.uri}/ms/{_util.sanitize_key(measurement_name)}/{_util.sanitize_key(collection_name)}"
        else:
            coll_uri = _util.uri_joinpath(meas.uri, _util.sanitize_key(collection_name))

        if collection_name in meas:
            coll = cast("Collection[RawHandle]", meas[collection_name])
        else:
            coll = _create_or_open_collection(
                Collection,
                coll_uri,
                ingestion_params=ingestion_params,
                context=context,
            )
            _maybe_set(meas, collection_name, coll, use_relative_uri=use_relative_uri)
        with coll:
            matrix_uri = _util.uri_joinpath(coll_uri, _util.sanitize_key(matrix_name))

            with _create_from_matrix(
                SparseNDArray,
                matrix_uri,
                matrix_data,
                ingestion_params=ingestion_params,
                context=context,
                axis_0_mapping=AxisIDMapping.identity(matrix_data.shape[0]),
                axis_1_mapping=AxisIDMapping.identity(matrix_data.shape[1]),
            ) as sparse_nd_array:
                _maybe_set(
                    coll,
                    matrix_name,
                    sparse_nd_array,
                    use_relative_uri=use_relative_uri,
                )


def _write_matrix_to_denseNDArray(
    soma_ndarray: DenseNDArray,
    matrix: Matrix | h5py.Dataset,
    tiledb_create_options: TileDBCreateOptions,
    tiledb_write_options: TileDBWriteOptions,  # noqa: ARG001
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
) -> None:
    """Write a matrix to an empty DenseNDArray."""
    add_metadata(soma_ndarray, additional_metadata)

    # TileDB does not support big-endian so coerce to little-endian
    if isinstance(matrix, np.ndarray) and matrix.dtype.byteorder == ">":
        matrix = matrix.byteswap().view(matrix.dtype.newbyteorder("<"))

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
        matrix_bounds = [(0, int(n - 1)) for n in matrix.shape]  # Cast for lint in case np.int64
        logging.log_io(
            None,
            f"Input bounds {tuple(matrix_bounds)} storage non-empty domain {storage_ned}",
        )
        if _chunk_is_contained_in(matrix_bounds, storage_ned):
            logging.log_io(f"Skipped {soma_ndarray.uri}", f"SKIPPED WRITING {soma_ndarray.uri}")
            return

    # Write all at once?
    if not tiledb_create_options.write_X_chunked:
        if not isinstance(matrix, np.ndarray):
            matrix = matrix.toarray()
        soma_ndarray.write((slice(None),), pa.Tensor.from_numpy(matrix))
        return

    # OR, write in chunks
    eta_tracker = eta.Tracker()
    if matrix.ndim == 2:
        nrow, ncol = matrix.shape
    elif matrix.ndim == 1:
        nrow = matrix.shape[0]
        ncol = 1
    else:
        raise ValueError(f"only 1D or 2D dense arrays are supported here; got {matrix.ndim}")

    # Number of rows to chunk by. These are dense writes, so this is loop-invariant.
    # * The goal_chunk_nnz is an older parameter. It's still important, as for backed AnnData,
    #   it controls how much is read into client RAM from the backing store on each chunk.
    # * The remote_cap_nbytes is an older parameter.
    # * Compute chunk sizes for both and take the minimum.
    chunk_size_using_nnz = math.ceil(tiledb_create_options.goal_chunk_nnz / ncol)

    try:
        # not scipy csr/csc
        itemsize = matrix.itemsize
    except AttributeError:
        # scipy csr/csc
        itemsize = matrix.data.itemsize

    total_nbytes = matrix.size * itemsize
    nbytes_num_chunks = math.ceil(total_nbytes / tiledb_create_options.remote_cap_nbytes)
    nbytes_num_chunks = min(1, nbytes_num_chunks)
    chunk_size_using_nbytes = math.floor(nrow / nbytes_num_chunks)

    chunk_size = min(chunk_size_using_nnz, chunk_size_using_nbytes)

    i = 0
    while i < nrow:
        t1 = time.time()
        i2 = i + chunk_size

        # Print doubly-inclusive lo..hi like 0..17 and 18..31.
        chunk_percent = min(100, 100 * (i2 - 1) / nrow)
        logging.log_io(
            None,
            "START  chunk rows %d..%d of %d (%.3f%%)" % (i, i2 - 1, nrow, chunk_percent),
        )

        chunk = matrix[i:i2, :] if matrix.ndim == 2 else matrix[i:i2]

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
                    "SKIP   chunk rows %d..%d of %d (%.3f%%)" % (i, i2 - 1, nrow, chunk_percent),
                )
                i = i2
                continue

        tensor = pa.Tensor.from_numpy(chunk) if isinstance(chunk, np.ndarray) else pa.Tensor.from_numpy(chunk.toarray())
        if matrix.ndim == 2:
            soma_ndarray.write((slice(i, i2 - 1), slice(0, ncol - 1)), tensor)
        else:
            soma_ndarray.write((slice(i, i2 - 1),), tensor)

        t2 = time.time()
        chunk_seconds = t2 - t1
        eta_seconds = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

        if chunk_percent < 100:
            logging.log_io(
                "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                "FINISH chunk in %.3f seconds, %7.3f%% done, ETA %s" % (chunk_seconds, chunk_percent, eta_seconds),
            )

        i = i2

    return


def _read_nonempty_domain(arr: SOMAArray) -> Any:  # noqa: ANN401
    try:
        return arr._handle.non_empty_domain()
    except (SOMAError, RuntimeError):
        # This means that we're open in write-only mode.
        # Reopen the array in read mode.
        # TODO macOS throws RuntimeError instead of SOMAError
        pass

    cls = type(arr)
    with cls.open(arr.uri, "r", platform_config=None, context=arr.context) as readarr:
        return readarr._handle.non_empty_domain()


def _find_sparse_chunk_size(
    matrix: SparseMatrix,
    start_index: int,
    axis: int,
    goal_chunk_nnz: int,
    mean_nnz: int,
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
        mean_nnz:
            Mean nnz along the desired axis. This is necessary for the backed
            case, and not needed in the non-backed case but easy to compute.
            It simplifies the internal API to simply pass it in unconditionally.
    """
    if isinstance(matrix, (sp.csr_matrix, sp.csc_matrix)):
        return _find_sparse_chunk_size_non_backed(
            matrix=matrix,
            start_index=start_index,
            axis=axis,
            goal_chunk_nnz=goal_chunk_nnz,
            mean_nnz=mean_nnz,
        )

    return _find_sparse_chunk_size_backed(
        matrix=matrix,
        start_index=start_index,
        axis=axis,
        goal_chunk_nnz=goal_chunk_nnz,
        mean_nnz=mean_nnz,
    )


def _find_sparse_chunk_size_non_backed(
    matrix: SparseMatrix,
    start_index: int,
    axis: int,
    goal_chunk_nnz: int,
    mean_nnz: int,  # noqa: ARG001
) -> int:
    """Helper routine for ``_find_sparse_chunk_size`` for when we're operating on AnnData
    matrices in non-backed mode. Here, unlike in backed mode, it's performant to exactly
    sum up nnz values from all matrix rows.
    """
    chunk_size = 0
    sum_nnz = 0
    coords: list[slice | int] = [slice(None), slice(None)]
    for index in range(start_index, matrix.shape[axis]):
        coords[axis] = index
        candidate_sum_nnz = sum_nnz + matrix[tuple(coords)].nnz
        if candidate_sum_nnz > goal_chunk_nnz:
            break
        sum_nnz = candidate_sum_nnz
        chunk_size += 1
        # The logger we use doesn't have a TRACE level. If it did, we'd use it here.
        # logging.logger.trace(
        #     f"non-backed: index={index} chunk_size={chunk_size} sum_nnz={sum_nnz} goal_nnz={goal_chunk_nnz}"
        # )
    return chunk_size


def _find_mean_nnz(matrix: Matrix, axis: int) -> int:
    """Helper for _find_sparse_chunk_size_backed."""
    extent = matrix.shape[axis]
    if extent == 0:
        return 0

    # Dense inputs don't have a .nnz
    if isinstance(matrix, (np.ndarray, h5py._hl.dataset.Dataset)):
        nr, nc = matrix.shape
        total = nr * nc
        return int(total // matrix.shape[axis])

    # This takes about as long but uses more RAM:
    #   total_nnz = matrix[:, :].nnz
    # So instead we break it up. Testing over a variety of H5AD sizes
    # shows that the performance is fine here.
    coords: list[slice] = [slice(None), slice(None)]  # type: ignore[unreachable]
    bsz = 1000
    total_nnz = 0
    for lo in range(0, extent, bsz):
        hi = min(extent, lo + bsz)
        coords[axis] = slice(lo, hi)
        total_nnz += matrix[tuple(coords)].nnz
    return math.ceil(total_nnz / extent)


def _find_sparse_chunk_size_backed(
    matrix: SparseMatrix,
    start_index: int,
    axis: int,
    goal_chunk_nnz: int,
    mean_nnz: int,
) -> int:
    """Helper routine for ``_find_sparse_chunk_size`` for when we're operating
    on AnnData matrices in backed mode.

    Empirically we find:

    * If the input matrix is sp.csr_matrix or sp.csc_matrix then getting all
      these nnz values is quick.  This happens when the input is AnnData via
      anndata.read_h5ad(name_of_h5ad) without the second backing-mode argument.

    * If the input matrix is ``anndata.abc.CSCDataset`` or
      ``anndata.abc.CSRDataset`` -- which happens with out-of-core anndata reads
      -- then getting all these nnz values is prohibitively expensive.  This
      happens when the input is AnnData via anndata.read_h5ad(name_of_h5ad, "r")
      with the second backing-mode argument, which is necessary for being able
      to ingest larger H5AD files.

    Say there are 100,000 rows, each with possibly quite different nnz values.
    Then in the non-backed case we simply check each row's nnz value. But for
    the backed case, that absolutely tanks performance.

    Another option is getting a sample of nnz values and using them as an
    estimator for nnz values of other rows. This often works, but there's enough
    diversity in the row nnzs that this estimator can result in us sending
    too-big chunks to remote services.

    It turns out, though, as a peculiarity of AnnData backed matrices, that
    while it's prohibitively expensive to ask the 10,000 questions
    matrix[0,:].nnz, matrix[1,:].nnz, ...  matrix[9999,:].nnz, it's quite quick
    to ask for matrix[0:10000,:].nnz.

    This is our way to thread the needle on good runtime performance with the
    AnnData backed-matrix API, while respecting remote resource limits:

    * Get the mean nnz for the entire matrix, along the desired axis.
      This needs to be computed only once, so we take it as an argument.

    * Set our initial estimate of chunk size to be the goal_chunk_nnz
      over the mean_nnz. E.g. if our goal is 100M nnz per chunk, and the matrix
      has average 4000 nnz per row, then we try chunk size to be
      100,000,000/4,000 = 25,000.

    * While tasking for the nnzs of each of those 25,000 rows is prohibitively
      expensive, we can quickly ask for the nnz of the contiguous region of
      25,000 rows.

    * Now, those 25,000 rows' sum nnz may or may not be reflective of the mean.
      If the guess is too high -- if those 25,000 rows' nnz is over the goal --
      we must reduce it; if it's too low, we must raise it. In either case we
      just divide by the overshoot/undershoot ratio. E.g. if the goal was 100M
      but the chunk nnz was 125M, then we divide 25,000 by 1.25 and try again.
      That second guess may not be quite right either, so we try a few times.
      Experiments across a variety of file sizes shows that we converge to
      between 70% and 100% of goal chunk nnz usually on the first try, and
      occasionally on the second and maybe the third. Again, _exact_ counting of
      each row's nnz is prohibitively expensive, so we adapt our algorithm to
      the constraints presented to us.

    Minor details on the basic theme:

    * We aim for between 70% and 100% of goal chunk nnz. Over is bad;
      under is less bad -- we do want to not make too many small requests of
      remote services, though.

    * On the very last bit of the matrix, there can be no sizing up --
      if there are only 7 rows remaining, that's that, end of story.
    """
    # Parameters as noted above and below.
    lower_ratio = 0.7
    upper_ratio = 1.0
    adjustment_ratio = 0.95
    max_iterations = 5

    if mean_nnz == 0:  # completely empty sparse array (corner case)
        return -1

    # This is num_rows or num_cols.
    extent = int(matrix.shape[axis])

    # Initial guess of chunk size.
    #
    # The goal_chunk_nnz is important, as for backed AnnData, it controls how much is read into
    # client RAM from the backing store on each chunk. We also subdivide chunks by
    # remote_cap_nbytes, if necessary, within _write_arrow_table in order to accommodate remote
    # object stores, which is a different ceiling.
    chunk_size = max(1, math.floor(goal_chunk_nnz / mean_nnz))
    if chunk_size > extent:
        chunk_size = extent

    # This is just matrix[0:m, :] or matrix[:, 0:m], but spelt out flexibly
    # given that the axis (0 or 1) is a variable.
    coords: list[slice] = [slice(None), slice(None)]
    coords[axis] = slice(start_index, start_index + chunk_size)
    chunk_nnz = int(matrix[tuple(coords)].nnz)

    # End of matrix, end of story.
    if start_index + chunk_size >= extent:
        return chunk_size

    # The entire chunk is unoccupied sparse; the caller can skip over it
    if chunk_nnz == 0:
        return chunk_size

    # Compare initial estimate chunk nnz against the goal
    ratio = chunk_nnz / goal_chunk_nnz

    # Close enough to full capacity without going over; not worth tweaking
    if ratio > lower_ratio and ratio <= upper_ratio:
        return chunk_size

    # This is a corner case, with artificially low goal_chunk_nnz far below
    # remote resource requirements. In this case, keep estimate to avoid the
    # iteration loop below.
    if chunk_size == 1:
        return chunk_size

    # Small chunk: try to enlarge for efficiency.
    # Large chunk: must split for local/remote resource limits.
    # In either case, we divide by the ration of actual vs goal.
    for _i in range(max_iterations):
        # Omitting the adjustment ratio would mean our next try would be likely
        # to be say 1.000002x the goal, resulting in an unnecessary extra
        # iteration.  Since the "just-right zone" is between 0.7 and 1.0 times
        # the goal, don't aim for the very edge of that zone: aim well within
        # it.
        chunk_size = max(1, int(adjustment_ratio * math.floor(chunk_size / ratio)))
        coords[axis] = slice(start_index, start_index + chunk_size)
        chunk_nnz = matrix[tuple(coords)].nnz
        ratio = chunk_nnz / goal_chunk_nnz
        if ratio > lower_ratio and ratio <= upper_ratio:
            return chunk_size

    # If the result is still smallish, so be it -- we have tried hard enough.
    # If the result is too large despite efforts, let remote resource-limit
    # errors be the hard fail, rather than trying to hard-enforce one of our
    # own.
    return chunk_size


def _write_matrix_to_sparseNDArray(
    soma_ndarray: SparseNDArray,
    matrix: Matrix,
    tiledb_create_options: TileDBCreateOptions,
    tiledb_write_options: TileDBWriteOptions,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata,
    axis_0_mapping: AxisIDMapping,
    axis_1_mapping: AxisIDMapping,
) -> None:
    """Write a matrix to an empty DenseNDArray."""
    # TileDB does not support big-endian so coerce to little-endian
    if isinstance(matrix, np.ndarray) and matrix.dtype.byteorder == ">":
        matrix = matrix.byteswap().view(matrix.dtype.newbyteorder("<"))

    def _coo_to_table(
        mat_coo: sp.coo_matrix,
        axis_0_mapping: AxisIDMapping,
        axis_1_mapping: AxisIDMapping,
        axis: int = 0,
        base: int = 0,
    ) -> pa.Table:
        soma_dim_0 = mat_coo.row + base if base > 0 and axis == 0 else mat_coo.row
        soma_dim_1 = mat_coo.col + base if base > 0 and axis == 1 else mat_coo.col

        # Apply registration mappings: e.g. columns 0,1,2,3 in an AnnData file might
        # have been assigned gene-ID labels 22,197,438,988.
        soma_dim_0 = axis_0_mapping.data[soma_dim_0]
        soma_dim_1 = axis_1_mapping.data[soma_dim_1]
        pydict = {
            "soma_data": mat_coo.data,
            "soma_dim_0": soma_dim_0,
            "soma_dim_1": soma_dim_1,
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
        matrix_bounds = [(0, int(n - 1)) for n in matrix.shape]  # Cast for lint in case np.int64
        logging.log_io(
            None,
            f"Input bounds {tuple(matrix_bounds)} storage non-empty domain {storage_ned}",
        )
        if _chunk_is_contained_in(matrix_bounds, storage_ned):
            logging.log_io(f"Skipped {soma_ndarray.uri}", f"SKIPPED WRITING {soma_ndarray.uri}")
            return

    add_metadata(soma_ndarray, additional_metadata)

    # Write all at once?
    if not tiledb_create_options.write_X_chunked:
        soma_ndarray.write(_coo_to_table(sp.coo_matrix(matrix), axis_0_mapping, axis_1_mapping))
        return

    # Or, write in chunks, striding across the most efficient slice axis

    stride_axis = 0
    if sp.isspmatrix_csc(matrix):
        # E.g. if we used anndata.X[:]
        stride_axis = 1
    if isinstance(matrix, CSCDataset):
        # E.g. if we used anndata.X without the [:]
        stride_axis = 1

    dim_max_size = matrix.shape[stride_axis]

    eta_tracker = eta.Tracker()
    goal_chunk_nnz = tiledb_create_options.goal_chunk_nnz
    mean_nnz = _find_mean_nnz(matrix, stride_axis)

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
            chunk_size = math.ceil(goal_chunk_nnz / matrix.shape[non_stride_axis])
        else:
            chunk_size = _find_sparse_chunk_size(  # type: ignore[unreachable]
                matrix,
                i,
                stride_axis,
                goal_chunk_nnz,
                mean_nnz,
            )
        if chunk_size == -1:  # completely empty array; nothing to write
            if i > 0:
                break
            chunk_size = 1

        # Don't display something like '0..100000 out of 98765' as this looks wrong.
        # Cap the part after the '..' at the dim_max_size.
        i2 = i + chunk_size
        if i2 > dim_max_size:
            i2 = dim_max_size

        coords[stride_axis] = slice(i, i2)
        chunk_coo = sp.coo_matrix(matrix[tuple(coords)])

        # As noted above, we support reading AnnData in backed mode which
        # provides opportunities as well as challenges.
        #
        # * Backed mode is crucial for our ability to ingest larger H5AD files
        #   -- those whose file sizes compete with host RAM.
        # * Semantics of AnnData backed-mode matrices is such that it is
        #   prohibitive (i.e. it ruins performance, and by a significant amount)
        #   to compute row nnz for every row -- which is precisely the
        #   information that any chunk-sizing algorithm requires.
        # * Therefore we use a sampling algorithm to estimate the chunk size
        #   (row-count) to satisfy a goal chunk nnz value. However, the nnz
        #   values in the sample do not always well represent the values in the
        #   population.
        #
        # Therefore as a performance optimization we do the following:
        #
        # * Use a sample of a small number of rows (100 as of this writing)
        #   to estimate a chunk size that satisfies the goal chunk nnz.
        # * Acquire the full chunk nnz (which is much more performant at the
        #   AnnData level than getting each of the individual row nnz values)
        # * Adapt that downward.
        #
        # In a future C++-only matrix-writer implementation where buffer sizes
        # are exposed to the implementation langauge -- a benefit we do not
        # enjoy here in Python -- it will be easier to simply fill buffers and
        # send them off, with simplified logic.
        num_tries = 0
        max_tries = 20
        while chunk_coo.nnz > tiledb_create_options.goal_chunk_nnz:
            num_tries += 1
            # The logger we use doesn't have a TRACE level. If it did, we'd use it here.
            # logging.logger.trace(
            #    f"Adapt: {num_tries}/{max_tries} {chunk_coo.nnz}/{tiledb_create_options.goal_chunk_nnz}"
            # )
            if num_tries > max_tries:
                raise SOMAError(
                    f"Unable to accommodate goal_chunk_nnz {goal_chunk_nnz}. "
                    "This may be reduced in TileDBCreateOptions.",
                )

            ratio = chunk_coo.nnz / tiledb_create_options.goal_chunk_nnz
            chunk_size = math.floor(0.9 * (i2 - i) / ratio)
            if chunk_size < 1:
                raise SOMAError(
                    f"Unable to accommodate a single row at goal_chunk_nnz {goal_chunk_nnz}. "
                    "This may be reduced in TileDBCreateOptions.",
                )
            i2 = i + chunk_size
            coords[stride_axis] = slice(i, i2)
            chunk_coo = sp.coo_matrix(matrix[tuple(coords)])

        chunk_percent = min(100, 100 * i2 / dim_max_size)

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
                    "SKIP   chunk rows %d..%d of %d (%.3f%%), nnz=%d, goal=%d"
                    % (
                        i,
                        i2 - 1,
                        dim_max_size,
                        chunk_percent,
                        chunk_coo.nnz,
                        tiledb_create_options.goal_chunk_nnz,
                    ),
                )
                i = i2
                continue

        # Print doubly inclusive lo..hi like 0..17 and 18..31.
        logging.log_io(
            None,
            "START  chunk rows %d..%d of %d (%.3f%%), nnz=%d, goal=%d"
            % (
                i,
                i2 - 1,
                dim_max_size,
                chunk_percent,
                chunk_coo.nnz,
                tiledb_create_options.goal_chunk_nnz,
            ),
        )

        arrow_table = _coo_to_table(chunk_coo, axis_0_mapping, axis_1_mapping, stride_axis, i)
        _write_arrow_table(arrow_table, soma_ndarray, tiledb_create_options, tiledb_write_options)

        t2 = time.time()
        chunk_seconds = t2 - t1
        eta_seconds = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

        if chunk_percent < 100:
            logging.log_io(
                "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                "FINISH chunk in %.3f seconds, %7.3f%% done, ETA %s" % (chunk_seconds, chunk_percent, eta_seconds),
            )

        i = i2


def _chunk_is_contained_in(
    chunk_bounds: Sequence[tuple[int, int]],
    storage_nonempty_domain: Sequence[tuple[int | None, int | None]],
) -> bool:
    """Determines if a dim range is included within the array's non-empty domain.  Ranges are inclusive
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
    if len(storage_nonempty_domain) == 0:
        return False

    if len(chunk_bounds) != len(storage_nonempty_domain):
        raise SOMAError(
            f"internal error: ingest data ndim {len(chunk_bounds)} != storage ndim {len(storage_nonempty_domain)}",
        )
    return all(_chunk_is_contained_in_axis(chunk_bounds, storage_nonempty_domain, i) for i in range(len(chunk_bounds)))


def _chunk_is_contained_in_axis(
    chunk_bounds: Sequence[tuple[int, int]],
    storage_nonempty_domain: Sequence[tuple[int | None, int | None]],
    stride_axis: int,
) -> bool:
    """Helper function for ``_chunk_is_contained_in``."""
    if len(storage_nonempty_domain) == 0:
        return False

    storage_lo, storage_hi = storage_nonempty_domain[stride_axis]
    if storage_lo is None or storage_hi is None:
        # E.g. an array has had its schema created but no data written yet
        return False

    chunk_lo, chunk_hi = chunk_bounds[stride_axis]
    if chunk_lo < storage_lo or chunk_lo > storage_hi:
        return False
    return not (chunk_hi < storage_lo or chunk_hi > storage_hi)


def _maybe_ingest_uns(
    m: Measurement,
    uns: UnsMapping,
    *,
    platform_config: PlatformConfig | None,
    context: SOMATileDBContext | None,
    ingestion_params: IngestionParams,
    use_relative_uri: bool | None,
    uns_keys: Sequence[str] | None = None,
    additional_metadata: AdditionalMetadata = None,
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
        uns_keys=uns_keys,
        additional_metadata=additional_metadata,
    )


def _ingest_uns_dict(
    parent: AnyTileDBCollection,
    parent_key: str,
    dct: UnsMapping,
    *,
    platform_config: PlatformConfig | None,
    context: SOMATileDBContext | None,
    ingestion_params: IngestionParams,
    use_relative_uri: bool | None,
    uns_keys: Sequence[str] | None = None,
    level: int = 0,
    additional_metadata: AdditionalMetadata = None,
) -> None:
    with _create_or_open_collection(
        Collection,
        _util.uri_joinpath(parent.uri, _util.sanitize_key(parent_key)),
        ingestion_params=ingestion_params,
        context=context,
        additional_metadata=additional_metadata,
    ) as coll:
        _maybe_set(parent, parent_key, coll, use_relative_uri=use_relative_uri)
        coll.metadata[_TILEDBSOMA_TYPE] = "uns"
        for key, value in dct.items():
            if level == 0 and uns_keys is not None and key not in uns_keys:
                continue
            _ingest_uns_node(
                coll,
                key,
                value,
                platform_config=platform_config,
                context=context,
                ingestion_params=ingestion_params,
                additional_metadata=additional_metadata,
                use_relative_uri=use_relative_uri,
                level=level + 1,
            )

    msg = f"Wrote   {coll.uri} (uns collection)"
    logging.log_io(msg, msg)


def _ingest_uns_node(
    coll: AnyTileDBCollection,
    key: str,
    value: UnsNode,
    *,
    platform_config: PlatformConfig | None,
    context: SOMATileDBContext | None,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
    use_relative_uri: bool | None,
    level: int,
) -> None:
    if isinstance(value, np.generic):
        # This is some kind of numpy scalar value. Metadata entries
        # only accept native Python types, so unwrap it.
        value = value.item()

    if isinstance(value, (int, float, str)):
        # Primitives get set on the metadata.
        coll.metadata[key] = value
        return

    ingest_platform_ctx: IngestPlatformCtx = {
        "platform_config": platform_config,
        "context": context,
        "ingestion_params": ingestion_params,
        "additional_metadata": additional_metadata,
    }
    if isinstance(value, Mapping):
        # Mappings are represented as sub-dictionaries.
        _ingest_uns_dict(
            coll,
            key,
            value,
            **ingest_platform_ctx,
            use_relative_uri=use_relative_uri,
            level=level + 1,
        )
        return

    if isinstance(value, pd.DataFrame):
        num_rows = value.shape[0]
        with _write_dataframe(
            _util.uri_joinpath(coll.uri, _util.sanitize_key(key)),
            # _write_dataframe modifies passed DataFrame in-place (adding a "soma_joinid" index)
            value.copy(),
            None,
            axis_mapping=AxisIDMapping.identity(num_rows),
            **ingest_platform_ctx,
        ) as df:
            _maybe_set(coll, key, df, use_relative_uri=use_relative_uri)
        return

    if isinstance(value, list) or "numpy" in str(type(value)):
        value = np.asarray(value)
    if isinstance(value, np.ndarray):
        _ingest_uns_array(
            coll,
            key,
            value,
            use_relative_uri=use_relative_uri,
            ingest_platform_ctx=ingest_platform_ctx,
        )
        return

    msg = f"Skipped {coll.uri}[{key!r}] (uns object): unrecognized type {type(value)}"
    logging.log_io(msg, msg)


def _ingest_uns_array(
    coll: AnyTileDBCollection,
    key: str,
    value: NPNDArray,
    use_relative_uri: bool | None,
    ingest_platform_ctx: IngestPlatformCtx,
) -> None:
    """Ingest an uns Numpy array.

    Delegates to :func:`_ingest_uns_string_array` for string arrays, and to
    :func:`_ingest_uns_ndarray` for numeric arrays.
    """
    if value.dtype.names is not None:
        # This is a structured array, which we do not support.
        logging.log_io_same(f"Skipped {coll.uri}[{key!r}] (uns): unsupported structured array")

    if value.dtype.char in ("U", "O"):
        # In the wild it's quite common to see arrays of strings in uns data.
        # Frequent example: uns["louvain_colors"].
        _ingest_uns_string_array(
            coll,
            key,
            value,
            use_relative_uri=use_relative_uri,
            **ingest_platform_ctx,
        )
    else:
        _ingest_uns_ndarray(
            coll,
            key,
            value,
            use_relative_uri=use_relative_uri,
            **ingest_platform_ctx,
        )


def _ingest_uns_string_array(
    coll: AnyTileDBCollection,
    key: str,
    value: NPNDArray,
    platform_config: PlatformConfig | None,
    context: SOMATileDBContext | None,
    *,
    use_relative_uri: bool | None,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
) -> None:
    """Ingest an uns string array. In the SOMA data model, we have NDArrays _of number only_ ...
    so we need to make this a SOMADataFrame.

    Ideally we don't want to add an index column "soma_joinid" -- "index", maybe.
    However, ``SOMADataFrame`` _requires_ that soma_joinid be present, either
    as an index column, or as a data column. The former is less confusing.
    """
    if len(value.shape) == 1:
        helper = _ingest_uns_1d_string_array
    elif len(value.shape) == 2:
        helper = _ingest_uns_2d_string_array
    else:
        msg = f"Skipped {coll.uri}[{key!r}] (uns object): string array is neither one-dimensional nor two-dimensional"
        logging.log_io(msg, msg)
        return

    helper(
        coll=coll,
        key=key,
        value=value,
        platform_config=platform_config,
        context=context,
        use_relative_uri=use_relative_uri,
        ingestion_params=ingestion_params,
        additional_metadata=additional_metadata,
    )


def _ingest_uns_1d_string_array(
    coll: AnyTileDBCollection,
    key: str,
    value: NPNDArray,
    platform_config: PlatformConfig | None,
    context: SOMATileDBContext | None,
    *,
    use_relative_uri: bool | None,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
) -> None:
    """Helper for ``_ingest_uns_string_array``."""
    n = len(value)
    # An array like ["a", "b", "c"] becomes a DataFrame like
    # soma_joinid value
    # 0           a
    # 1           b
    # 2           c
    df = pd.DataFrame(
        data={
            SOMA_JOINID: np.arange(n, dtype=np.int64),
            _UNS_OUTGEST_COLUMN_NAME_1D: [str(e) if e else "" for e in value],
        },
    )
    df.set_index("soma_joinid", inplace=True)

    df_uri = _util.uri_joinpath(coll.uri, _util.sanitize_key(key))
    with _write_dataframe_impl(
        df,
        df_uri,
        None,
        shape=df.shape[0],
        ingestion_params=ingestion_params,
        platform_config=platform_config,
        context=context,
        additional_metadata=additional_metadata,
    ) as soma_df:
        _maybe_set(coll, key, soma_df, use_relative_uri=use_relative_uri)
        # Since ND arrays in the SOMA data model are arrays _of number_,
        # we write these as dataframes. The metadata hint reminds us to
        # convert back to numpy.ndarray on outgest, since the numpy
        # data model supports ND arrays of string.
        soma_df.metadata[_UNS_OUTGEST_HINT_KEY] = _UNS_OUTGEST_HINT_1D


def _ingest_uns_2d_string_array(
    coll: AnyTileDBCollection,
    key: str,
    value: NPNDArray,
    platform_config: PlatformConfig | None,
    context: SOMATileDBContext | None,
    *,
    use_relative_uri: bool | None,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
) -> None:
    """Helper for ``_ingest_uns_string_array``. Even if the 2D array is 1xN or Nx1, we
    must nonetheless keep this as 2D rather than flattening to length-N 1D. That's because
    this ``uns`` data is solely of interest for AnnData ingest/outgest, and it must go
    back out the way it came in.
    """
    num_rows, num_cols = value.shape
    data: dict[str, Any] = {"soma_joinid": np.arange(num_rows, dtype=np.int64)}
    # An array like [["a", "b", "c"], ["d", "e", "f"]] becomes a DataFrame like
    # soma_joinid values_0 values_1 values_2
    # 0           a        b        c
    # 1           d        e        f
    for j in range(num_cols):
        column_name = f"values_{j}"
        data[column_name] = [str(e) if e else "" for e in value[:, j]]
    df = pd.DataFrame(data=data)
    df.set_index("soma_joinid", inplace=True)

    df_uri = _util.uri_joinpath(coll.uri, _util.sanitize_key(key))
    with _write_dataframe_impl(
        df,
        df_uri,
        None,
        shape=df.shape[0],
        ingestion_params=ingestion_params,
        additional_metadata=additional_metadata,
        platform_config=platform_config,
        context=context,
    ) as soma_df:
        _maybe_set(coll, key, soma_df, use_relative_uri=use_relative_uri)
        # Since ND arrays in the SOMA data model are arrays _of number_,
        # we write these as dataframes. The metadata hint reminds us to
        # convert back to numpy.ndarray on outgest, since the numpy
        # data model supports ND arrays of string.
        soma_df.metadata[_UNS_OUTGEST_HINT_KEY] = _UNS_OUTGEST_HINT_2D


def _ingest_uns_ndarray(
    coll: AnyTileDBCollection,
    key: str,
    value: NPNDArray,
    platform_config: PlatformConfig | None,
    context: SOMATileDBContext | None,
    *,
    use_relative_uri: bool | None,
    ingestion_params: IngestionParams,
    additional_metadata: AdditionalMetadata = None,
) -> None:
    arr_uri = _util.uri_joinpath(coll.uri, _util.sanitize_key(key))

    if any(e <= 0 for e in value.shape):
        msg = f"Skipped {arr_uri} (uns ndarray): zero in shape {value.shape}"
        logging.log_io(msg, msg)
        return

    try:
        pa_dtype = pa.from_numpy_dtype(value.dtype)
    except pa.ArrowNotImplementedError:
        msg = f"Skipped {arr_uri} (uns ndarray): unsupported dtype {value.dtype!r} ({value.dtype})"
        logging.log_io(msg, msg)
        return
    try:
        soma_arr = DenseNDArray.create(
            arr_uri,
            type=pa_dtype,
            shape=value.shape,
            platform_config=platform_config,
            context=context,
        )
    except (AlreadyExistsError, NotCreateableError):
        soma_arr = DenseNDArray.open(arr_uri, "w", context=context)

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

        _write_matrix_to_denseNDArray(
            soma_arr,
            value,
            tiledb_create_options=TileDBCreateOptions.from_platform_config(platform_config),
            tiledb_write_options=TileDBWriteOptions.from_platform_config(platform_config),
            ingestion_params=ingestion_params,
            additional_metadata=additional_metadata,
        )

    msg = f"Wrote   {soma_arr.uri} (uns ndarray)"
    logging.log_io(msg, msg)


def _concurrency_level(context: SOMATileDBContext) -> int:
    """Private helper function to determine appropriate concurrency level for
    ingestion of H5AD when use_multiprocessing is enabled.

    Functionally, this just allows the user to control concurrency via the
    context configuration ``soma.compute_concurrency_level``, with error checking.
    """
    concurrency_level: int = os.cpu_count() or 1
    if context is not None:
        concurrency_level = min(
            concurrency_level,
            int(context.tiledb_config.get("soma.compute_concurrency_level", concurrency_level)),
        )
    return concurrency_level


@deprecated(
    """The 'resume' ingest_mode is deprecated and will be removed in a future version. The current implementation has a known issue that can can cause multi-dataset appends to not resume correctly.

The recommended and safest approach for recovering from a failed ingestion is to delete the partially written SOMA Experiment and restart the ingestion process from the original input files or a known-good backup.""",
    stacklevel=3,
)
def _resume_mode_is_deprecated() -> None:
    pass


def _check_for_deprecated_modes(ingest_mode: str) -> None:
    if ingest_mode == "resume":
        _resume_mode_is_deprecated()
