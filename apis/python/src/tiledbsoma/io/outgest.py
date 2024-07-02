# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Outgestion methods.

This module contains methods to export SOMA artifacts to other formats.
Currently only ``.h5ad`` (`AnnData <https://anndata.readthedocs.io/>`_) is supported.
"""

import json
from typing import (
    Any,
    Dict,
    KeysView,
    Optional,
    Sequence,
    Union,
    cast,
)

import anndata as ad
import numpy as np
import pandas as pd

from .. import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    Measurement,
    SparseNDArray,
    _util,
    logging,
)
from .._constants import SOMA_JOINID
from .._exception import SOMAError
from .._types import NPNDArray, Path
from . import conversions
from ._common import (
    _DATAFRAME_ORIGINAL_INDEX_NAME_JSON,
    _UNS_OUTGEST_COLUMN_NAME_1D,
    _UNS_OUTGEST_COLUMN_PREFIX_2D,
    _UNS_OUTGEST_HINT_1D,
    _UNS_OUTGEST_HINT_2D,
    _UNS_OUTGEST_HINT_KEY,
    Matrix,
)


# ----------------------------------------------------------------
def to_h5ad(
    experiment: Experiment,
    h5ad_path: Path,
    measurement_name: str,
    *,
    X_layer_name: Optional[str] = "data",
    obs_id_name: Optional[str] = None,
    var_id_name: Optional[str] = None,
    obsm_varm_width_hints: Optional[Dict[str, Dict[str, int]]] = None,
    uns_keys: Optional[Sequence[str]] = None,
) -> None:
    """Converts the experiment group to `AnnData <https://anndata.readthedocs.io/>`_
    format and writes it to the specified ``.h5ad`` file.

    Arguments are as in ``to_anndata``.

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
        obsm_varm_width_hints=obsm_varm_width_hints,
        uns_keys=uns_keys,
    )

    s2 = _util.get_start_stamp()
    logging.log_io(None, f"START  write {h5ad_path}")

    anndata.write_h5ad(h5ad_path)

    logging.log_io(None, _util.format_elapsed(s2, f"FINISH write {h5ad_path}"))

    logging.log_io(
        None, _util.format_elapsed(s, f"FINISH Experiment.to_h5ad -> {h5ad_path}")
    )


# ----------------------------------------------------------------
def _extract_X_key(
    measurement: Measurement,
    X_layer_name: str,
    nobs: int,
    nvar: int,
) -> Matrix:
    """Helper function for to_anndata"""

    if X_layer_name not in measurement.X:
        raise ValueError(
            f"X_layer_name {X_layer_name} not found in data: {measurement.X.keys()}"
        )

    # Acquire handle to TileDB-SOMA data
    soma_X_data_handle = measurement.X[X_layer_name]

    # Read data from SOMA into memory
    if isinstance(soma_X_data_handle, DenseNDArray):
        data = soma_X_data_handle.read((slice(None), slice(None))).to_numpy()
    elif isinstance(soma_X_data_handle, SparseNDArray):
        X_mat = soma_X_data_handle.read().tables().concat().to_pandas()
        data = conversions.csr_from_tiledb_df(X_mat, nobs, nvar)
    else:
        raise TypeError(f"Unexpected NDArray type {type(soma_X_data_handle)}")

    return data


# ----------------------------------------------------------------
def to_anndata(
    experiment: Experiment,
    measurement_name: str,
    *,
    X_layer_name: Optional[str] = "data",
    extra_X_layer_names: Optional[Union[Sequence[str], KeysView[str]]] = None,
    obs_id_name: Optional[str] = None,
    var_id_name: Optional[str] = None,
    obsm_varm_width_hints: Optional[Dict[str, Dict[str, int]]] = None,
    uns_keys: Optional[Sequence[str]] = None,
) -> ad.AnnData:
    """Converts the experiment group to `AnnData <https://anndata.readthedocs.io/>`_
    format. Choice of matrix formats is following what we often see in input
    ``.h5ad`` files:

    * ``X`` as ``scipy.sparse.csr_matrix``
    * ``obs``,``var`` as ``pandas.dataframe``
    * ``obsm``,``varm`` arrays as ``numpy.ndarray``
    * ``obsp``,``varp`` arrays as ``scipy.sparse.csr_matrix``

    The ``X_layer_name`` is the name of the TileDB-SOMA measurement's ``X``
    collection which will be outgested to the resulting AnnData object's
    ``adata.X``. If this is ``None``, then the return value's ``adata.X`` will
    be None, and ``adata.layers`` will be unpopulated. If this is not ``None``,
    then ``adata.X`` will be taken from this layer name within the input
    measurement.

    The ``extra_X_layer_names`` are used to specify how the output ``adata``
    object's ``adata.layers`` is populated.  The default behavior --
    ``extra_X_layer_names`` being ``None`` -- means that ``adata.layers`` will be
    empty.  If ``extra_X_layer_names`` is a provided list these will be used for
    populating ``adata.layers``. If you want all the layers to be outgested,
    without having to name them individually, you can use
    ``extra_X_layer_names=experiment.ms[measurement_name].X.keys()``.  To make
    this low-friction for you, we introduce one more feature: we'll ignore
    ``X_layer_name`` when populating ``adata.layers``.  For example, if X keys
    are ``"a"``, ``"b"``, ``"c"``, ``"d"``, and you say ``X_layer_name="b"`` and
    ``extra_X_layer_names=experiment.ms[measurement_name].X.keys()``, we'll not
    write ``"b"`` to ``adata.layers``.

    The ``obs_id_name`` and ``var_id_name`` are columns within the TileDB-SOMA
    experiment which will become index names within the resulting AnnData
    object's ``obs``/``var`` dataframes. If not specified as arguments, the
    TileDB-SOMA's dataframes will be checked for an original-index-name key.
    When that also is unavailable, these default to ``"obs_id"`` and
    ``"var_id"``, respectively.

    The ``obsm_varm_width_hints`` is optional. If provided, it should be of the form
    ``{"obsm":{"X_tSNE":2}}`` to aid with export errors.

    If ``uns_keys`` is provided, only the specified top-level ``uns`` keys
    are extracted.  The default is to extract them all.  Use ``uns_keys=[]``
    to not ingest any ``uns`` keys.

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

    # How to choose index name for AnnData obs and var dataframes:
    # * If the desired names are passed in, use them.
    # * Else if the names used at ingest time are available, use them.
    # * Else use the default/fallback name.

    # Restore the original index name for outgest. We use JSON for elegant indication of index
    # name being None (in Python anyway). It may be 'null' which maps to Pyhton None.
    obs_id_name = obs_id_name or json.loads(
        experiment.obs.metadata.get(_DATAFRAME_ORIGINAL_INDEX_NAME_JSON, '"obs_id"')
    )
    var_id_name = var_id_name or json.loads(
        measurement.var.metadata.get(_DATAFRAME_ORIGINAL_INDEX_NAME_JSON, '"var_id"')
    )

    obs_df = experiment.obs.read().concat().to_pandas()
    obs_df.drop([SOMA_JOINID], axis=1, inplace=True)
    if obs_id_name is not None:
        if obs_id_name not in obs_df.keys():
            raise ValueError(
                f"requested obs IDs column name {obs_id_name} not found in input: {obs_df.keys()}"
            )
        obs_df.set_index(obs_id_name, inplace=True)
    else:
        # There are multiple cases to be handled here, all tested in CI.
        # This else-block handle this one:
        #
        #                 orig.ident  nCount_RNA  ...
        # ATGCCAGAACGACT           0        70.0  ...
        # CATGGCCTGTGCAT           0        85.0  ...
        # GAACCTGATGAACC           0        87.0  ...
        #
        # Namely:
        # * The input AnnData dataframe had an index with no name
        # * In the SOMA experiment we name that column "obs_id" and our index is "soma_joinid"
        # * On outgest we drop "soma_joinid"
        # * The thing we named "obs_id" needs to become the index again ...
        # * ... and it needs to be nameless.
        if "obs_id" in obs_df:
            obs_df.set_index("obs_id", inplace=True)
            obs_df.index.name = None

    var_df = measurement.var.read().concat().to_pandas()

    var_df.drop([SOMA_JOINID], axis=1, inplace=True)
    if var_id_name is not None:
        if var_id_name not in var_df.keys():
            raise ValueError(
                f"requested var IDs column name {var_id_name} not found in input: {var_df.keys()}"
            )
        var_df.set_index(var_id_name, inplace=True)
    else:
        if "var_id" in var_df:
            var_df.set_index("var_id", inplace=True)
            var_df.index.name = None

    nobs = len(obs_df.index)
    nvar = len(var_df.index)

    anndata_X = None
    anndata_X_dtype = None  # some datasets have no X
    anndata_layers = {}

    # Let them use
    #   extra_X_layer_names=exp.ms["RNA"].X.keys()
    # while avoiding
    #   TypeError: 'ABCMeta' object is not subscriptable
    if isinstance(extra_X_layer_names, KeysView):
        extra_X_layer_names = list(extra_X_layer_names)

    if X_layer_name is None and extra_X_layer_names:
        # The latter boolean check covers both not None and not []
        raise ValueError(
            "If X_layer_name is None, extra_X_layer_names must not be provided"
        )

    if X_layer_name is not None:
        anndata_X = _extract_X_key(measurement, X_layer_name, nobs, nvar)
        anndata_X_dtype = anndata_X.dtype

    if extra_X_layer_names is not None:
        for extra_X_layer_name in extra_X_layer_names:
            if extra_X_layer_name == X_layer_name:
                continue
            assert extra_X_layer_name is not None  # appease linter; already checked
            data = _extract_X_key(measurement, extra_X_layer_name, nobs, nvar)
            anndata_layers[extra_X_layer_name] = data

    if obsm_varm_width_hints is None:
        obsm_varm_width_hints = {}

    obsm = {}
    if "obsm" in measurement:
        obsm_width_hints = obsm_varm_width_hints.get("obsm", {})
        for key in measurement.obsm.keys():
            obsm[key] = _extract_obsm_or_varm(
                measurement.obsm[key], "obsm", key, nobs, obsm_width_hints
            )

    varm = {}
    if "varm" in measurement:
        varm_width_hints = obsm_varm_width_hints.get("obsm", {})
        for key in measurement.varm.keys():
            varm[key] = _extract_obsm_or_varm(
                measurement.varm[key], "varm", key, nvar, varm_width_hints
            )

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

    uns = {}
    if "uns" in measurement:
        s = _util.get_start_stamp()
        logging.log_io(None, f'Start  writing uns for {measurement["uns"].uri}')
        uns = _extract_uns(
            cast(Collection[Any], measurement["uns"]),
            uns_keys=uns_keys,
        )
        logging.log_io(
            None,
            _util.format_elapsed(s, f'Finish writing uns for {measurement["uns"].uri}'),
        )

    anndata = ad.AnnData(
        X=anndata_X,
        layers=anndata_layers,
        obs=obs_df,
        var=var_df,
        obsm=obsm,
        varm=varm,
        obsp=obsp,
        varp=varp,
        uns=uns,
        dtype=anndata_X_dtype,
    )

    logging.log_io(None, _util.format_elapsed(s, "FINISH Experiment.to_anndata"))

    return anndata


def _extract_obsm_or_varm(
    soma_nd_array: Union[SparseNDArray, DenseNDArray],
    collection_name: str,
    element_name: str,
    num_rows: int,
    width_configs: Dict[str, int],
) -> Matrix:
    """
    This is a helper function for ``to_anndata`` of ``obsm`` and ``varm`` elements.
    """

    # SOMA shape is capacity/domain -- not what AnnData wants.
    # But here do check the array is truly 2D.
    shape = soma_nd_array.shape
    if len(shape) != 2:
        raise ValueError(f"expected shape == 2; got {shape}")

    if isinstance(soma_nd_array, DenseNDArray):
        matrix = soma_nd_array.read().to_numpy()
        # The spelling ``sp.csr_array`` is more idiomatic but doesn't exist until Python
        # 3.8 and we still support Python 3.7
        return matrix

    matrix = soma_nd_array.read().tables().concat().to_pandas()

    # Problem to solve: whereas for other sparse arrays we have:
    #
    # * X    matrices are always nobs x nvar
    # * obsp matrices are always nobs x nobs
    # * varp matrices are always nvar x nvar
    #
    # but:
    #
    # * obsm is nobs x some number -- number of PCA components, etc.
    # * varm is nvar x some number -- number of PCs, etc.
    #
    # Three ways to get the number of columns for obsm/varm sparse matrices:
    #
    # * Explicit user specification
    # * Bounding-box metadata, if present
    # * Try arithmetic on nnz / num_rows, for densely occupied sparse matrices
    # Beyond that, we have to throw.

    description = f'{collection_name}["{element_name}"]'

    num_cols = width_configs.get(element_name, None)

    if num_cols is None:
        try:
            used_shape = soma_nd_array.used_shape()
            num_cols = used_shape[1][1] + 1
        except SOMAError:
            pass  # We tried; moving on to next option

    if num_cols is None:
        num_rows_times_width, coo_column_count = matrix.shape

        if coo_column_count != 3:
            raise SOMAError(
                f"internal error: expect COO width of 3; got {coo_column_count} for {description}"
            )

        if num_rows_times_width % num_rows == 0:
            num_cols = num_rows_times_width // coo_column_count

    if num_cols is None:
        raise SOMAError(
            f"could not determine outgest width for {description}: please try to_anndata's obsm_varm_width_hints option"
        )

    return conversions.csr_from_tiledb_df(matrix, num_rows, num_cols).toarray()


def _extract_uns(
    collection: Collection[Any],
    uns_keys: Optional[Sequence[str]] = None,
    level: int = 0,
) -> Dict[str, Any]:
    """
    This is a helper function for ``to_anndata`` of ``uns`` elements.
    """

    extracted: Dict[str, Any] = {}
    for key, element in collection.items():
        if level == 0 and uns_keys is not None and key not in uns_keys:
            continue

        if isinstance(element, Collection):
            extracted[key] = _extract_uns(element, level=level + 1)
        elif isinstance(element, DataFrame):
            hint = element.metadata.get(_UNS_OUTGEST_HINT_KEY)
            pdf = element.read().concat().to_pandas()
            if hint is None:
                extracted[key] = pdf
            elif hint == _UNS_OUTGEST_HINT_1D:
                extracted[key] = _outgest_uns_1d_string_array(pdf, element.uri)
            elif hint == _UNS_OUTGEST_HINT_2D:
                extracted[key] = _outgest_uns_2d_string_array(pdf, element.uri)
            else:
                msg = (
                    f"Warning: uns {collection.uri}[{key!r}] has "
                    + "{_UNS_OUTGEST_HINT_KEY} as unrecognized {hint}: leaving this as Pandas DataFrame"
                )
                logging.log_io_same(msg)
                extracted[key] = pdf
        elif isinstance(element, SparseNDArray):
            extracted[key] = element.read().tables().concat().to_pandas()
        elif isinstance(element, DenseNDArray):
            extracted[key] = element.read().to_numpy()
        else:
            logging.log_io_same(
                f"Skipping uns key {key} with unhandled type {element.soma_type}"
            )

    # Primitives got set on the SOMA-experiment uns metadata.
    for key, value in collection.metadata.items():
        if level == 0 and uns_keys is not None and key not in uns_keys:
            continue
        if not key.startswith("soma_"):
            extracted[key] = value

    return extracted


def _outgest_uns_1d_string_array(pdf: pd.DataFrame, uri_for_logging: str) -> NPNDArray:
    """Helper methods for _extract_uns"""
    num_rows, num_cols = pdf.shape
    # An array like ["a", "b", "c"] had become a DataFrame like
    # soma_joinid value
    # 0           a
    # 1           b
    if num_cols != 2:
        raise SOMAError(f"Expected 2 columns in {uri_for_logging}; got {num_cols}")
    for column_name in [SOMA_JOINID, _UNS_OUTGEST_COLUMN_NAME_1D]:
        if column_name not in pdf:
            raise SOMAError(f"Expected {column_name} column in {uri_for_logging}")
    return np.asarray(list(pdf[_UNS_OUTGEST_COLUMN_NAME_1D]))


def _outgest_uns_2d_string_array(pdf: pd.DataFrame, uri_for_logging: str) -> NPNDArray:
    """Helper methods for _extract_uns"""
    num_rows, num_cols = pdf.shape
    if num_cols < 2:
        raise SOMAError(f"Expected 2 columns in {uri_for_logging}; got {num_cols}")
    if SOMA_JOINID not in pdf:
        raise SOMAError(f"Expected {SOMA_JOINID} column in {uri_for_logging}")
    num_cols -= 1
    columns = []
    # An array like [["a", "b", "c"], ["d", "e", "f"]] had become a DataFrame like
    # soma_joinid values_0 values_1 values_2
    # 0           a        b        c
    # 1           d        e        f
    for i in range(num_cols):
        column_name = _UNS_OUTGEST_COLUMN_PREFIX_2D + str(i)
        if column_name not in pdf:
            raise SOMAError(f"Expected {column_name} column in {uri_for_logging}")
        columns.append(list(pdf[column_name]))
    return np.asarray(columns).transpose()
