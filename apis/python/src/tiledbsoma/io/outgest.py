# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Outgestion methods.

This module contains methods to export SOMA artifacts to other formats.
Currently only ``.h5ad`` (`AnnData <https://anndata.readthedocs.io/>`_) is supported.
"""

from __future__ import annotations

import json
from concurrent.futures import Future
from typing import (
    Any,
    KeysView,
    Sequence,
    Union,
    cast,
)

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp

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
from .._util import MISSING, Sentinel, _resolve_futures
from . import conversions
from ._common import (
    _DATAFRAME_ORIGINAL_INDEX_NAME_JSON,
    _UNS_OUTGEST_COLUMN_NAME_1D,
    _UNS_OUTGEST_COLUMN_PREFIX_2D,
    _UNS_OUTGEST_HINT_1D,
    _UNS_OUTGEST_HINT_2D,
    _UNS_OUTGEST_HINT_KEY,
    Matrix,
    UnsDict,
    UnsLeaf,
)

FutureUnsLeaf = Union[UnsLeaf, Future[UnsLeaf]]
FutureUnsDictNode = Union[FutureUnsLeaf, dict[str, "FutureUnsDictNode"]]


# ----------------------------------------------------------------
def to_h5ad(
    experiment: Experiment,
    h5ad_path: Path,
    measurement_name: str,
    *,
    X_layer_name: str | None = "data",
    obs_id_name: str | None = None,
    var_id_name: str | None = None,
    obsm_varm_width_hints: dict[str, dict[str, int]] | None = None,
    uns_keys: Sequence[str] | None = None,
) -> None:
    """Converts the experiment group to AnnData format and writes it to the specified ``.h5ad`` file.

    Arguments are as in ``to_anndata``.

    Lifecycle:
        Maturing.
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


def _read_partitioned_sparse(X: SparseNDArray, d0_size: int) -> pa.Table:
    # Partition dimension 0 based on a target point count and average
    # density of matrix. Magic number determined empirically, as a tradeoff
    # between concurrency and fixed query overhead.
    tgt_point_count = 96 * 1024**2
    nnz = X.nnz
    partition_sz = (
        max(1024 * round(d0_size * tgt_point_count / nnz / 1024), 1024)
        if nnz > 0
        else d0_size
    )
    partitions = [
        slice(st, min(st + partition_sz - 1, d0_size - 1))
        for st in range(0, d0_size, partition_sz)
    ]
    n_partitions = len(partitions)

    def _read_sparse_X(A: SparseNDArray, row_slc: slice) -> pa.Table:
        return A.read(coords=(row_slc,)).tables().concat()

    if n_partitions > 1:  # don't consume threads unless there is a reason to do so
        return pa.concat_tables(
            X.context.threadpool.map(_read_sparse_X, (X,) * n_partitions, partitions)
        )
    else:
        return _read_sparse_X(X, partitions[0])


def _extract_X_key(
    measurement: Measurement, X_layer_name: str, nobs: int, nvar: int
) -> Future[Matrix]:
    """Helper function for to_anndata"""
    if X_layer_name not in measurement.X:
        raise ValueError(
            f"X_layer_name {X_layer_name} not found in data: {measurement.X.keys()}"
        )

    # Acquire handle to TileDB-SOMA data
    X = measurement.X[X_layer_name]
    tp = X.context.threadpool

    # Read data from SOMA into memory
    if isinstance(X, DenseNDArray):

        def _read_dense_X(A: DenseNDArray) -> Matrix:
            return A.read((slice(None), slice(None))).to_numpy()

        return tp.submit(_read_dense_X, X)

    elif isinstance(X, SparseNDArray):

        def _read_X_partitions() -> Matrix:
            stk_of_coo = _read_partitioned_sparse(X, nobs)
            return conversions.csr_from_coo_table(
                stk_of_coo, nobs, nvar, context=X.context
            )

        return X.context.threadpool.submit(_read_X_partitions)

    else:
        raise TypeError(f"Unexpected NDArray type {type(X)}")


def _read_dataframe(
    df: DataFrame,
    default_index_name: str | None = None,
    fallback_index_name: str | None = None,
) -> pd.DataFrame:
    """Outgest a SOMA DataFrame to Pandas, including restoring the original index{,.name}.

    An `OriginalIndexMetadata` attached to the DataFrame, if present, contains a string (or `null`), indicating a
    SOMADataFrame column which should be restored as the pd.DataFrame index.

    `default_index_name`, if provided, overrides the stored `OriginalIndexMetadata`; a column with this name will be
    verified to exist, and set as index of the returned `pd.DataFrame`.

    If neither `default_index_name` nor `OriginalIndexMetadata` are provided, the `fallback_index_name` will be used.
    `to_anndata` passes "obs_id" / "var_id" for obs/var, matching `from_anndata`'s default `{obs,var}_id_name` values.

    NOTE: several edge cases result in the outgested DataFrame not matching the original DataFrame; see
    `test_dataframe_io_roundtrips.py` / https://github.com/single-cell-data/TileDB-SOMA/issues/2829.
    """
    # Read and validate the "original index metadata" stored alongside this SOMA DataFrame.
    original_index_metadata = json.loads(
        df.metadata.get(_DATAFRAME_ORIGINAL_INDEX_NAME_JSON, "null")
    )
    if not (
        original_index_metadata is None or isinstance(original_index_metadata, str)
    ):
        raise ValueError(
            f"{df.uri}: invalid {_DATAFRAME_ORIGINAL_INDEX_NAME_JSON} metadata: {original_index_metadata}"
        )

    pdf: pd.DataFrame = df.read().concat().to_pandas()
    # SOMA DataFrames always have a `soma_joinid` added, as part of the ingest process, which we remove on outgest.
    pdf.drop(columns=SOMA_JOINID, inplace=True)

    default_index_name = default_index_name or original_index_metadata
    if default_index_name is not None:
        # One or both of the following was true:
        # - Original DataFrame had an index name (other than "index") ⇒ that name was written as `OriginalIndexMetadata`
        # - `default_index_name` was provided (e.g. `{obs,var}_id_name` args to `to_anndata`)
        #
        # ⇒ Verify a column with that name exists, and set it as index (keeping its name).
        if default_index_name not in pdf.keys():
            raise ValueError(
                f"Requested ID column name {default_index_name} not found in input: {pdf.keys()}"
            )
        pdf.set_index(default_index_name, inplace=True)
    else:
        # The assumption here is that the original index was unnamed, and was given a "fallback name" (e.g. "obs_id",
        # "var_id") during ingest that matches the `fallback_index_name` arg here. In this case, we restore that column
        # as index, and remove the name.
        #
        # NOTE: several edge cases result in the outgested DF not matching the original DF; see
        # https://github.com/single-cell-data/TileDB-SOMA/issues/2829.
        if fallback_index_name is not None and fallback_index_name in pdf:
            pdf.set_index(fallback_index_name, inplace=True)
            pdf.index.name = None

    return pdf


def to_anndata(
    experiment: Experiment,
    measurement_name: str,
    *,
    X_layer_name: str | Sentinel | None = MISSING,
    extra_X_layer_names: Sequence[str] | KeysView[str] | None = None,
    obs_id_name: str | None = None,
    var_id_name: str | None = None,
    obsm_varm_width_hints: dict[str, dict[str, int]] | None = None,
    uns_keys: Sequence[str] | None = None,
) -> ad.AnnData:
    """Converts the experiment group to AnnData format.

    The choice of matrix formats is following what we often see in input ``.h5ad`` files:

    * ``X`` as ``scipy.sparse.csr_matrix``
    * ``obs``, ``var`` as ``pandas.dataframe``
    * ``obsm``, ``varm`` arrays as ``numpy.ndarray``
    * ``obsp``, ``varp`` arrays as ``scipy.sparse.csr_matrix``

    The ``X_layer_name`` is the name of the TileDB-SOMA measurement's ``X``
    collection which will be outgested to the resulting AnnData object's
    ``adata.X``.  If ``X_layer_name`` is unspecified, and the Measurement
    contains an X layer named "data", it will be returned.  If ``X_layer_name``
    is ``None``, then the return value's ``adata.X`` will be None, and
    ``adata.layers`` will be unpopulated.  If ``X_layer_name`` is a string, then
    ``adata.X`` will be taken from this layer name within the input measurement,
    and it will be an error if the measurement's ``X`` does not contain that
    layer name.

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
    to not outgest any ``uns`` keys.

    Lifecycle:
        Maturing.
    """

    s = _util.get_start_stamp()
    logging.log_io(None, "START  Experiment.to_anndata")

    if measurement_name not in experiment.ms.keys():
        raise ValueError(
            f"requested measurement name {measurement_name} not found in input: {experiment.ms.keys()}"
        )
    measurement = experiment.ms[measurement_name]
    tp = experiment.context.threadpool

    # How to choose index name for AnnData obs and var dataframes:
    # * If the desired names are passed in, use them.
    # * Else if the names used at ingest time are available, use them.
    # * Else use the default/fallback name.

    obs_df, var_df = tp.map(
        _read_dataframe,
        (experiment.obs, measurement.var),
        (obs_id_name, var_id_name),
        ("obs_id", "var_id"),
    )

    nobs = len(obs_df.index)
    nvar = len(var_df.index)

    anndata_layers_futures = {}

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

    # Problem: when the specified layer name does not exist.
    # * If they didn't specify X_layer_name at all:
    #   o Default to "data" if the measurement has that X layer name
    # * Else if they did specify X_layer_name as None:
    #   o They don't want one to be outgested
    # * Else, they specified an X_layer_name as a string:
    #   o It's an error if it doesn't exist

    # * This might be an error: if they explicitly requested X_layer_name "foo"
    #   and the X doesn't have that layer.
    # * But it might not be an error: if they didn't explicitly specify
    #   X_layer_name, and it defaulted to "data", and the experiment doesn't have that.
    # How to detect the latter?
    # * We could use **kwargs -- but that would bork the online help docs.
    # * Our consolation: check if the layer name is the _default_,
    #   and the experiment doesn't have it.
    anndata_X_future: Future[Matrix] | None = None

    if X_layer_name == MISSING:
        if "data" in measurement.X:
            anndata_X_future = _extract_X_key(measurement, "data", nobs, nvar)
    elif X_layer_name is not None:
        if X_layer_name not in measurement.X:
            raise ValueError(
                f"X_layer_name '{X_layer_name}' not found in measurement: {measurement.X.keys()}"
            )
        anndata_X_future = _extract_X_key(
            measurement, cast(str, X_layer_name), nobs, nvar
        )

    if extra_X_layer_names is not None:
        for extra_X_layer_name in extra_X_layer_names:
            if extra_X_layer_name == X_layer_name:
                continue
            assert extra_X_layer_name is not None  # appease linter; already checked
            data = _extract_X_key(measurement, extra_X_layer_name, nobs, nvar)
            anndata_layers_futures[extra_X_layer_name] = data

    if obsm_varm_width_hints is None:
        obsm_varm_width_hints = {}

    obsm = {}
    if "obsm" in measurement:
        obsm_width_hints = obsm_varm_width_hints.get("obsm", {})
        for key in measurement.obsm.keys():
            obsm[key] = tp.submit(
                _extract_obsm_or_varm,
                measurement.obsm[key],
                "obsm",
                key,
                nobs,
                obsm_width_hints,
            )

    varm = {}
    if "varm" in measurement:
        varm_width_hints = obsm_varm_width_hints.get("obsm", {})
        for key in measurement.varm.keys():
            varm[key] = tp.submit(
                _extract_obsm_or_varm,
                measurement.varm[key],
                "varm",
                key,
                nvar,
                varm_width_hints,
            )

    obsp = {}
    if "obsp" in measurement:

        def load_obsp(measurement: Measurement, key: str, nobs: int) -> sp.csr_matrix:
            A = measurement.obsp[key]
            return conversions.csr_from_coo_table(
                _read_partitioned_sparse(A, nobs),
                nobs,
                nobs,
                A.context,
            )

        for key in measurement.obsp.keys():
            obsp[key] = tp.submit(load_obsp, measurement, key, nobs)

    varp = {}
    if "varp" in measurement:

        def load_varp(measurement: Measurement, key: str, nvar: int) -> sp.csr_matrix:
            A = measurement.varp[key]
            return conversions.csr_from_coo_table(
                _read_partitioned_sparse(A, nvar),
                nvar,
                nvar,
                A.context,
            )

        for key in measurement.varp.keys():
            varp[key] = tp.submit(load_varp, measurement, key, nvar)

    uns_future: Future[dict[str, FutureUnsDictNode]] | None = None
    if "uns" in measurement:
        s = _util.get_start_stamp()
        uns_coll = cast(Collection[Any], measurement["uns"])
        logging.log_io(None, f"Start  writing uns for {uns_coll.uri}")
        uns_future = tp.submit(_extract_uns, uns_coll, uns_keys=uns_keys)
        logging.log_io(
            None,
            _util.format_elapsed(s, f"Finish writing uns for {uns_coll.uri}"),
        )

    # Resolve all futures
    obsm = _resolve_futures(obsm)
    varm = _resolve_futures(varm)
    obsp = _resolve_futures(obsp)
    varp = _resolve_futures(varp)
    anndata_X = anndata_X_future.result() if anndata_X_future else None
    anndata_layers = _resolve_futures(anndata_layers_futures)
    uns: UnsDict = (
        _resolve_futures(uns_future.result(), deep=True)
        if uns_future is not None
        else {}
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
        dtype=anndata_X.dtype if anndata_X is not None else None,
    )

    logging.log_io(None, _util.format_elapsed(s, "FINISH Experiment.to_anndata"))

    return anndata


def _extract_obsm_or_varm(
    soma_nd_array: Union[SparseNDArray, DenseNDArray],
    collection_name: str,
    element_name: str,
    num_rows: int,
    width_configs: dict[str, int],
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

    matrix_tbl = _read_partitioned_sparse(soma_nd_array, num_rows)

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
    # Four ways to get the number of columns for obsm/varm sparse matrices:
    #
    # * Explicit user specification.
    # * Shape information, if present -- this was introduced in TileDB-SOMA 1.15
    #   and solves the problem completely. However, we don't have that available
    #   for experiments created before 1.15 that have not been upgraded to have
    #   the new-shape feature.
    # * Try arithmetic on nnz / num_rows, for densely occupied sparse matrices.
    # * Try non-empty domain.
    #
    # Beyond that, we have to throw.

    description = f'{collection_name}["{element_name}"]'

    # First, try width config
    num_cols = width_configs.get(element_name, None)

    # Second, try the shape feature introduced in TileDB-SOMA 1.15
    if num_cols is None and soma_nd_array.tiledbsoma_has_upgraded_shape:
        num_cols = soma_nd_array.shape[1]

    # Third, try arithmetic on nnz / num_rows
    if num_cols is None:
        # Example:
        # * True num_rows is 100 -- this is known
        # * True num_cols is  56 -- this is unknown and needs to be solved for
        # * COO data is 5600 x 3
        # * 5600 / 100 is 56
        # This only works if the matrix is entirely occupied
        num_rows_times_width, coo_column_count = matrix_tbl.shape

        if coo_column_count != 3:
            raise SOMAError(
                f"internal error: expect COO width of 3; got {coo_column_count} for {description}"
            )

        if num_rows_times_width % num_rows == 0:
            num_cols = num_rows_times_width // num_rows

    # Fourth, try non-empty domain
    if num_cols is None:
        ned = soma_nd_array.non_empty_domain()
        num_rows = ned[0][1] + 1
        num_cols = ned[1][1] + 1

    return conversions.csr_from_coo_table(
        matrix_tbl, num_rows, num_cols, soma_nd_array.context
    ).toarray()


def _extract_uns(
    collection: Collection[Any],
    uns_keys: Sequence[str] | None = None,
    level: int = 0,
) -> dict[str, FutureUnsDictNode]:
    """
    This is a helper function for ``to_anndata`` of ``uns`` elements.
    """
    extracted: dict[str, FutureUnsDictNode] = {}
    tp = collection.context.threadpool
    for key in collection.keys():
        if level == 0 and uns_keys is not None and key not in uns_keys:
            continue

        element = collection[key]
        if isinstance(element, Collection):
            extracted[key] = _extract_uns(element, level=level + 1)
        elif isinstance(element, DataFrame):
            hint = element.metadata.get(_UNS_OUTGEST_HINT_KEY)

            def _outgest_df(
                element: Any, hint: Any, key: Any, collection: Collection[Any]
            ) -> NPNDArray | pd.DataFrame:
                if hint == _UNS_OUTGEST_HINT_1D:
                    pdf = element.read().concat().to_pandas()
                    return _outgest_uns_1d_string_array(pdf, element.uri)
                elif hint == _UNS_OUTGEST_HINT_2D:
                    pdf = element.read().concat().to_pandas()
                    return _outgest_uns_2d_string_array(pdf, element.uri)
                else:
                    if hint is not None:
                        logging.log_io_same(
                            f"Warning: uns {collection.uri}[{key!r}] has {_UNS_OUTGEST_HINT_KEY} as unrecognized {hint}: leaving this as Pandas DataFrame"
                        )
                    return _read_dataframe(element, fallback_index_name="index")

            extracted[key] = tp.submit(_outgest_df, element, hint, key, collection)

        elif isinstance(element, SparseNDArray):
            extracted[key] = tp.submit(
                lambda e: e.read().tables().concat().to_pandas(), element
            )
        elif isinstance(element, DenseNDArray):
            extracted[key] = tp.submit(lambda e: e.read().to_numpy(), element)
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
