import math
import time
from typing import Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb

import tiledbsoma.eta as eta
import tiledbsoma.util as util
import tiledbsoma.util_ann as util_ann
import tiledbsoma.util_tiledb as util_tiledb
from tiledbsoma import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    Measurement,
    SparseNDArray,
    logging,
    util_scipy,
)

from .constants import SOMA_JOINID
from .types import Path, PlatformConfig


# ----------------------------------------------------------------
def from_h5ad(
    experiment: Experiment,
    input_path: Path,
    measurement_name: str,
    ctx: Optional[tiledb.Ctx] = None,
    platform_config: Optional[PlatformConfig] = None,
) -> None:
    """
    Reads an .h5ad file and writes to a TileDB group structure.
    """
    if isinstance(input_path, ad.AnnData):
        raise TypeError("Input path is an AnnData object -- did you want from_anndata?")

    s = util.get_start_stamp()
    logging.log_io(
        None,
        f"START  Experiment.from_h5ad {input_path}",
    )

    logging.log_io(None, f"START  READING {input_path}")

    anndata = ad.read_h5ad(input_path, backed="r")

    logging.log_io(
        None,
        util.format_elapsed(s, f"FINISH READING {input_path}"),
    )

    from_anndata(experiment, anndata, measurement_name, ctx, platform_config)

    logging.log_io(
        None,
        util.format_elapsed(
            s,
            f"FINISH Experiment.from_h5ad {input_path}",
        ),
    )


# ----------------------------------------------------------------
def from_anndata(
    experiment: Experiment,
    anndata: ad.AnnData,
    measurement_name: str,
    ctx: Optional[tiledb.Ctx] = None,
    platform_config: Optional[PlatformConfig] = None,
) -> None:
    """
    Top-level writer method for creating a TileDB group for a ``Experiment`` object.
    """
    if not isinstance(anndata, ad.AnnData):
        raise TypeError(
            "Second argument is not an AnnData object -- did you want from_h5ad?"
        )

    # Without _at least_ an index, there is nothing to indicate the dimension indices.
    if anndata.obs.index.empty or anndata.var.index.empty:
        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

    s = util.get_start_stamp()
    logging.log_io(None, "START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()

    logging.log_io(
        None,
        util.format_elapsed(s, "FINISH DECATEGORICALIZING"),
    )

    s = util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {experiment.uri}")

    # Must be done first, to create the parent directory.
    if not experiment.exists():
        experiment.create()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # OBS
    obs = DataFrame(uri=util.uri_joinpath(experiment.uri, "obs"))
    _write_dataframe(
        obs,
        util_ann._decategoricalize_obs_or_var(anndata.obs),
        id_column_name="obs_id",
        platform_config=platform_config,
    )
    experiment.set("obs", obs)

    experiment.set(
        "ms",
        Collection(uri=util.uri_joinpath(experiment.uri, "ms")).create(),
    )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS
    measurement = Measurement(uri=f"{experiment.ms.uri}/{measurement_name}", ctx=ctx)
    measurement.create()
    experiment.ms.set(measurement_name, measurement)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/VAR
    var = DataFrame(uri=util.uri_joinpath(measurement.uri, "var"))
    _write_dataframe(
        var,
        util_ann._decategoricalize_obs_or_var(anndata.var),
        id_column_name="var_id",
        platform_config=platform_config,
    )
    measurement["var"] = var

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/X/DATA
    measurement["X"] = Collection(uri=util.uri_joinpath(measurement.uri, "X")).create()

    # TODO: more types to check?
    if isinstance(anndata.X, np.ndarray):
        ddata = DenseNDArray(uri=util.uri_joinpath(measurement.X.uri, "data"), ctx=ctx)
        # Code here and in else-block duplicated for linter appeasement
        create_from_matrix(ddata, anndata.X[:], platform_config=platform_config)
        measurement.X.set("data", ddata)
    else:
        sdata = SparseNDArray(uri=util.uri_joinpath(measurement.X.uri, "data"), ctx=ctx)
        create_from_matrix(sdata, anndata.X[:], platform_config=platform_config)
        measurement.X.set("data", sdata)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if len(anndata.obsm.keys()) > 0:  # do not create an empty collection
        measurement["obsm"] = Collection(
            uri=util.uri_joinpath(measurement.uri, "obsm")
        ).create()
        for key in anndata.obsm.keys():
            arr = DenseNDArray(
                uri=util.uri_joinpath(measurement.obsm.uri, key), ctx=ctx
            )
            create_from_matrix(
                arr,
                util_tiledb.to_tiledb_supported_array_type(anndata.obsm[key]),
                platform_config=platform_config,
            )
            measurement.obsm.set(key, arr)

    if len(anndata.varm.keys()) > 0:  # do not create an empty collection
        measurement["varm"] = Collection(
            uri=util.uri_joinpath(measurement.uri, "varm")
        ).create()
        for key in anndata.varm.keys():
            darr = DenseNDArray(
                uri=util.uri_joinpath(measurement.varm.uri, key), ctx=ctx
            )
            create_from_matrix(
                darr,
                util_tiledb.to_tiledb_supported_array_type(anndata.varm[key]),
                platform_config=platform_config,
            )
            measurement.varm.set(key, darr)

    if len(anndata.obsp.keys()) > 0:  # do not create an empty collection
        measurement["obsp"] = Collection(
            uri=util.uri_joinpath(measurement.uri, "obsp")
        ).create()
        for key in anndata.obsp.keys():
            sarr = SparseNDArray(
                uri=util.uri_joinpath(measurement.obsp.uri, key), ctx=ctx
            )
            create_from_matrix(
                sarr,
                util_tiledb.to_tiledb_supported_array_type(anndata.obsp[key]),
                platform_config=platform_config,
            )
            measurement.obsp.set(key, sarr)

    if len(anndata.varp.keys()) > 0:  # do not create an empty collection
        measurement["varp"] = Collection(
            uri=util.uri_joinpath(measurement.uri, "varp")
        ).create()
        for key in anndata.varp.keys():
            sarr = SparseNDArray(
                uri=util.uri_joinpath(measurement.varp.uri, key), ctx=ctx
            )
            create_from_matrix(
                sarr,
                util_tiledb.to_tiledb_supported_array_type(anndata.varp[key]),
                platform_config=platform_config,
            )
            measurement.varp.set(key, sarr)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RAW
    if anndata.raw is not None:
        raw_measurement = Measurement(
            uri=util.uri_joinpath(experiment.ms.uri, "raw"), ctx=ctx
        )
        raw_measurement.create()
        experiment.ms.set("raw", raw_measurement)

        var = DataFrame(uri=util.uri_joinpath(raw_measurement.uri, "var"))
        _write_dataframe(
            var,
            util_ann._decategoricalize_obs_or_var(anndata.raw.var),
            id_column_name="var_id",
            platform_config=platform_config,
        )
        raw_measurement.set("var", var)

        raw_measurement["X"] = Collection(
            uri=util.uri_joinpath(raw_measurement.uri, "X")
        ).create()

        rawXdata = SparseNDArray(
            uri=util.uri_joinpath(raw_measurement.X.uri, "data"), ctx=ctx
        )
        create_from_matrix(rawXdata, anndata.raw.X[:], platform_config=platform_config)
        raw_measurement.X.set("data", rawXdata)

    # TODO: port uns from main-old to main
    #    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #    if anndata.uns is not None:
    #        experiment.uns.from_anndata_uns(anndata.uns)
    #        experiment.set(experiment.uns)

    logging.log_io(
        f"Wrote {experiment.uri}",
        util.format_elapsed(s, f"FINISH WRITING {experiment.uri}"),
    )


def _write_dataframe(
    soma_df: DataFrame,
    df: pd.DataFrame,
    id_column_name: Optional[str],
    platform_config: Optional[PlatformConfig] = None,
) -> None:
    s = util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {soma_df.uri}")

    assert not soma_df.exists()

    df[SOMA_JOINID] = np.asarray(range(len(df)), dtype=np.int64)

    df.reset_index(inplace=True)
    if id_column_name is not None:
        df.rename(columns={"index": id_column_name}, inplace=True)
    df.set_index(SOMA_JOINID, inplace=True)  # XXX MAYBE NOT?

    # Categoricals are not yet well supported, so we must flatten
    for k in df:
        if df[k].dtype == "category":
            df[k] = df[k].astype(df[k].cat.categories.dtype)
    arrow_table = pa.Table.from_pandas(df)
    soma_df.create(arrow_table.schema, platform_config=platform_config)
    soma_df.write(arrow_table)

    logging.log_io(
        f"Wrote {soma_df.uri}",
        util.format_elapsed(s, f"FINISH WRITING {soma_df.uri}"),
    )


def create_from_matrix(
    soma_ndarray: Union[DenseNDArray, SparseNDArray],
    src_matrix: Union[np.ndarray, sp.csr_matrix, sp.csc_matrix],
    platform_config: Optional[PlatformConfig] = None,
) -> None:
    """
    Create and populate the ``soma_matrix`` from the contents of ``src_matrix``.
    """
    assert not soma_ndarray.exists()
    assert src_matrix.ndim == 2
    assert soma_ndarray.soma_type in ("SOMADenseNDArray", "SOMASparseNDArray")

    s = util.get_start_stamp()
    logging.log_io(None, f"START  WRITING {soma_ndarray.uri}")

    soma_ndarray.create(
        type=pa.from_numpy_dtype(src_matrix.dtype),
        shape=src_matrix.shape,
        platform_config=platform_config,
    )

    if isinstance(soma_ndarray, DenseNDArray):
        _write_matrix_to_denseNDArray(soma_ndarray, src_matrix)
    else:  # SOMmASparseNDArray
        _write_matrix_to_sparseNDArray(soma_ndarray, src_matrix)

    logging.log_io(
        f"Wrote {soma_ndarray.uri}",
        util.format_elapsed(s, f"FINISH WRITING {soma_ndarray.uri}"),
    )


def _write_matrix_to_denseNDArray(
    soma_ndarray: DenseNDArray,
    src_matrix: Union[np.ndarray, sp.csr_matrix, sp.csc_matrix],
) -> None:
    """Write a matrix to an empty DenseNDArray"""

    # Write all at once?
    if soma_ndarray._tiledb_platform_config.write_X_chunked:
        if isinstance(src_matrix, np.ndarray):
            nd_array = src_matrix
        else:
            nd_array = src_matrix.toarray()
        soma_ndarray.write((slice(None),), pa.Tensor.from_numpy(nd_array))
        return

    # OR, write in chunks
    s = util.get_start_stamp()
    logging.log_io(None, "START  ingest")

    eta_tracker = eta.Tracker()
    nrow, ncol = src_matrix.shape
    i = 0
    # number of rows to chunk by. Dense writes, so this is a constant.
    chunk_size = int(
        math.ceil(soma_ndarray._tiledb_platform_config.goal_chunk_nnz / ncol)
    )
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

        chunk = src_matrix[i:i2, :]
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

    logging.log_io(None, util.format_elapsed(s, "FINISH ingest"))
    return


def _write_matrix_to_sparseNDArray(
    soma_ndarray: SparseNDArray,
    src_matrix: Union[np.ndarray, sp.csr_matrix, sp.csc_matrix],
) -> None:
    """Write a matrix to an empty DenseNDArray"""

    def _coo_to_table(mat_coo: sp.coo_matrix, axis: int = 0, base: int = 0) -> pa.Table:
        pydict = {
            "soma_data": mat_coo.data,
            "soma_dim_0": mat_coo.row + base if base > 0 and axis == 0 else mat_coo.row,
            "soma_dim_1": mat_coo.col + base if base > 0 and axis == 1 else mat_coo.col,
        }
        return pa.Table.from_pydict(pydict)

    def _find_chunk_size(
        mat: Union[np.ndarray, sp.csr_matrix, sp.csc_matrix],
        start_index: int,
        axis: int,
        goal_chunk_nnz: int,
    ) -> int:
        if isinstance(mat, np.ndarray):
            return int(math.ceil(goal_chunk_nnz / mat.shape[axis]))
        else:
            return util_scipy.find_sparse_chunk_size(
                mat, start_index, axis, goal_chunk_nnz
            )

    # Write all at once?
    if not soma_ndarray._tiledb_platform_config.write_X_chunked:
        soma_ndarray.write(_coo_to_table(sp.coo_matrix(src_matrix)))
        return

    # Or, write in chunks, striding across the most efficient slice axis
    stride_axis = 0 if not sp.isspmatrix_csc(src_matrix) else 1
    dim_max_size = src_matrix.shape[stride_axis]

    s = util.get_start_stamp()
    logging.log_io(None, "START  ingest")

    eta_tracker = eta.Tracker()
    goal_chunk_nnz = soma_ndarray._tiledb_platform_config.goal_chunk_nnz
    coords = [slice(None), slice(None)]
    i = 0
    while i < dim_max_size:
        t1 = time.time()

        # chunk size on the stride axis
        chunk_size = _find_chunk_size(src_matrix, i, stride_axis, goal_chunk_nnz)
        i2 = i + chunk_size

        coords[stride_axis] = slice(i, i2)
        chunk_coo = sp.coo_matrix(src_matrix[tuple(coords)])

        # print doubly-inclusive lo..hi like 0..17 and 18..31.
        chunk_percent = min(100, 100 * (i2 - 1) / dim_max_size)
        logging.log_io(
            None,
            "START  chunk rows %d..%d of %d (%.3f%%), nnz=%d"
            % (
                i,
                i2 - 1,
                dim_max_size,
                chunk_percent,
                chunk_coo.nnz,
            ),
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

    logging.log_io(None, util.format_elapsed(s, "FINISH ingest"))


# ----------------------------------------------------------------
def to_h5ad(
    experiment: Experiment,
    h5ad_path: Path,
    measurement_name: str,
    ctx: Optional[tiledb.Ctx] = None,
) -> None:
    """
    Converts the experiment group to anndata format and writes it to the specified .h5ad file.
    """

    s = util.get_start_stamp()
    logging.log_io(None, f"START  Experiment.to_h5ad -> {h5ad_path}")

    anndata = to_anndata(experiment, measurement_name=measurement_name)

    s2 = util.get_start_stamp()
    logging.log_io(None, f"START  write {h5ad_path}")

    anndata.write_h5ad(h5ad_path)

    logging.log_io(
        None,
        util.format_elapsed(s2, f"FINISH write {h5ad_path}"),
    )

    logging.log_io(
        None,
        util.format_elapsed(s, f"FINISH Experiment.to_h5ad -> {h5ad_path}"),
    )


# ----------------------------------------------------------------
def to_anndata(
    experiment: Experiment,
    *,
    # TODO: set a better name as capitalized-const
    # TODO: maybe if there are multiple measurements, default to the first one not named 'raw'
    measurement_name: str,
    ctx: Optional[tiledb.Ctx] = None,
) -> ad.AnnData:
    """
    Converts the experiment group to anndata. Choice of matrix formats is following what we often see in input .h5ad files:

    * X as ``scipy.sparse.csr_matrix``
    * obs,var as ``pandas.dataframe``
    * obsm,varm arrays as ``numpy.ndarray``
    * obsp,varp arrays as ``scipy.sparse.csr_matrix``
    """

    s = util.get_start_stamp()
    logging.log_io(None, "START  Experiment.to_anndata")

    measurement = experiment.ms[measurement_name]

    obs_df = experiment.obs.read().concat().to_pandas()
    obs_df.drop([SOMA_JOINID], axis=1, inplace=True)
    obs_df.set_index("obs_id", inplace=True)

    var_df = measurement.var.read().concat().to_pandas()
    var_df.drop([SOMA_JOINID], axis=1, inplace=True)
    var_df.set_index("var_id", inplace=True)

    nobs = len(obs_df.index)
    nvar = len(var_df.index)

    X_data = measurement.X["data"]
    X_csr = None
    assert X_data is not None
    X_dtype = None  # some datasets have no X
    if isinstance(X_data, DenseNDArray):
        X_ndarray = X_data.read((slice(None), slice(None))).to_numpy()
        X_dtype = X_ndarray.dtype
    elif isinstance(X_data, SparseNDArray):
        X_mat = X_data.read().tables().concat().to_pandas()  # TODO: CSR/CSC options ...
        X_csr = util_scipy.csr_from_tiledb_df(X_mat, nobs, nvar)
        X_dtype = X_csr.dtype
    else:
        raise TypeError(f"Unexpected NDArray type {type(X_data)}")

    obsm = {}
    if "obsm" in measurement and measurement.obsm.exists():
        for key in measurement.obsm.keys():
            shape = measurement.obsm[key].shape
            assert len(shape) == 2
            mat = measurement.obsm[key].read((slice(None),) * len(shape)).to_numpy()
            # The spelling `sp.csr_array` is more idiomatic but doesn't exist until Python 3.8
            obsm[key] = sp.csr_matrix(mat)

    varm = {}
    if "varm" in measurement and measurement.varm.exists():
        for key in measurement.varm.keys():
            shape = measurement.varm[key].shape
            assert len(shape) == 2
            mat = measurement.varm[key].read((slice(None),) * len(shape)).to_numpy()
            # The spelling `sp.csr_array` is more idiomatic but doesn't exist until Python 3.8
            varm[key] = sp.csr_matrix(mat)

    obsp = {}
    if "obsp" in measurement and measurement.obsp.exists():
        for key in measurement.obsp.keys():
            mat = measurement.obsp[key].read().tables().concat().to_pandas()
            obsp[key] = util_scipy.csr_from_tiledb_df(mat, nobs, nobs)

    varp = {}
    if "varp" in measurement and measurement.varp.exists():
        for key in measurement.varp.keys():
            mat = measurement.varp[key].read().tables().concat().to_pandas()
            varp[key] = util_scipy.csr_from_tiledb_df(mat, nvar, nvar)

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

    logging.log_io(
        None,
        util.format_elapsed(s, "FINISH Experiment.to_anndata"),
    )

    return anndata
