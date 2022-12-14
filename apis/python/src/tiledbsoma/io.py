import math
import time
from typing import Callable, Optional, Union

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb

import tiledbsoma.eta as eta
import tiledbsoma.util_ann as util_ann
from tiledbsoma import (
    Collection,
    DataFrame,
    DenseNDArray,
    Experiment,
    Measurement,
    SparseNDArray,
    logging,
    util,
    util_scipy,
)

from .constants import SOMA_JOINID
from .types import Path, PlatformConfig
from .util import uri_joinpath


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
    _from_h5ad_common(
        experiment,
        input_path,
        measurement_name,
        from_anndata,
        ctx=ctx,
        platform_config=platform_config,
    )


# ----------------------------------------------------------------
def _from_h5ad_common(
    experiment: Experiment,
    input_path: Path,
    measurement_name: str,
    handler_func: Callable[
        [Experiment, ad.AnnData, str, tiledb.Ctx, Optional[PlatformConfig]], None
    ],
    ctx: Optional[tiledb.Ctx] = None,
    platform_config: Optional[PlatformConfig] = None,
) -> None:
    """
    Common code for things we do when processing a .h5ad file for ingest/update.
    """
    if isinstance(input_path, ad.AnnData):
        raise TypeError("Input path is an AnnData object -- did you want from_anndata?")

    s = util.get_start_stamp()
    logging.log_io(
        None,
        f"START  Experiment.from_h5ad {input_path}",
    )

    logging.log_io(None, f"{experiment._indent}START  READING {input_path}")

    anndata = ad.read_h5ad(input_path)

    logging.log_io(
        None,
        util.format_elapsed(s, f"{experiment._indent}FINISH READING {input_path}"),
    )

    handler_func(experiment, anndata, measurement_name, ctx, platform_config)

    logging.log_io(
        None,
        util.format_elapsed(
            s,
            f"FINISH Experiment.from_h5ad {input_path}",
        ),
    )


def _write_dataframe(
    soma_df: DataFrame,
    df: pd.DataFrame,
    id_column_name: Optional[str],
    platform_config: Optional[PlatformConfig] = None,
) -> None:
    s = util.get_start_stamp()
    logging.log_io(None, f"{soma_df._indent}START  WRITING {soma_df.uri}")

    assert not soma_df.exists()

    df[SOMA_JOINID] = np.asarray(range(len(df)), dtype=np.int64)

    df.reset_index(inplace=True)
    if id_column_name is not None:
        df.rename(columns={"index": id_column_name}, inplace=True)
    df.set_index(SOMA_JOINID, inplace=True)  # XXX MAYBE NOT?

    # TODO: This is a proposed replacement for use of tiledb.from_pandas,
    # behind a feature flag.
    #
    if soma_df._tiledb_platform_config.from_anndata_write_pandas_using_arrow:
        # categoricals are not yet well supported, so we must flatten
        for k in df:
            if df[k].dtype == "category":
                df[k] = df[k].astype(df[k].cat.categories.dtype)
        arrow_table = pa.Table.from_pandas(df)
        soma_df.create(arrow_table.schema, platform_config=platform_config)
        soma_df.write(arrow_table)

    else:
        # Legacy cut & paste - to be removed if the above code works

        offsets_filters = tiledb.FilterList(
            [tiledb.PositiveDeltaFilter(), tiledb.ZstdFilter(level=-1)]
        )
        dim_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])
        attr_filters = tiledb.FilterList([tiledb.ZstdFilter(level=-1)])

        # Force ASCII storage if string, in order to make obs/var columns queryable.
        # TODO: when UTF-8 attributes are fully supported we can remove this.
        column_types = {}
        for column_name in df.keys():
            dfc = df[column_name]
            if len(dfc) > 0 and isinstance(dfc[0], str):
                column_types[column_name] = "ascii"
            if len(dfc) > 0 and isinstance(dfc[0], bytes):
                column_types[column_name] = "bytes"

        tiledb.from_pandas(
            uri=soma_df.uri,
            dataframe=df,
            sparse=True,
            allows_duplicates=False,
            offsets_filters=offsets_filters,
            attr_filters=attr_filters,
            dim_filters=dim_filters,
            capacity=100000,
            column_types=column_types,
            ctx=soma_df._ctx,
            mode="ingest",
        )

        soma_df._common_create()  # object-type metadata etc

    logging.log_io(
        f"Wrote {soma_df.uri}",
        util.format_elapsed(s, f"{soma_df._indent}FINISH WRITING {soma_df.uri}"),
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
        soma_ndarray.write_numpy((slice(None),), nd_array)
        return

    # OR, write in chunks
    s = util.get_start_stamp()
    logging.log_io(None, f"{soma_ndarray._indent}START  ingest")

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
            "%sSTART  chunk rows %d..%d of %d (%.3f%%)"
            % (soma_ndarray._indent, i, i2 - 1, nrow, chunk_percent),
        )

        chunk = src_matrix[i:i2, :]
        if isinstance(chunk, np.ndarray):
            tensor = pa.Tensor.from_numpy(chunk)
        else:
            tensor = pa.Tensor.from_numpy(chunk.toarray())
        soma_ndarray.write_tensor((slice(i, i2), slice(None)), tensor)

        t2 = time.time()
        chunk_seconds = t2 - t1
        eta_seconds = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

        if chunk_percent < 100:
            logging.log_io(
                "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                % (soma_ndarray._indent, chunk_seconds, chunk_percent, eta_seconds),
            )

        i = i2

    logging.log_io(None, util.format_elapsed(s, f"{soma_ndarray._indent}FINISH ingest"))
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
        soma_ndarray.write_table(_coo_to_table(sp.coo_matrix(src_matrix)))
        return

    # Or, write in chunks, striding across the most efficient slice axis
    stride_axis = 0 if not sp.isspmatrix_csc(src_matrix) else 1
    dim_max_size = src_matrix.shape[stride_axis]

    s = util.get_start_stamp()
    logging.log_io(None, f"{soma_ndarray._indent}START  ingest")

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
            "%sSTART  chunk rows %d..%d of %d (%.3f%%), nnz=%d"
            % (
                soma_ndarray._indent,
                i,
                i2 - 1,
                dim_max_size,
                chunk_percent,
                chunk_coo.nnz,
            ),
        )

        soma_ndarray.write_table(_coo_to_table(chunk_coo, stride_axis, i))

        t2 = time.time()
        chunk_seconds = t2 - t1
        eta_seconds = eta_tracker.ingest_and_predict(chunk_percent, chunk_seconds)

        if chunk_percent < 100:
            logging.log_io(
                "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
                "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                % (soma_ndarray._indent, chunk_seconds, chunk_percent, eta_seconds),
            )

        i = i2

    logging.log_io(None, util.format_elapsed(s, f"{soma_ndarray._indent}FINISH ingest"))


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
    logging.log_io(None, f"{soma_ndarray._indent}START  WRITING {soma_ndarray.uri}")

    soma_ndarray.create(
        type=pa.from_numpy_dtype(src_matrix.dtype),
        shape=src_matrix.shape,
        platform_config=platform_config,
    )

    if soma_ndarray.soma_type == "SOMADenseNDArray":
        _write_matrix_to_denseNDArray(soma_ndarray, src_matrix)
    else:  # SOMmASparseNDArray
        _write_matrix_to_sparseNDArray(soma_ndarray, src_matrix)

    logging.log_io(
        f"Wrote {soma_ndarray.uri}",
        util.format_elapsed(
            s, f"{soma_ndarray._indent}FINISH WRITING {soma_ndarray.uri}"
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
    logging.log_io(None, f"{experiment._indent}START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()
    anndata = util_ann._decategoricalize(anndata)

    logging.log_io(
        None,
        util.format_elapsed(s, f"{experiment._indent}FINISH DECATEGORICALIZING"),
    )

    s = util.get_start_stamp()
    logging.log_io(None, f"{experiment._indent}START  WRITING {experiment.uri}")

    # Must be done first, to create the parent directory.
    if not experiment.exists():
        experiment.create()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # OBS
    obs = DataFrame(uri=uri_joinpath(experiment.uri, "obs"))
    _write_dataframe(
        obs, anndata.obs, id_column_name="obs_id", platform_config=platform_config
    )
    experiment.set("obs", obs)

    experiment.set(
        "ms",
        Collection(uri=uri_joinpath(experiment.uri, "ms")).create(),
    )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS
    measurement = Measurement(uri=f"{experiment.ms.uri}/{measurement_name}", ctx=ctx)
    measurement.create()
    experiment.ms.set(measurement_name, measurement)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/VAR
    var = DataFrame(uri=uri_joinpath(measurement.uri, "var"))
    _write_dataframe(
        var, anndata.var, id_column_name="var_id", platform_config=platform_config
    )
    measurement["var"] = var

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/X/DATA
    measurement["X"] = Collection(uri=uri_joinpath(measurement.uri, "X")).create()

    # TODO: more types to check?
    if isinstance(anndata.X, np.ndarray):
        ddata = DenseNDArray(uri=uri_joinpath(measurement.X.uri, "data"), ctx=ctx)
        # Code here and in else-block duplicated for linter appeasement
        create_from_matrix(ddata, anndata.X, platform_config=platform_config)
        measurement.X.set("data", ddata)
    else:
        sdata = SparseNDArray(uri=uri_joinpath(measurement.X.uri, "data"), ctx=ctx)
        create_from_matrix(sdata, anndata.X, platform_config=platform_config)
        measurement.X.set("data", sdata)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if len(anndata.obsm.keys()) > 0:  # do not create an empty collection
        measurement["obsm"] = Collection(
            uri=uri_joinpath(measurement.uri, "obsm")
        ).create()
        for key in anndata.obsm.keys():
            arr = DenseNDArray(uri=uri_joinpath(measurement.obsm.uri, key), ctx=ctx)
            create_from_matrix(arr, anndata.obsm[key], platform_config=platform_config)
            measurement.obsm.set(key, arr)

    if len(anndata.varm.keys()) > 0:  # do not create an empty collection
        measurement["varm"] = Collection(
            uri=uri_joinpath(measurement.uri, "varm")
        ).create()
        for key in anndata.varm.keys():
            darr = DenseNDArray(uri=uri_joinpath(measurement.varm.uri, key), ctx=ctx)
            create_from_matrix(darr, anndata.varm[key], platform_config=platform_config)
            measurement.varm.set(key, darr)

    if len(anndata.obsp.keys()) > 0:  # do not create an empty collection
        measurement["obsp"] = Collection(
            uri=uri_joinpath(measurement.uri, "obsp")
        ).create()
        for key in anndata.obsp.keys():
            sarr = SparseNDArray(uri=uri_joinpath(measurement.obsp.uri, key), ctx=ctx)
            create_from_matrix(sarr, anndata.obsp[key], platform_config=platform_config)
            measurement.obsp.set(key, sarr)

    if len(anndata.varp.keys()) > 0:  # do not create an empty collection
        measurement["varp"] = Collection(
            uri=uri_joinpath(measurement.uri, "varp")
        ).create()
        for key in anndata.varp.keys():
            sarr = SparseNDArray(uri=uri_joinpath(measurement.varp.uri, key), ctx=ctx)
            create_from_matrix(sarr, anndata.varp[key], platform_config=platform_config)
            measurement.varp.set(key, sarr)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RAW
    if anndata.raw is not None:
        raw_measurement = Measurement(
            uri=uri_joinpath(experiment.ms.uri, "raw"), ctx=ctx
        )
        raw_measurement.create()
        experiment.ms.set("raw", raw_measurement)

        var = DataFrame(uri=uri_joinpath(raw_measurement.uri, "var"))
        _write_dataframe(
            var,
            anndata.raw.var,
            id_column_name="var_id",
            platform_config=platform_config,
        )
        raw_measurement.set("var", var)

        raw_measurement["X"] = Collection(
            uri=uri_joinpath(raw_measurement.uri, "X")
        ).create()

        rawXdata = SparseNDArray(
            uri=uri_joinpath(raw_measurement.X.uri, "data"), ctx=ctx
        )
        create_from_matrix(rawXdata, anndata.raw.X, platform_config=platform_config)
        raw_measurement.X.set("data", rawXdata)

    # TODO: port uns from v0 to v1
    #    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #    if anndata.uns is not None:
    #        experiment.uns.from_anndata_uns(anndata.uns)
    #        experiment.set(experiment.uns)

    logging.log_io(
        f"Wrote {experiment.uri}",
        util.format_elapsed(s, f"{experiment._indent}FINISH WRITING {experiment.uri}"),
    )


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
    logging.log_io(None, f"{experiment._indent}START  write {h5ad_path}")

    anndata.write_h5ad(h5ad_path)

    logging.log_io(
        None,
        util.format_elapsed(s2, f"{experiment._indent}FINISH write {h5ad_path}"),
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

    obs_df = experiment.obs.read_as_pandas_all()
    obs_df.reset_index(inplace=True)
    obs_df.set_index("obs_id", inplace=True)
    var_df = measurement.var.read_as_pandas_all()
    var_df.reset_index(inplace=True)
    var_df.set_index("var_id", inplace=True)

    nobs = len(obs_df.index)
    nvar = len(var_df.index)

    X_data = measurement.X["data"]
    assert X_data is not None
    X_dtype = None  # some datasets have no X
    if type(X_data) == DenseNDArray:
        X_ndarray = X_data.read_numpy((slice(None), slice(None)))
        X_dtype = X_ndarray.dtype
    elif type(X_data) == SparseNDArray:
        X_mat = X_data.read_as_pandas_all()  # TODO: CSR/CSC options ...
        X_csr = util_scipy.csr_from_tiledb_df(X_mat, nobs, nvar)
        X_dtype = X_csr.dtype

    # XXX FIX OBSM/VARM SHAPES

    obsm = {}
    if measurement.obsm.exists():
        for key in measurement.obsm.keys():
            shape = measurement.obsm[key].shape
            assert len(shape) == 2
            mat = measurement.obsm[key].read_numpy((slice(None),) * len(shape))
            obsm[key] = sp.csr_array(mat)

    varm = {}
    if measurement.varm.exists():
        for key in measurement.varm.keys():
            shape = measurement.varm[key].shape
            assert len(shape) == 2
            mat = measurement.varm[key].read_numpy((slice(None),) * len(shape))
            varm[key] = sp.csr_array(mat)

    obsp = {}
    if measurement.obsp.exists():
        for key in measurement.obsp.keys():
            mat = measurement.obsp[key].read_as_pandas_all()
            obsp[key] = util_scipy.csr_from_tiledb_df(mat, nobs, nobs)

    varp = {}
    if measurement.varp.exists():
        for key in measurement.varp.keys():
            mat = measurement.varp[key].read_as_pandas_all()
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

    #    #   raw = None
    #    #   if experiment.raw.exists():
    #    #       (raw_X, raw_var_df, raw_varm) = experiment.ms['raw'].to_anndata_raw(obs_df.index)
    #    #       raw = ad.Raw(
    #    #           anndata,
    #    #           X=raw_X,
    #    #           var=raw_var_df,
    #    #           varm=raw_varm,
    #    #       )
    #
    #    # TODO: PORT FROM V0 TO V1
    #    # uns = experiment.uns.to_dict_of_matrices()

    #    anndata = ad.AnnData(
    #        #       X=anndata.X,
    #        #       dtype=None if anndata.X is None else anndata.X.dtype,  # some datasets have no X
    #        #       obs=anndata.obs,
    #        #       var=anndata.var,
    #        #       obsm=anndata.obsm,
    #        #       obsp=anndata.obsp,
    #        #       varm=anndata.varm,
    #        #       varp=anndata.varp,
    #        #       raw=raw,
    #        #       uns=uns,
    #    )

    logging.log_io(
        None,
        util.format_elapsed(s, "FINISH Experiment.to_anndata"),
    )

    return anndata
