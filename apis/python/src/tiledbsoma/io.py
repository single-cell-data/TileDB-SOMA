from typing import Callable, Optional

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb

import tiledbsoma.util_ann as util_ann
from tiledbsoma import (
    SOMACollection,
    SOMADataFrame,
    SOMADenseNdArray,
    SOMAExperiment,
    SOMAMeasurement,
    SOMASparseNdArray,
    logging,
    util,
    util_scipy,
)

from .constants import SOMA_JOINID, SOMA_ROWID
from .types import Path
from .util import uri_joinpath


# ----------------------------------------------------------------
def from_h5ad(
    experiment: SOMAExperiment,
    input_path: Path,
    measurement_name: str,
    ctx: Optional[tiledb.Ctx] = None,
) -> None:
    """
    Reads an .h5ad file and writes to a TileDB group structure.
    """
    _from_h5ad_common(experiment, input_path, measurement_name, from_anndata, ctx=ctx)


# ----------------------------------------------------------------
def _from_h5ad_common(
    experiment: SOMAExperiment,
    input_path: Path,
    measurement_name: str,
    handler_func: Callable[[SOMAExperiment, ad.AnnData, str, tiledb.Ctx], None],
    ctx: Optional[tiledb.Ctx] = None,
) -> None:
    """
    Common code for things we do when processing a .h5ad file for ingest/update.
    """
    if isinstance(input_path, ad.AnnData):
        raise Exception("Input path is an AnnData object -- did you want from_anndata?")

    s = util.get_start_stamp()
    logging.log_io(
        None,
        f"START  SOMAExperiment.from_h5ad {input_path}",
    )

    logging.log_io(None, f"{experiment._indent}START  READING {input_path}")

    anndata = ad.read_h5ad(input_path)

    logging.log_io(
        None,
        util.format_elapsed(s, f"{experiment._indent}FINISH READING {input_path}"),
    )

    handler_func(experiment, anndata, measurement_name, ctx)

    logging.log_io(
        None,
        util.format_elapsed(
            s,
            f"FINISH SOMAExperiment.from_h5ad {input_path}",
        ),
    )


def _write_dataframe(
    soma_df: SOMADataFrame, df: pd.DataFrame, id_column_name: Optional[str]
) -> None:
    s = util.get_start_stamp()
    logging.log_io(None, f"{soma_df._indent}START  WRITING {soma_df.uri}")

    assert not soma_df.exists()

    df[SOMA_ROWID] = np.asarray(range(len(df)), dtype=np.int64)
    df[SOMA_JOINID] = np.asarray(range(len(df)), dtype=np.int64)

    df.reset_index(inplace=True)
    if id_column_name is not None:
        df.rename(columns={"index": id_column_name}, inplace=True)
    df.set_index(SOMA_ROWID, inplace=True)

    # TODO: This is a proposed replacement for use of tiledb.from_pandas,
    # behind a feature flag.
    #
    if soma_df._tiledb_platform_config.from_anndata_write_pandas_using_arrow:
        # categoricals are not yet well supported, so we must flatten
        for k in df:
            if df[k].dtype == "category":
                df[k] = df[k].astype(df[k].cat.categories.dtype)
        arrow_table = pa.Table.from_pandas(df)
        soma_df.create(arrow_table.schema)
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
            allows_duplicates=soma_df._tiledb_platform_config.allows_duplicates,
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


# ----------------------------------------------------------------
def from_anndata(
    experiment: SOMAExperiment,
    anndata: ad.AnnData,
    measurement_name: str,
    ctx: Optional[tiledb.Ctx] = None,
) -> None:
    """
    Top-level writer method for creating a TileDB group for a ``SOMAExperiment`` object.
    """
    if not isinstance(anndata, ad.AnnData):
        raise Exception(
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
    obs = SOMADataFrame(uri=uri_joinpath(experiment.uri, "obs"))
    _write_dataframe(obs, anndata.obs, id_column_name="obs_id")
    experiment.set("obs", obs)

    experiment.set(
        "ms", SOMACollection(uri=uri_joinpath(experiment.uri, "ms")).create()
    )

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS
    measurement = SOMAMeasurement(
        uri=f"{experiment.ms.uri}/{measurement_name}", ctx=ctx
    )
    measurement.create()
    experiment.ms.set(measurement_name, measurement)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/VAR
    var = SOMADataFrame(uri=uri_joinpath(measurement.uri, "var"))
    _write_dataframe(var, anndata.var, id_column_name="var_id")
    measurement["var"] = var

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/X/DATA
    measurement["X"] = SOMACollection(uri=uri_joinpath(measurement.uri, "X")).create()

    # TODO: more types to check?
    if isinstance(anndata.X, np.ndarray):
        ddata = SOMADenseNdArray(uri=uri_joinpath(measurement.X.uri, "data"), ctx=ctx)
        # Code here and in else-block duplicated for linter appeasement
        ddata.from_matrix(anndata.X)
        measurement.X.set("data", ddata)
    else:
        sdata = SOMASparseNdArray(uri=uri_joinpath(measurement.X.uri, "data"), ctx=ctx)
        sdata.from_matrix(anndata.X)
        measurement.X.set("data", sdata)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # TODO: port from v0

    if len(anndata.obsm.keys()) > 0:  # do not create an empty collection
        measurement["obsm"] = SOMACollection(
            uri=uri_joinpath(measurement.uri, "obsm")
        ).create()
        for key in anndata.obsm.keys():
            arr = SOMADenseNdArray(uri=uri_joinpath(measurement.obsm.uri, key), ctx=ctx)
            arr.from_matrix(anndata.obsm[key])
            measurement.obsm.set(key, arr)

    if len(anndata.varm.keys()) > 0:  # do not create an empty collection
        measurement["varm"] = SOMACollection(
            uri=uri_joinpath(measurement.uri, "varm")
        ).create()
        for key in anndata.varm.keys():
            darr = SOMADenseNdArray(
                uri=uri_joinpath(measurement.varm.uri, key), ctx=ctx
            )
            darr.from_matrix(anndata.varm[key])
            measurement.varm.set(key, darr)

    if len(anndata.obsp.keys()) > 0:  # do not create an empty collection
        measurement["obsp"] = SOMACollection(
            uri=uri_joinpath(measurement.uri, "obsp")
        ).create()
        for key in anndata.obsp.keys():
            sarr = SOMASparseNdArray(
                uri=uri_joinpath(measurement.obsp.uri, key), ctx=ctx
            )
            sarr.from_matrix(anndata.obsp[key])
            measurement.obsp.set(key, sarr)

    if len(anndata.varp.keys()) > 0:  # do not create an empty collection
        measurement["varp"] = SOMACollection(
            uri=uri_joinpath(measurement.uri, "varp")
        ).create()
        for key in anndata.varp.keys():
            sarr = SOMASparseNdArray(
                uri=uri_joinpath(measurement.varp.uri, key), ctx=ctx
            )
            sarr.from_matrix(anndata.varp[key])
            measurement.varp.set(key, sarr)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RAW
    if anndata.raw is not None:
        raw_measurement = SOMAMeasurement(
            uri=uri_joinpath(experiment.ms.uri, "raw"), ctx=ctx
        )
        raw_measurement.create()
        experiment.ms.set("raw", raw_measurement)

        var = SOMADataFrame(uri=uri_joinpath(raw_measurement.uri, "var"))
        _write_dataframe(var, anndata.raw.var, id_column_name="var_id")
        raw_measurement.set("var", var)

        raw_measurement["X"] = SOMACollection(
            uri=uri_joinpath(raw_measurement.uri, "X")
        ).create()

        rawXdata = SOMASparseNdArray(
            uri=uri_joinpath(raw_measurement.X.uri, "data"), ctx=ctx
        )
        rawXdata.from_matrix(anndata.raw.X)
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
    experiment: SOMAExperiment,
    h5ad_path: Path,
    measurement_name: str,
    ctx: Optional[tiledb.Ctx] = None,
) -> None:
    """
    Converts the experiment group to anndata format and writes it to the specified .h5ad file.
    """

    s = util.get_start_stamp()
    logging.log_io(None, f"START  SOMAExperiment.to_h5ad -> {h5ad_path}")

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
        util.format_elapsed(s, f"FINISH SOMAExperiment.to_h5ad -> {h5ad_path}"),
    )


# ----------------------------------------------------------------
def to_anndata(
    experiment: SOMAExperiment,
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
    logging.log_io(None, "START  SOMAExperiment.to_anndata")

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
    if type(X_data) == SOMADenseNdArray:
        X_ndarray = X_data.read_numpy((slice(None), slice(None)))
        X_dtype = X_ndarray.dtype
    elif type(X_data) == SOMASparseNdArray:
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
        util.format_elapsed(s, "FINISH SOMAExperiment.to_anndata"),
    )

    return anndata


# ----------------------------------------------------------------
# def to_anndata_from_raw(experiment: SOMAExperiment) -> ad.AnnData:
#    """
#    Extract only the raw parts as a new AnnData object.
#    """
#
#    obs_df = experiment.obs.df()
#    var_df = experiment.raw.var.df()
#    data = experiment.raw.X["data"]
#    assert data is not None
#    X_mat = data.to_csr_matrix(obs_df.index, var_df.index)
#
#    return ad.AnnData(X=X_mat, obs=obs_df, var=var_df)
