from typing import Callable, Optional

import anndata as ad
import numpy as np
import tiledb

import tiledbsoma.util_ann as util_ann
from tiledbsoma import (
    SOMADenseNdArray,
    SOMAExperiment,
    SOMAMeasurement,
    SOMASparseNdArray,
    logging,
    util,
    util_scipy,
)

from .types import Path

# import scanpy
# import tiledb

# import tiledbsoma
# import tiledbsoma.logging
# import tiledbsoma.util


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
        f"START  SOMAExperiment.from_h5ad {input_path} -> {experiment._nested_name}",
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
            f"FINISH SOMAExperiment.from_h5ad {input_path} -> {experiment._nested_name}",
        ),
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
    logging.log_io(
        None, f"{experiment._indent}START  WRITING {experiment._nested_name}"
    )

    # Must be done first, to create the parent directory
    if not experiment.exists():
        experiment.create()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # OBS
    experiment.obs.write_all_from_pandas(
        dataframe=anndata.obs, extent=256, id_column_name="obs_id"
    )
    experiment.set(experiment.obs)

    experiment.ms.create()
    experiment.set(experiment.ms)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS
    measurement = SOMAMeasurement(
        uri=f"{experiment.ms.get_uri()}/{measurement_name}", ctx=ctx
    )
    measurement.create()
    experiment.ms.set(measurement)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/VAR
    measurement.var.write_all_from_pandas(
        dataframe=anndata.var, extent=2048, id_column_name="var_id"
    )
    measurement.set(measurement.var)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/X/DATA
    measurement.X.create()
    measurement.set(measurement.X)

    # TODO: more types to check?
    if isinstance(anndata.X, np.ndarray):
        ddata = SOMADenseNdArray(uri=f"{measurement.X.get_uri()}/data", ctx=ctx)
        # Code here and in else-block duplicated for linter appeasement
        ddata.from_matrix(anndata.X)
        measurement.X.set(ddata)
    else:
        sdata = SOMASparseNdArray(uri=f"{measurement.X.get_uri()}/data", ctx=ctx)
        sdata.from_matrix(anndata.X)
        measurement.X.set(sdata)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # TODO: port from v0

    if len(anndata.obsm.keys()) > 0:  # do not create an empty collection
        measurement.obsm.create()
        for key in anndata.obsm.keys():
            arr = SOMADenseNdArray(uri=f"{measurement.obsm.get_uri()}/{key}", ctx=ctx)
            arr.from_matrix(anndata.obsm[key])
            measurement.obsm.set(arr)
        measurement.set(measurement.obsm)

    if len(anndata.varm.keys()) > 0:  # do not create an empty collection
        measurement.varm.create()
        for key in anndata.varm.keys():
            darr = SOMADenseNdArray(uri=f"{measurement.varm.get_uri()}/{key}", ctx=ctx)
            darr.from_matrix(anndata.varm[key])
            measurement.varm.set(darr)
        measurement.set(measurement.varm)

    if len(anndata.obsp.keys()) > 0:  # do not create an empty collection
        measurement.obsp.create()
        for key in anndata.obsp.keys():
            sarr = SOMASparseNdArray(uri=f"{measurement.obsp.get_uri()}/{key}", ctx=ctx)
            sarr.from_matrix(anndata.obsp[key])
            measurement.obsp.set(sarr)
        measurement.set(measurement.obsp)

    if len(anndata.varp.keys()) > 0:  # do not create an empty collection
        measurement.varp.create()
        for key in anndata.varp.keys():
            sarr = SOMASparseNdArray(uri=f"{measurement.varp.get_uri()}/{key}", ctx=ctx)
            sarr.from_matrix(anndata.varp[key])
            measurement.varp.set(sarr)
        measurement.set(measurement.varp)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RAW
    if anndata.raw is not None:
        raw_measurement = SOMAMeasurement(uri=f"{experiment.ms.get_uri()}/raw", ctx=ctx)
        raw_measurement.create()
        experiment.ms.set(raw_measurement)

        raw_measurement.var.write_all_from_pandas(
            dataframe=anndata.raw.var, extent=2048, id_column_name="var_id"
        )
        raw_measurement.set(raw_measurement.var)

        raw_measurement.X.create()
        raw_measurement.set(raw_measurement.X)

        rawXdata = SOMASparseNdArray(uri=f"{raw_measurement.X.get_uri()}/data", ctx=ctx)
        rawXdata.from_matrix(anndata.raw.X)
        raw_measurement.X.set(rawXdata)

    # TODO: port uns from v0 to v1
    #    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #    if anndata.uns is not None:
    #        experiment.uns.from_anndata_uns(anndata.uns)
    #        experiment.set(experiment.uns)

    logging.log_io(
        f"Wrote {experiment._nested_name}",
        util.format_elapsed(
            s, f"{experiment._indent}FINISH WRITING {experiment._nested_name}"
        ),
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
    logging.log_io(
        None, f"START  SOMAExperiment.to_h5ad {experiment._nested_name} -> {h5ad_path}"
    )

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
        util.format_elapsed(
            s, f"FINISH SOMAExperiment.to_h5ad {experiment._nested_name} -> {h5ad_path}"
        ),
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
    logging.log_io(None, f"START  SOMAExperiment.to_anndata {experiment._nested_name}")

    measurement = experiment.ms[measurement_name]

    obs_df = experiment.obs.read_as_pandas_all(id_column_name="obs_id")
    var_df = measurement.var.read_as_pandas_all(id_column_name="var_id")

    nobs = len(obs_df.index)
    nvar = len(var_df.index)

    X_data = measurement.X["data"]
    assert X_data is not None
    X_mat = X_data.read_as_pandas_all()  # TODO: CSR/CSC options ...
    X_csr = util_scipy.csr_from_tiledb_df(X_mat, nobs, nvar)

    # XXX FIX OBSM/VARM SHAPES

    obsm = {}
    if measurement.obsm.exists():
        for key in measurement.obsm.keys():
            shape = measurement.obsm[key].shape
            assert len(shape) == 2
            ncols = shape[1]
            mat = measurement.obsm[key].read_as_pandas_all()
            print("OBSM", key, mat.shape)
            obsm[key] = util_scipy.csr_from_tiledb_df(mat, nobs, ncols)

    varm = {}
    if measurement.varm.exists():
        for key in measurement.varm.keys():
            shape = measurement.varm[key].shape
            assert len(shape) == 2
            ncols = shape[1]
            mat = measurement.varm[key].read_as_pandas_all()
            print("VARM", key, mat.shape)
            varm[key] = util_scipy.csr_from_tiledb_df(mat, nvar, ncols)

    obsp = {}
    if measurement.obsp.exists():
        for key in measurement.obsp.keys():
            mat = measurement.obsp[key].read_as_pandas_all()
            print("OBSP", key, mat.shape)
            obsp[key] = util_scipy.csr_from_tiledb_df(mat, nobs, nobs)

    varp = {}
    if measurement.varp.exists():
        for key in measurement.varp.keys():
            mat = measurement.varp[key].read_as_pandas_all()
            print("VARP", key, mat.shape)
            varp[key] = util_scipy.csr_from_tiledb_df(mat, nvar, nvar)

    anndata = ad.AnnData(
        X=X_csr,
        obs=obs_df,
        var=var_df,
        obsm=obsm,
        varm=varm,
        obsp=obsp,
        varp=varp,
        dtype=None if X_csr is None else X_csr.dtype,  # some datasets have no X
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
        util.format_elapsed(
            s, f"FINISH SOMAExperiment.to_anndata {experiment._nested_name}"
        ),
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
