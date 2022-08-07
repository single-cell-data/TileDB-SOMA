from typing import Callable, Optional

import anndata as ad

import tiledbsc.v1.util_ann as util_ann
from tiledbsc.v1 import (
    SOMADenseNdArray,
    SOMAExperiment,
    SOMAMeasurement,
    SOMASparseNdArray,
    logging,
    util,
)

from .types import Path

# import scanpy
# import tiledb

# import tiledbsc
# import tiledbsc.logging
# import tiledbsc.util


# ----------------------------------------------------------------
def from_h5ad(experiment: SOMAExperiment, input_path: Path) -> None:
    """
    Reads an .h5ad file and writes to a TileDB group structure.
    """
    _from_h5ad_common(experiment, input_path, from_anndata)


# ----------------------------------------------------------------
# def from_h5ad_update_obs_and_var(experiment: SOMAExperiment, input_path: Path) -> None:
#    """
#    Rewrites obs and var from the specified .h5ad file, leaving all other data in place. Useful for
#    updating schema/compression/etc. within an existing dataset.
#    """
#    _from_h5ad_common(experiment, input_path, from_anndata_update_obs_and_var)
#

# ----------------------------------------------------------------
# TODO TEMP WIP
def _from_h5ad_common(
    experiment: SOMAExperiment,
    input_path: Path,
    handler_func: Callable[[SOMAExperiment, ad.AnnData], None],
) -> None:
    """
    Common code for things we do when processing a .h5ad file for ingest/update.
    """
    s = util.get_start_stamp()
    logging.log_io(
        None,
        f"START  v1.SOMAExperiment.from_h5ad {input_path} -> {experiment.get_uri()}",
    )

    logging.log_io(None, f"{experiment._indent}START  READING {input_path}")

    anndata = ad.read_h5ad(input_path)

    logging.log_io(
        None,
        util.format_elapsed(s, f"{experiment._indent}FINISH READING {input_path}"),
    )

    handler_func(experiment, anndata)

    logging.log_io(
        None,
        util.format_elapsed(
            s,
            f"FINISH v1.SOMAExperiment.from_h5ad {input_path} -> {experiment.get_uri()}",
        ),
    )


# ----------------------------------------------------------------
# def from_10x(experiment: SOMAExperiment, input_path: Path) -> None:
#    """
#    Reads a 10X file and writes to a TileDB group structure.
#    """
#    s = util.get_start_stamp()
#    logging.log_io(None, f"START  v1.SOMAExperiment.from_10x {input_path} -> {experiment.get_uri()}")
#
#    logging.log_io(None, f"{experiment._indent}START  READING {input_path}")
#
#    anndata = scanpy.read_10x_h5(input_path)
#
#    logging.log_io(
#        None,
#        util.format_elapsed(s, f"{experiment._indent}FINISH READING {input_path}"),
#    )
#
#    from_anndata(experiment, anndata)
#
#    logging.log_io(
#        None,
#        util.format_elapsed(
#            s, f"FINISH v1.SOMAExperiment.from_10x {input_path} -> {experiment.get_uri()}"
#        ),
#    )


# ----------------------------------------------------------------
# TODO TEMP WIP
def from_anndata(experiment: SOMAExperiment, anndata: ad.AnnData) -> None:
    """
    Top-level writer method for creating a TileDB group for a v1.SOMAExperiment object.
    """

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
    logging.log_io(None, f"{experiment._indent}START  WRITING {experiment.get_uri()}")

    # Must be done first, to create the parent directory
    if not experiment.exists():
        experiment.create()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # OBS
    experiment.obs.from_dataframe(
        dataframe=anndata.obs, extent=256, id_column_name="obs_id"
    )
    experiment.set(experiment.obs)

    experiment.ms.create()
    experiment.set(experiment.ms)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS
    measurement = SOMAMeasurement(uri=f"{experiment.ms.get_uri()}/mRNA")
    measurement.create()
    experiment.ms.set(measurement)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/VAR
    measurement.var.from_dataframe(
        dataframe=anndata.var, extent=2048, id_column_name="var_id"
    )
    measurement.set(measurement.var)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # MS/meas/X/DATA
    measurement.X.create()
    measurement.set(measurement.X)

    Xdata = SOMASparseNdArray(uri=f"{measurement.X.get_uri()}/data")
    Xdata.from_matrix(anndata.X)
    measurement.X.set(Xdata)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # TODO: port from v0

    if len(anndata.obsm.keys()) > 0:  # do not create an empty collection
        measurement.obsm.create()
        for key in anndata.obsm.keys():
            arr = SOMADenseNdArray(uri=f"{measurement.obsm.get_uri()}/{key}")
            arr.from_matrix(anndata.obsm[key])
            measurement.obsm.set(arr)
        measurement.set(measurement.obsm)

    if len(anndata.varm.keys()) > 0:  # do not create an empty collection
        measurement.varm.create()
        for key in anndata.varm.keys():
            arr = SOMADenseNdArray(uri=f"{measurement.varm.get_uri()}/{key}")
            arr.from_matrix(anndata.varm[key])
            measurement.varm.set(arr)
        measurement.set(measurement.varm)

    if len(anndata.obsp.keys()) > 0:  # do not create an empty collection
        measurement.obsp.create()
        for key in anndata.obsp.keys():
            arr = SOMASparseNdArray(uri=f"{measurement.obsp.get_uri()}/{key}")
            arr.from_matrix(anndata.obsp[key])
            measurement.obsp.set(arr)
        measurement.set(measurement.obsp)

    if len(anndata.varp.keys()) > 0:  # do not create an empty collection
        measurement.varp.create()
        for key in anndata.varp.keys():
            arr = SOMASparseNdArray(uri=f"{measurement.varp.get_uri()}/{key}")
            arr.from_matrix(anndata.varp[key])
            measurement.varp.set(arr)
        measurement.set(measurement.varp)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # RAW
    if anndata.raw is not None:
        raw_measurement = SOMAMeasurement(uri=f"{experiment.ms.get_uri()}/raw")
        raw_measurement.create()
        experiment.ms.set(raw_measurement)

        raw_measurement.var.from_dataframe(
            dataframe=anndata.raw.var, extent=2048, id_column_name="var_id"
        )
        raw_measurement.set(raw_measurement.var)

        raw_measurement.X.create()
        raw_measurement.set(raw_measurement.X)

        rawXdata = SOMASparseNdArray(uri=f"{raw_measurement.X.get_uri()}/data")
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
            s, f"{experiment._indent}FINISH WRITING {experiment.get_uri()}"
        ),
    )


# ----------------------------------------------------------------

# def from_anndata_update_obs_and_var(experiment: SOMAExperiment, anndata: ad.AnnData) -> None:
#    """
#    Rewrites obs and var from anndata, leaving all other data in place. Useful
#    for updating schema/compression/etc. within an existing dataset.
#    """
#
#    # Without _at least_ an index, there is nothing to indicate the dimension indices.
#    if anndata.obs.index.empty or anndata.var.index.empty:
#        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")
#
#    s = util.get_start_stamp()
#    logging.log_io(None, f"{experiment._indent}START  DECATEGORICALIZING")
#
#    anndata.obs_names_make_unique()
#    anndata.var_names_make_unique()
#    anndata = util_ann._decategoricalize(anndata)
#
#    logging.log_io(
#        None,
#        util.format_elapsed(s, f"{experiment._indent}FINISH DECATEGORICALIZING"),
#    )
#
#    s = util.get_start_stamp()
#    logging.log_io(None, f"{experiment._indent}START  WRITING {experiment.get_uri()}")
#
#    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#    experiment._remove_object(experiment.obs)
#    experiment.obs.from_dataframe(dataframe=anndata.obs, extent=256)
#    experiment.set(experiment.obs)
#    tiledb.consolidate(experiment.obs.get_uri())
#    tiledb.vacuum(experiment.obs.get_uri())
#
#    experiment._remove_object(experiment.var)
#    experiment.var.from_dataframe(dataframe=anndata.var, extent=2048)
#    experiment.set(experiment.var)
#    tiledb.consolidate(experiment.var.get_uri())
#    tiledb.vacuum(experiment.var.get_uri())
#
#    logging.log_io(
#        None,
#        util.format_elapsed(s, f"{experiment._indent}FINISH WRITING {experiment.get_uri()}"),
#    )


# ----------------------------------------------------------------
def to_h5ad(experiment: SOMAExperiment, h5ad_path: Path) -> None:
    """
    Converts the experiment group to anndata format and writes it to the specified .h5ad file.
    As of 2022-05-05 this is an incomplete prototype.
    """

    s = util.get_start_stamp()
    logging.log_io(
        None, f"START  v1.SOMAExperiment.to_h5ad {experiment.get_uri()} -> {h5ad_path}"
    )

    anndata = to_anndata(experiment)

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
            s, f"FINISH v1.SOMAExperiment.to_h5ad {experiment.get_uri()} -> {h5ad_path}"
        ),
    )


# ----------------------------------------------------------------
def to_anndata(
    experiment: SOMAExperiment,
    *,
    # TODO: set a better name as capitalized-const
    # TODO: maybe if there are multiple measurements, default to the first one not named 'raw'
    measurement_name: Optional[str] = "mRNA",
) -> ad.AnnData:
    """
    Converts the experiment group to anndata. Choice of matrix formats is following
    what we often see in input .h5ad files:
    * X as scipy.sparse.csr_matrix
    * obs,var as pandas.dataframe
    * obsm,varm arrays as numpy.ndarray
    * obsp,varp arrays as scipy.sparse.csr_matrix
    """

    raise Exception("to_anndata not implemented yet")


#    s = util.get_start_stamp()
#    logging.log_io(None, f"START  v1.SOMAExperiment.to_anndata {experiment.get_uri()}")
#
#    # TODO: FINISH PORTING
#
#    # TODO: need an index-converter ... inside the class maybe?
#    # sdf.from_dataframe takes an optional id_column_name; so should sdf.to_dataframe
#
#    obs_df = experiment.obs.to_dataframe(id_column_name="obs_id")
#    var_df = experiment.ms["mRNA"].var.to_dataframe(id_column_name="var_id")
#
#    #   data = experiment.ms["mRNA"].X["data"]
#    #   assert data is not None
#    #   X_mat = data.to_csr_matrix(obs_df.index, var_df.index)
#
#    #   obsm = experiment.ms['mRNA'].obsm.to_dict_of_npnda()
#    #   varm = experiment.ms['mRNA'].varm.to_dict_of_npndacsr()
#
#    #   obsp = experiment.ms['mRNA'].obsp.to_dict_of_csr(obs_df.index, obs_df.index)
#    #   varp = experiment.ms['mRNA'].varp.to_dict_of_csr(var_df.index, var_df.index)
#
#    anndata = ad.AnnData(
#        # X=X_mat,
#        obs=obs_df,
#        var=var_df,
#        # obsm=obsm,
#        # varm=varm,
#        # obsp=obsp,
#        # varp=varp
#    )
#
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
#
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
#
#    logging.log_io(
#        None,
#        util.format_elapsed(
#            s, f"FINISH v1.SOMAExperiment.to_anndata {experiment.get_uri()}"
#        ),
#    )
#
#    return anndata


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
