import anndata as ad
import scanpy
import tiledb

import tiledbsc
import tiledbsc.util
import tiledbsc.util_ann

from .logging import logger


# ----------------------------------------------------------------
def from_h5ad(soma: tiledbsc.SOMA, input_path: str) -> None:
    """
    Reads an .h5ad file and writes to a TileDB group structure.
    """
    _from_h5ad_common(soma, input_path, from_anndata)


# ----------------------------------------------------------------
def from_h5ad_update_obs_and_var(soma: tiledbsc.SOMA, input_path: str) -> None:
    """
    Rewrites obs and var from the specified .h5ad file, leaving all other data in place. Useful for
    updating schema/compression/etc. within an existing dataset.
    """
    _from_h5ad_common(soma, input_path, from_anndata_update_obs_and_var)


# ----------------------------------------------------------------
def _from_h5ad_common(soma: tiledbsc.SOMA, input_path: str, handler_func) -> None:
    """
    Common code for things we do when processing a .h5ad file for ingest/update.
    """
    s = tiledbsc.util.get_start_stamp()
    logger.info(f"START  SOMA.from_h5ad {input_path} -> {soma.uri}")

    s = tiledbsc.util.get_start_stamp()
    logger.info(f"{soma._indent}START  READING {input_path}")

    anndata = ad.read_h5ad(input_path)

    logger.info(
        tiledbsc.util.format_elapsed(s, f"{soma._indent}FINISH READING {input_path}")
    )

    handler_func(soma, anndata)

    logger.info(
        tiledbsc.util.format_elapsed(
            s, f"FINISH SOMA.from_h5ad {input_path} -> {soma.uri}"
        )
    )


# ----------------------------------------------------------------
def from_10x(soma: tiledbsc.SOMA, input_path: str) -> None:
    """
    Reads a 10X file and writes to a TileDB group structure.
    """
    s = tiledbsc.util.get_start_stamp()
    logger.info(f"START  SOMA.from_10x {input_path} -> {soma.uri}")

    logger.info(f"{soma._indent}START  READING {input_path}")

    anndata = scanpy.read_10x_h5(input_path)

    logger.info(
        tiledbsc.util.format_elapsed(s, f"{soma._indent}FINISH READING {input_path}")
    )

    from_anndata(soma, anndata)

    logger.info(
        tiledbsc.util.format_elapsed(
            s, f"FINISH SOMA.from_10x {input_path} -> {soma.uri}"
        )
    )

    return anndata


# ----------------------------------------------------------------
def from_anndata(soma: tiledbsc.SOMA, anndata: ad.AnnData) -> None:
    """
    Top-level writer method for creating a TileDB group for a SOMA object.
    """

    # Without _at least_ an index, there is nothing to indicate the dimension indices.
    if anndata.obs.index.empty or anndata.var.index.empty:
        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

    s = tiledbsc.util.get_start_stamp()
    logger.info(f"{soma._indent}START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()
    anndata = tiledbsc.util_ann._decategoricalize(anndata)

    logger.info(
        tiledbsc.util.format_elapsed(s, f"{soma._indent}FINISH DECATEGORICALIZING")
    )

    s = tiledbsc.util.get_start_stamp()
    logger.info(f"{soma._indent}START  WRITING {soma.uri}")

    # Must be done first, to create the parent directory
    soma.create_unless_exists()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma.obs.from_dataframe(dataframe=anndata.obs, extent=256)
    soma._add_object(soma.obs)

    soma.var.from_dataframe(dataframe=anndata.var, extent=2048)
    soma._add_object(soma.var)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma.X.add_layer_from_matrix_and_dim_values(
        matrix=anndata.X,
        row_names=anndata.obs.index,
        col_names=anndata.var.index,
        layer_name="data",
    )
    soma._add_object(soma.X)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma.obsm.create_unless_exists()
    for key in anndata.obsm.keys():
        soma.obsm.add_matrix_from_matrix_and_dim_values(
            anndata.obsm[key], anndata.obs_names, key
        )
    soma._add_object(soma.obsm)

    soma.varm.create_unless_exists()
    for key in anndata.varm.keys():
        soma.varm.add_matrix_from_matrix_and_dim_values(
            anndata.varm[key], anndata.var_names, key
        )
    soma._add_object(soma.varm)

    soma.obsp.create_unless_exists()
    for key in anndata.obsp.keys():
        soma.obsp.add_matrix_from_matrix_and_dim_values(
            anndata.obsp[key], anndata.obs_names, key
        )
    soma._add_object(soma.obsp)

    soma.varp.create_unless_exists()
    for key in anndata.varp.keys():
        soma.varp.add_matrix_from_matrix_and_dim_values(
            anndata.varp[key], anndata.var_names, key
        )
    soma._add_object(soma.varp)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if anndata.raw is not None:
        soma.raw.from_anndata(anndata)
        soma._add_object(soma.raw)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if anndata.uns is not None:
        soma.uns.from_anndata_uns(anndata.uns)
        soma._add_object(soma.uns)

    logger.info(
        tiledbsc.util.format_elapsed(s, f"{soma._indent}FINISH WRITING {soma.uri}")
    )


# ----------------------------------------------------------------
def from_anndata_update_obs_and_var(soma: tiledbsc.SOMA, anndata: ad.AnnData) -> None:
    """
    Rewrites obs and var from anndata, leaving all other data in place. Useful
    for updating schema/compression/etc. within an existing dataset.
    """

    # Without _at least_ an index, there is nothing to indicate the dimension indices.
    if anndata.obs.index.empty or anndata.var.index.empty:
        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

    s = tiledbsc.util.get_start_stamp()
    logger.info(f"{soma._indent}START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()
    anndata = tiledbsc.util_ann._decategoricalize(anndata)

    logger.info(
        tiledbsc.util.format_elapsed(s, f"{soma._indent}FINISH DECATEGORICALIZING")
    )

    s = tiledbsc.util.get_start_stamp()
    logger.info(f"{soma._indent}START  WRITING {soma.uri}")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma._remove_object(soma.obs)
    soma.obs.from_dataframe(dataframe=anndata.obs, extent=256)
    soma._add_object(soma.obs)
    tiledb.consolidate(soma.obs.uri)
    tiledb.vacuum(soma.obs.uri)

    soma._remove_object(soma.var)
    soma.var.from_dataframe(dataframe=anndata.var, extent=2048)
    soma._add_object(soma.var)
    tiledb.consolidate(soma.var.uri)
    tiledb.vacuum(soma.var.uri)

    logger.info(
        tiledbsc.util.format_elapsed(s, f"{soma._indent}FINISH WRITING {soma.uri}")
    )


# ----------------------------------------------------------------
def to_h5ad(soma: tiledbsc.SOMA, h5ad_path: str) -> None:
    """
    Converts the soma group to anndata format and writes it to the specified .h5ad file.
    As of 2022-05-05 this is an incomplete prototype.
    """

    s = tiledbsc.util.get_start_stamp()
    logger.info(f"START  SOMA.to_h5ad {soma.uri} -> {h5ad_path}")

    anndata = to_anndata(soma)

    s2 = tiledbsc.util.get_start_stamp()
    logger.info(f"{soma._indent}START  write {h5ad_path}")

    anndata.write_h5ad(h5ad_path)

    logger.info(
        tiledbsc.util.format_elapsed(s2, f"{soma._indent}FINISH write {h5ad_path}")
    )

    logger.info(
        tiledbsc.util.format_elapsed(
            s, f"FINISH SOMA.to_h5ad {soma.uri} -> {h5ad_path}"
        )
    )


# ----------------------------------------------------------------
def to_anndata(soma: tiledbsc.SOMA) -> ad.AnnData:
    """
    Converts the soma group to anndata. Choice of matrix formats is following
    what we often see in input .h5ad files:
    * X as scipy.sparse.csr_matrix
    * obs,var as pandas.dataframe
    * obsm,varm arrays as numpy.ndarray
    * obsp,varp arrays as scipy.sparse.csr_matrix
    As of 2022-05-05 this is an incomplete prototype.
    """

    s = tiledbsc.util.get_start_stamp()
    logger.info(f"START  SOMA.to_anndata {soma.uri}")

    obs_df = soma.obs.df()
    var_df = soma.var.df()

    X_mat = soma.X["data"].to_csr_matrix(obs_df.index, var_df.index)

    obsm = soma.obsm.to_dict_of_csr()
    varm = soma.varm.to_dict_of_csr()

    obsp = soma.obsp.to_dict_of_csr(obs_df.index, obs_df.index)
    varp = soma.varp.to_dict_of_csr(var_df.index, var_df.index)

    anndata = ad.AnnData(
        X=X_mat, obs=obs_df, var=var_df, obsm=obsm, varm=varm, obsp=obsp, varp=varp
    )

    raw = None
    if soma.raw.exists():
        (raw_X, raw_var_df, raw_varm) = soma.raw.to_anndata_raw(obs_df.index)
        raw = ad.Raw(
            anndata,
            X=raw_X,
            var=raw_var_df,
            varm=raw_varm,
        )

    uns = soma.uns.to_dict_of_matrices()

    anndata = ad.AnnData(
        X=anndata.X,
        dtype=None if anndata.X is None else anndata.X.dtype,  # some datasets have no X
        obs=anndata.obs,
        var=anndata.var,
        obsm=anndata.obsm,
        obsp=anndata.obsp,
        varm=anndata.varm,
        varp=anndata.varp,
        raw=raw,
        uns=uns,
    )

    logger.info(tiledbsc.util.format_elapsed(s, f"FINISH SOMA.to_anndata {soma.uri}"))

    return anndata


# ----------------------------------------------------------------
def to_anndata_from_raw(soma: tiledbsc.SOMA) -> ad.AnnData:
    """
    Extract only the raw parts as a new AnnData object.
    """

    obs_df = soma.obs.df()
    var_df = soma.raw.var.df()
    X_mat = soma.raw.X["data"].to_csr_matrix(obs_df.index, var_df.index)

    return ad.AnnData(
        X=X_mat,
        obs=obs_df,
        var=var_df,
    )
