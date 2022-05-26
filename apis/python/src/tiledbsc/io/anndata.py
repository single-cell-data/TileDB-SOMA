import tiledbsc
import tiledbsc.util
import tiledbsc.util_ann
import anndata as ad

# ----------------------------------------------------------------
def from_h5ad(soma: tiledbsc.SOMA, input_path: str) -> None:
    """
    Reads an .h5ad file and writes to a TileDB group structure.
    """
    if soma._verbose:
        s = tiledbsc.util.get_start_stamp()
        print(f"START  SOMA.from_h5ad {input_path} -> {soma.uri}")

    if soma._verbose:
        s = tiledbsc.util.get_start_stamp()
        print(f"{soma._indent}START  READING {input_path}")
    anndata = ad.read_h5ad(input_path)
    if soma._verbose:
        print(
            tiledbsc.util.format_elapsed(
                s, f"{soma._indent}FINISH READING {input_path}"
            )
        )

    from_anndata(soma, anndata)

    if soma._verbose:
        print(
            tiledbsc.util.format_elapsed(
                s, f"FINISH SOMA.from_h5ad {input_path} -> {soma.uri}"
            )
        )


# ----------------------------------------------------------------
def from_10x(soma: tiledbsc.SOMA, input_path: str) -> None:
    """
    Reads a 10X file and writes to a TileDB group structure.
    """
    if soma._verbose:
        s = tiledbsc.util.get_start_stamp()
        print(f"START  SOMA.from_10x {input_path} -> {soma.uri}")

    if soma._verbose:
        s = tiledbsc.util.get_start_stamp()
        print(f"{soma._indent}START  READING {input_path}")
    anndata = scanpy.read_10x_h5(input_path)
    if soma._verbose:
        print(
            tiledbsc.util.format_elapsed(
                s, f"{soma._indent}FINISH READING {input_path}"
            )
        )

    from_anndata(soma, anndata)

    if soma._verbose:
        print(
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

    if soma._verbose:
        s = tiledbsc.util.get_start_stamp()
        print(f"{soma._indent}START  DECATEGORICALIZING")
    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()
    anndata = tiledbsc.util_ann._decategoricalize(anndata)
    if soma._verbose:
        print(
            tiledbsc.util.format_elapsed(s, f"{soma._indent}FINISH DECATEGORICALIZING")
        )

    if soma._verbose:
        s = tiledbsc.util.get_start_stamp()
        print(f"{soma._indent}START  WRITING {soma.uri}")

    # Must be done first, to create the parent directory
    soma._create()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma.X.from_matrix_and_dim_values(anndata.X, anndata.obs.index, anndata.var.index)
    soma._add_object(soma.X)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma.obs.from_dataframe(dataframe=anndata.obs, extent=256)
    soma._add_object(soma.obs)

    soma.var.from_dataframe(dataframe=anndata.var, extent=2048)
    soma._add_object(soma.var)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma.obsm.from_matrices_and_dim_values(anndata.obsm, anndata.obs_names)
    soma._add_object(soma.obsm)

    soma.varm.from_matrices_and_dim_values(anndata.varm, anndata.var_names)
    soma._add_object(soma.varm)

    soma.obsp.from_matrices_and_dim_values(anndata.obsp, anndata.obs_names)
    soma._add_object(soma.obsp)

    soma.varp.from_matrices_and_dim_values(anndata.varp, anndata.var_names)
    soma._add_object(soma.varp)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if anndata.raw != None:
        soma.raw.from_anndata(anndata)
        soma._add_object(soma.raw)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if anndata.uns != None:
        soma.uns.from_anndata_uns(anndata.uns)
        soma._add_object(soma.uns)

    if soma._verbose:
        print(
            tiledbsc.util.format_elapsed(s, f"{soma._indent}FINISH WRITING {soma.uri}")
        )


# ----------------------------------------------------------------
def to_h5ad(soma: tiledbsc.SOMA, h5ad_path: str) -> None:
    """
    Converts the soma group to anndata format and writes it to the specified .h5ad file.
    As of 2022-05-05 this is an incomplete prototype.
    """

    if soma._verbose:
        s = tiledbsc.util.get_start_stamp()
        print(f"START  SOMA.to_h5ad {soma.uri} -> {h5ad_path}")

    anndata = to_anndata(soma)

    if soma._verbose:
        s2 = tiledbsc.util.get_start_stamp()
        print(f"{soma._indent}START  write {h5ad_path}")
    anndata.write_h5ad(h5ad_path)
    if soma._verbose:
        print(
            tiledbsc.util.format_elapsed(s2, f"{soma._indent}FINISH write {h5ad_path}")
        )

    if soma._verbose:
        print(
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

    if soma._verbose:
        s = tiledbsc.util.get_start_stamp()
        print(f"START  SOMA.to_anndata {soma.uri}")

    obs_df = soma.obs.to_dataframe()
    var_df = soma.var.to_dataframe()

    X_mat = soma.X.data.to_csr_matrix(obs_df.index, var_df.index)

    obsm = soma.obsm.to_dict_of_csr()
    varm = soma.varm.to_dict_of_csr()

    # TODO
    print("  OBSP OUTGEST NOT WORKING YET")
    # obsp = soma.obsp.to_dict_of_csr()
    print("  VARP OUTGEST NOT WORKING YET")
    # varp = soma.varp.to_dict_of_csr()

    (raw_X, raw_var_df, raw_varm) = soma.raw.to_anndata_raw(obs_df.index)

    anndata = ad.AnnData(X=X_mat, obs=obs_df, var=var_df, obsm=obsm, varm=varm)

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

    if soma._verbose:
        print(tiledbsc.util.format_elapsed(s, f"FINISH SOMA.to_anndata {soma.uri}"))

    return anndata


# ----------------------------------------------------------------
def to_anndata_from_raw(soma: tiledbsc.SOMA) -> ad.AnnData:
    """
    Extract only the raw parts as a new AnnData object.
    """

    obs_df = soma.obs.to_dataframe()
    var_df = soma.raw.var.to_dataframe()
    X_mat = soma.raw.X.data.to_csr_matrix(obs.index, var.index)

    return ad.AnnData(
        X=X_mat,
        obs=obs_df,
        var=var_df,
    )
