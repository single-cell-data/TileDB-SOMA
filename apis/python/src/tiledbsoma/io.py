from concurrent.futures import ThreadPoolExecutor
from typing import Callable

import anndata as ad
import scanpy
import tiledb

import tiledbsoma
import tiledbsoma.logging
import tiledbsoma.util
import tiledbsoma.util_ann

from .logging import log_io
from .types import Path


# ----------------------------------------------------------------
def from_h5ad_unless_exists(
    soma: tiledbsoma.SOMA, input_path: str, X_layer_name: str = "data"
) -> None:
    """
    Skips the ingest if the SOMA is already there. A convenient keystroke-saver
    so users don't need to replicate the if-test.
    """
    if tiledbsoma.util.is_soma(soma.uri, ctx=soma._ctx):
        tiledbsoma.logging.logger.info(
            f"Already exists, skipping ingest: {soma.nested_name}"
        )
    else:
        from_h5ad(soma, input_path, X_layer_name)


# ----------------------------------------------------------------
def from_h5ad(
    soma: tiledbsoma.SOMA,
    input_path: Path,
    X_layer_name: str = "data",
    ingest_mode: str = "write",
) -> None:
    """
    Reads an ``.h5ad`` local-disk file and writes to a TileDB SOMA structure.

    The "write" ingest_mode (which is the default) writes all data, creating new layers if the soma already exists.

    The "resume" ingest_mode skips data writes if data are within dimension ranges of the existing soma.
    This is useful for continuing after a partial/interrupted previous upload.

    The "schema_only" ingest_mode creates groups and array schema, without writing array data.
    This is useful as a prep-step for parallel append-ingest of multiple H5ADs to a single soma.
    """
    assert ingest_mode in tiledbsoma.util.INGEST_MODES

    if isinstance(input_path, ad.AnnData):
        raise Exception("Input path is an AnnData object -- did you want from_anndata?")
    assert ingest_mode in ["write", "schema_only", "resume"]
    _from_h5ad_common(soma, input_path, _from_anndata_aux, X_layer_name, ingest_mode)


# ----------------------------------------------------------------
def from_h5ad_update_obs_and_var(soma: tiledbsoma.SOMA, input_path: Path) -> None:
    """
    Rewrites obs and var from the specified .h5ad file, leaving all other data in place. Useful for
    updating schema/compression/etc. within an existing dataset.
    """
    _from_h5ad_common(
        soma,
        input_path,
        from_anndata_update_obs_and_var,
        "unused",
        "write",
    )


# ----------------------------------------------------------------
def _from_h5ad_common(
    soma: tiledbsoma.SOMA,
    input_path: Path,
    handler_func: Callable[[tiledbsoma.SOMA, ad.AnnData, str, str], None],
    X_layer_name: str,
    ingest_mode: str,
) -> None:
    """
    Common code for things we do when processing a .h5ad file for ingest/update.
    """
    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"START  SOMA.from_h5ad {input_path} -> {soma.nested_name}")

    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"{soma._indent}START  READING {input_path}")

    # No gains found from using "r+" AFAICT
    anndata = ad.read_h5ad(input_path, backed="r")

    log_io(
        None,
        tiledbsoma.util.format_elapsed(s, f"{soma._indent}FINISH READING {input_path}"),
    )

    handler_func(soma, anndata, X_layer_name, ingest_mode)

    log_io(
        None,
        tiledbsoma.util.format_elapsed(
            s, f"FINISH SOMA.from_h5ad {input_path} -> {soma.nested_name}"
        ),
    )


# ----------------------------------------------------------------
def from_10x_unless_exists(soma: tiledbsoma.SOMA, input_path: Path) -> None:
    """
    Skips the ingest if the SOMA is already there. A convenient keystroke-saver
    so users don't need to replicate the if-test.
    """
    if tiledbsoma.util.is_soma(soma.uri):
        tiledbsoma.logging.logger.info(
            f"Already exists, skipping ingest: {soma.nested_name}"
        )
    else:
        from_10x(soma, input_path)


def from_10x(
    soma: tiledbsoma.SOMA,
    input_path: Path,
    X_layer_name: str = "data",
    ingest_mode: str = "write",
) -> None:
    """
    Reads a 10X file and writes to a TileDB group structure.
    """
    assert ingest_mode in tiledbsoma.util.INGEST_MODES

    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"START  SOMA.from_10x {input_path} -> {soma.nested_name}")

    log_io(None, f"{soma._indent}START  READING {input_path}")

    anndata = scanpy.read_10x_h5(input_path)

    log_io(
        None,
        tiledbsoma.util.format_elapsed(s, f"{soma._indent}FINISH READING {input_path}"),
    )

    _from_anndata_aux(soma, anndata, X_layer_name, ingest_mode=ingest_mode)

    log_io(
        None,
        tiledbsoma.util.format_elapsed(
            s, f"FINISH SOMA.from_10x {input_path} -> {soma.nested_name}"
        ),
    )


# ----------------------------------------------------------------
def from_anndata_unless_exists(
    soma: tiledbsoma.SOMA,
    anndata: ad.AnnData,
    X_layer_name: str = "data",
    ingest_mode: str = "write",
) -> None:
    """
    Skips the ingest if the SOMA is already there. A convenient keystroke-saver
    so users don't need to replicate the if-test.
    """
    assert ingest_mode in tiledbsoma.util.INGEST_MODES

    assert ingest_mode in ["write", "schema_only", "resume"]
    if tiledbsoma.util.is_soma(soma.uri):
        tiledbsoma.logging.logger.info(
            f"Already exists, skipping ingest: {soma.nested_name}"
        )
    else:
        _from_anndata_aux(soma, anndata, X_layer_name, ingest_mode)


# ----------------------------------------------------------------
def from_anndata(
    soma: tiledbsoma.SOMA,
    anndata: ad.AnnData,
    X_layer_name: str = "data",
    ingest_mode: str = "write",
) -> None:
    """
    Given an in-memory ``AnnData`` object, writes to a TileDB SOMA structure.

    See comments on `from_h5ad` for details about the `ingest_mode` parameter.
    """
    assert ingest_mode in tiledbsoma.util.INGEST_MODES

    assert ingest_mode in ["write", "schema_only", "resume"]
    return _from_anndata_aux(soma, anndata, X_layer_name, ingest_mode)


def _from_anndata_aux(
    soma: tiledbsoma.SOMA,
    anndata: ad.AnnData,
    X_layer_name: str,
    ingest_mode: str,
) -> None:
    """
    Helper method for ``from_anndata``. This simplified type-checking using ``mypy`` with regard to
    callback functions -- this helper method as ``X_layer_name`` as non-optional, which confuses
    ``mypy`` less.
    """
    if not isinstance(anndata, ad.AnnData):
        raise Exception(
            "Second argument is not an AnnData object -- did you want from_h5ad?"
        )

    # Without _at least_ an index, there is nothing to indicate the dimension indices.
    if anndata.obs.index.empty or anndata.var.index.empty:
        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"{soma._indent}START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()

    log_io(
        None,
        tiledbsoma.util.format_elapsed(s, f"{soma._indent}FINISH DECATEGORICALIZING"),
    )

    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"{soma._indent}START  WRITING {soma.nested_name}")

    # Must be done first, to create the parent directory
    soma.create_unless_exists()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    max_workers = soma._soma_options.max_thread_pool_workers

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    futures = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures.append(
            executor.submit(
                soma.obs.from_dataframe,
                dataframe=tiledbsoma.util_ann._decategoricalize_obs_or_var(anndata.obs),
                extent=256,
                ingest_mode=ingest_mode,
            )
        )
        futures.append(
            executor.submit(
                soma.var.from_dataframe,
                dataframe=tiledbsoma.util_ann._decategoricalize_obs_or_var(anndata.var),
                extent=2048,
                ingest_mode=ingest_mode,
            )
        )
        futures.append(
            executor.submit(
                soma.X.add_layer_from_matrix_and_dim_values,
                # Note:
                # * Using anndata.X[:] loads the entire CSR (or whatever) into memory
                # * If we use just anndata.X then we get a super-cool object which
                #   we can index with things like matrix[0:10], matrix[10:20], etc
                #   which would allow us even further control over memory use --
                #   in particular this would allow us to do distributed ingest
                #   via UDFs ...
                # * ... but ... sadly there's a problem (in main-old).
                # * Context: while in main the dims are join-ids and we can H5AD-read
                #   and TileDB-write the matrices in their H5AD order (e.g. matrix[0:10] or
                #   whatever), in main-old we must write string dims in sorted order
                #   in order to avoid non-overlapping fragments. That in turn requires
                #   util._get_sort_and_permutation() which gets matrix[i].nnz. And
                #   in turn the anndata super-snazzy indexed-matrix object I was raving
                #   about ... doesn't do well _at all_ in performance here. Effectively
                #   the util._get_sort_and_permutation() never completes. :(
                # * We can sub-sample -- just ask at most say 100 rows for their nnz and
                #   extrapolate that up -- which works fine. But there's a deeper issue
                #   which is (from h5py):
                #     TypeError: Indexing elements must be in increasing order
                #   i.e. h5py is not designed to do these permuted reads.
                # * So for main we'll be able to drop the [:] here and it will be grand. :)
                matrix=anndata.X[:],
                row_names=anndata.obs.index,
                col_names=anndata.var.index,
                layer_name=X_layer_name,
                ingest_mode=ingest_mode,
            )
        )

        if anndata.layers is not None:
            for layer_name in anndata.layers:
                futures.append(
                    executor.submit(
                        soma.X.add_layer_from_matrix_and_dim_values,
                        # See comments on anndata.X[:], above
                        matrix=anndata.layers[layer_name][:],
                        row_names=anndata.obs.index,
                        col_names=anndata.var.index,
                        layer_name=layer_name,
                        ingest_mode=ingest_mode,
                    )
                )

    for future in futures:
        future.result()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma._add_object(soma.obs)
    soma._add_object(soma.var)
    soma._add_object(soma.X)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma.obsm.create_unless_exists()
    soma.varm.create_unless_exists()
    soma.obsp.create_unless_exists()
    soma.varp.create_unless_exists()

    futures = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for key in anndata.obsm.keys():
            futures.append(
                executor.submit(
                    soma.obsm.add_matrix_from_matrix_and_dim_values,
                    tiledbsoma.util._to_tiledb_supported_array_type(anndata.obsm[key]),
                    anndata.obs_names,
                    key,
                    ingest_mode=ingest_mode,
                )
            )

        for key in anndata.varm.keys():
            futures.append(
                executor.submit(
                    soma.varm.add_matrix_from_matrix_and_dim_values,
                    tiledbsoma.util._to_tiledb_supported_array_type(anndata.varm[key]),
                    anndata.var_names,
                    key,
                    ingest_mode=ingest_mode,
                )
            )

        for key in anndata.obsp.keys():
            futures.append(
                executor.submit(
                    soma.obsp.add_matrix_from_matrix_and_dim_values,
                    tiledbsoma.util._to_tiledb_supported_array_type(anndata.obsp[key]),
                    anndata.obs_names,
                    key,
                    ingest_mode=ingest_mode,
                )
            )

        for key in anndata.varp.keys():
            futures.append(
                executor.submit(
                    soma.varp.add_matrix_from_matrix_and_dim_values,
                    tiledbsoma.util._to_tiledb_supported_array_type(anndata.varp[key]),
                    anndata.var_names,
                    key,
                    ingest_mode=ingest_mode,
                )
            )

        if anndata.raw is not None:
            futures.append(
                executor.submit(
                    soma.raw.from_anndata,
                    anndata,
                    ingest_mode=ingest_mode,
                )
            )

    for future in futures:
        future.result()

    soma._add_object(soma.obsm)
    soma._add_object(soma.varm)
    soma._add_object(soma.obsp)
    soma._add_object(soma.varp)
    if anndata.raw is not None:
        soma._add_object(soma.raw)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Already parallelized recursively
    if ingest_mode != "schema_only":
        # Writing multiple H5ADs in append mode to the same SOMA is a supported mode.  However the
        # uns structures _cannot_ have all the same schema -- in particular there are dense arrays.
        # For append mode, users must set `anndata.uns = {}`, or "nest" each input anndata object's
        # `uns` as `anndata.uns = { "some_unique_name" : anndata.uns }`. In either case, there is
        # nothing to be done at the schema_only step. The uns objects _have_ no fixed schema -- as
        # indicated by the name `uns` for "unstructured".
        if anndata.uns is not None:
            soma.uns.from_anndata_uns(anndata.uns, ingest_mode)
            soma._add_object(soma.uns)

    log_io(
        f"Wrote {soma.nested_name}",
        tiledbsoma.util.format_elapsed(
            s, f"{soma._indent}FINISH WRITING {soma.nested_name}"
        ),
    )


# ----------------------------------------------------------------
def from_anndata_update_obs_and_var(
    soma: tiledbsoma.SOMA,
    anndata: ad.AnnData,
    _unused1: str = "",
    _unused2: str = "",
) -> None:
    """
    Rewrites obs and var from anndata, leaving all other data in place. Useful
    for updating schema/compression/etc. within an existing dataset.
    """

    # Without _at least_ an index, there is nothing to indicate the dimension indices.
    if anndata.obs.index.empty or anndata.var.index.empty:
        raise NotImplementedError("Empty AnnData.obs or AnnData.var unsupported.")

    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"{soma._indent}START  DECATEGORICALIZING")

    anndata.obs_names_make_unique()
    anndata.var_names_make_unique()

    obs = tiledbsoma.util_ann._decategoricalize_obs_or_var(anndata.obs)
    var = tiledbsoma.util_ann._decategoricalize_obs_or_var(anndata.var)

    log_io(
        None,
        tiledbsoma.util.format_elapsed(s, f"{soma._indent}FINISH DECATEGORICALIZING"),
    )

    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"{soma._indent}START  WRITING {soma.nested_name}")

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    soma._remove_object(soma.obs)
    soma.obs.from_dataframe(dataframe=obs, extent=256)
    soma._add_object(soma.obs)
    tiledb.consolidate(soma.obs.uri)
    tiledb.vacuum(soma.obs.uri)

    soma._remove_object(soma.var)
    soma.var.from_dataframe(dataframe=var, extent=2048)
    soma._add_object(soma.var)
    tiledb.consolidate(soma.var.uri)
    tiledb.vacuum(soma.var.uri)

    log_io(
        None,
        tiledbsoma.util.format_elapsed(
            s, f"{soma._indent}FINISH WRITING {soma.nested_name}"
        ),
    )


# ----------------------------------------------------------------
def to_h5ad(soma: tiledbsoma.SOMA, h5ad_path: Path, X_layer_name: str = "data") -> None:
    """
    Converts the soma group to anndata format and writes it to the specified .h5ad file.
    As of 2022-05-05 this is an incomplete prototype.
    """

    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"START  SOMA.to_h5ad {soma.nested_name} -> {h5ad_path}")

    anndata = to_anndata(soma, X_layer_name=X_layer_name)

    s2 = tiledbsoma.util.get_start_stamp()
    log_io(None, f"{soma._indent}START  write {h5ad_path}")

    anndata.write_h5ad(h5ad_path)

    log_io(
        None,
        tiledbsoma.util.format_elapsed(s2, f"{soma._indent}FINISH write {h5ad_path}"),
    )

    log_io(
        None,
        tiledbsoma.util.format_elapsed(
            s, f"FINISH SOMA.to_h5ad {soma.nested_name} -> {h5ad_path}"
        ),
    )


# ----------------------------------------------------------------
def to_anndata(soma: tiledbsoma.SOMA, X_layer_name: str = "data") -> ad.AnnData:
    """
    Converts the soma group to anndata. Choice of matrix formats is following what we often see in input ``.h5ad`` files:

    * X as ``scipy.sparse.csr_matrix``
    * ``obs``, ``var`` as ``pandas.dataframe``
    * ``obsm``, ``varm`` arrays as ``numpy.ndarray``
    * ``obsp``, ``varp`` arrays as ``scipy.sparse.csr_matrix``
    """

    s = tiledbsoma.util.get_start_stamp()
    log_io(None, f"START  SOMA.to_anndata {soma.nested_name}")

    obs_df = soma.obs.df()
    var_df = soma.var.df()

    data = soma.X[X_layer_name]
    assert data is not None
    X_mat = data.to_csr_matrix(obs_df.index, var_df.index)

    obsm = soma.obsm.to_dict_of_csr()
    varm = soma.varm.to_dict_of_csr()

    obsp = soma.obsp.to_dict_of_csr(obs_df.index, obs_df.index)
    varp = soma.varp.to_dict_of_csr(var_df.index, var_df.index)

    anndata = ad.AnnData(
        X=X_mat,
        obs=obs_df,
        var=var_df,
        obsm=obsm,
        varm=varm,
        obsp=obsp,
        varp=varp,
        dtype=None if X_mat is None else X_mat.dtype,  # some datasets have no X
    )

    raw = None
    if soma.raw.exists():
        (raw_X, raw_var_df, raw_varm) = soma.raw.to_anndata_raw(
            obs_df.index, X_layer_name=X_layer_name
        )
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

    log_io(
        None,
        tiledbsoma.util.format_elapsed(s, f"FINISH SOMA.to_anndata {soma.nested_name}"),
    )

    return anndata


# ----------------------------------------------------------------
def to_anndata_from_raw(
    soma: tiledbsoma.SOMA, X_layer_name: str = "data"
) -> ad.AnnData:
    """
    Extract only the raw parts as a new AnnData object.
    """

    obs_df = soma.obs.df()
    var_df = soma.raw.var.df()
    data = soma.raw.X[X_layer_name]
    assert data is not None
    X_mat = data.to_csr_matrix(obs_df.index, var_df.index)

    return ad.AnnData(X=X_mat, obs=obs_df, var=var_df)
