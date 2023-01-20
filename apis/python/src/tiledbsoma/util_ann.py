import os
from typing import Any, List, Mapping

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import util_tiledb
from .types import Path


# ----------------------------------------------------------------
def describe_ann_file(
    input_path: Path,
    show_summary: bool = True,
    show_types: bool = False,
    show_data: bool = False,
) -> None:
    """
    This is an anndata-describer that goes a bit beyond what ``h5ls`` does for us.  In particular, it shows us that for one HDF5 file we have ``anndata.X`` being of type ``numpy.ndarray`` while for another HDF5 file we have ``anndata.X`` being of type ``scipy.sparse.csr.csr_matrix``.  This is crucial information for building I/O logic that accepts a diversity of anndata HDF5 files.
    """
    anndata = ad.read_h5ad(input_path)
    anndata.var_names_make_unique()

    print()
    print(
        f"================================================================ {input_path}"
    )

    if show_summary:
        _describe_ann_file_show_summary(anndata)
    if show_types:
        _describe_ann_file_show_types(anndata)
    if show_data:
        _describe_ann_file_show_data(anndata)


# ----------------------------------------------------------------
def _describe_ann_file_show_summary(anndata: ad.AnnData) -> None:
    print()
    print("----------------------------------------------------------------")
    print("ANNDATA SUMMARY:")
    print(anndata)

    print("X SHAPE  ", anndata.X.shape)
    print("OBS  LEN ", len(anndata.obs))
    print("VAR  LEN ", len(anndata.var))

    print("OBS IS A", type(anndata.obs))
    for name in anndata.obs.keys():
        print("  ", name, anndata.obs[name].dtype)
    print("VAR IS A", type(anndata.var))
    for name in anndata.var.keys():
        print("  ", name, anndata.var[name].dtype)

    try:  # not all groups have raw X
        print("RAW X SHAPE  ", anndata.raw.X.shape)
        print(" RAW OBS  LEN ", len(anndata.raw.X.obs_names))
        print("RAW VAR  LEN ", len(anndata.raw.X.var_names))
    except AttributeError:
        pass

    print("OBS  KEYS", list(anndata.obs.keys()))
    print("VAR  KEYS", list(anndata.var.keys()))

    print("OBSM KEYS", list(anndata.obsm.keys()))
    print("VARM KEYS", list(anndata.varm.keys()))
    print("OBSP KEYS", list(anndata.obsp.keys()))
    print("VARP KEYS", list(anndata.varp.keys()))

    _describe_ann_file_show_uns_summary(anndata.uns, ["uns"])


# ----------------------------------------------------------------
def _describe_ann_file_show_types(anndata: ad.AnnData) -> None:
    print()
    print("----------------------------------------------------------------")
    print("ANNDATA FILE TYPES:")

    namewidth = 40

    X = anndata.X
    print("%-*s %s" % (namewidth, "X/data", type(X)))
    m, n = X.shape
    print("%-*s (%d, %d)" % (namewidth, "X/data shape", m, n))
    print("%-*s %s" % (namewidth, "X/data dtype", X.dtype))
    if isinstance(X, (sp.csr_matrix, sp.csc_matrix)):
        density = X.nnz / (m * n)
        print("%-*s %.4f" % (namewidth, "X/data density", density))

    has_raw = False
    try:  # not all groups have raw X
        X = anndata.raw.X
        has_raw = True
    except AttributeError:
        pass

    if has_raw:
        X = anndata.raw.X
        print("%-*s %s" % (namewidth, "X/raw", type(X)))
        m, n = X.shape
        print("%-*s (%d, %d)" % (namewidth, "X/raw shape", m, n))
        print("%-*s %s" % (namewidth, "X/data dtype", X.dtype))
        if isinstance(X, (sp.csr_matrix, sp.csc_matrix)):
            density = X.nnz / (m * n)
            print("%-*s %.4f" % (namewidth, "X/raw density", density))

    print("%-*s %s" % (namewidth, "obs", type(anndata.obs)))
    print("%-*s %s" % (namewidth, "var", type(anndata.var)))

    for k in anndata.obsm.keys():
        print("%-*s %s" % (namewidth, "obsm/" + k, type(anndata.obsm[k])))

    for k in anndata.varm.keys():
        print("%-*s %s" % (namewidth, "varm/" + k, type(anndata.varm[k])))

    for k in anndata.obsp.keys():
        print("%-*s %s" % (namewidth, "obsp/" + k, type(anndata.obsp[k])))

    for k in anndata.varp.keys():
        print("%-*s %s" % (namewidth, "varp/" + k, type(anndata.varp[k])))

    _describe_ann_file_show_uns_types(anndata.uns, ["uns"])


# ----------------------------------------------------------------
def _describe_ann_file_show_data(anndata: ad.AnnData) -> None:
    print()
    print("----------------------------------------------------------------")
    print("ANNDATA FILE DATA:")

    print("X DATA", type(anndata.X), anndata.X.shape)
    print(anndata.X)

    print()
    print("----------------------------------------------------------------")
    print("OBS DATA")
    print(anndata.obs)

    print()
    print("----------------------------------------------------------------")
    print("VAR DATA")
    print(anndata.var)

    print()
    print("----------------------------------------------------------------")
    for k in anndata.obsm.keys():
        print()
        d = anndata.obsm[k]
        print(f"OBSM/{k} DATA", type(d), d.shape)
        print(d)

    print()
    print("----------------------------------------------------------------")
    for k in anndata.varm.keys():
        print()
        d = anndata.varm[k]
        print(f"varm/{k} DATA", type(d), d.shape)
        print(d)

    print()
    print("----------------------------------------------------------------")
    for k in anndata.obsp.keys():
        print()
        d = anndata.obsp[k]
        print(f"obsp/{k} DATA", type(d), d.shape)
        print(d)

    print()
    print("----------------------------------------------------------------")
    for k in anndata.varp.keys():
        print()
        d = anndata.varp[k]
        print(f"varp/{k} DATA", type(d), d.shape)
        print(d)

    _describe_ann_file_show_uns_data(anndata.uns, ["uns"])


# ----------------------------------------------------------------
def _describe_ann_file_show_uns_summary(
    uns: Mapping[str, Any], parent_path_components: List[str]
) -> None:
    """
    Recursively shows summary information about the anndata.uns structure.
    """
    for key in uns.keys():
        current_path_components = parent_path_components + [key]
        value = uns[key]
        if isinstance(value, Mapping):
            _describe_ann_file_show_uns_summary(value, current_path_components)
        else:
            print(os.path.sep.join(current_path_components))


# ----------------------------------------------------------------
def _describe_ann_file_show_uns_types(
    uns: Mapping[str, Any], parent_path_components: List[str]
) -> None:
    """
    Recursively shows data-type information about the anndata.uns structure.
    """
    namewidth = 40
    for key in uns.keys():
        current_path_components = parent_path_components + [key]
        value = uns[key]
        display_name = os.path.sep.join(current_path_components)
        if isinstance(value, Mapping):
            _describe_ann_file_show_uns_types(value, current_path_components)
        else:
            info = ["%-*s" % (namewidth, display_name), type(value)]
            if isinstance(value, (np.ndarray, sp.csr_matrix, pd.DataFrame)):
                info.extend((value.shape, value.dtype))
            print(*info)


# ----------------------------------------------------------------
def _describe_ann_file_show_uns_data(
    uns: Mapping[str, Any], parent_path_components: List[str]
) -> None:
    """
    Recursively shows data contained within the anndata.uns structure.
    """
    namewidth = 40
    for key in uns.keys():
        current_path_components = parent_path_components + [key]
        value = uns[key]
        display_name = os.path.sep.join(current_path_components)
        if isinstance(value, Mapping):
            _describe_ann_file_show_uns_data(value, current_path_components)
        elif (
            isinstance(value, np.ndarray)
            or isinstance(value, sp.csr_matrix)
            or isinstance(value, pd.DataFrame)
        ):
            print()
            print("%-*s" % (namewidth, display_name), type(value), value.shape)
            print(value)
        else:
            print()
            print("%-*s" % (namewidth, display_name), type(value))
            print(value)


# ----------------------------------------------------------------
def _decategoricalize_obs_or_var(obs_or_var: pd.DataFrame) -> pd.DataFrame:
    """
    Performs a typecast into types that TileDB can persist.
    """
    if len(obs_or_var.columns) > 0:
        return pd.DataFrame.from_dict(
            {
                k: util_tiledb.to_tiledb_supported_array_type(v)
                for k, v in obs_or_var.items()
            },
        )
    else:
        return obs_or_var
