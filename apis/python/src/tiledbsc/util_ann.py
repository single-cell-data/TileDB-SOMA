import sys, os
import anndata as ad
import pandas as pd
import numpy as np
import tiledb

# ----------------------------------------------------------------
def describe_ann_file(input_path):
    """
    This is an anndata-describer that goes a bit beyond what h5ls does for us.
    In particular, it shows us that for one HDF5 file we have anndata.X being of type numpy.ndarray
    while for another HDF5 file we have anndata.X being of type 'scipy.sparse.csr.csr_matrix'.  This is
    crucial information for building I/O logic that accepts a diversity of anndata HDF5 files.
    """
    h5ad_data = ad.read_h5ad(input_path)
    h5ad_data.var_names_make_unique()

    print()
    print("================================================================ {input_path}")
    print("ANNDATA SUMMARY:")
    print(h5ad_data)

    print("X IS A   ", type(h5ad_data.X))
    print("  X SHAPE  ", h5ad_data.X.shape)
    print("  OBS  LEN ", len(h5ad_data.obs))
    print("  VAR  LEN ", len(h5ad_data.var))

    print('OBS IS A', type(h5ad_data.obs))
    print("  OBS  KEYS", [k for k in h5ad_data.obs.keys()])
    print('VAR IS A', type(h5ad_data.var))
    print("  VAR  KEYS", [k for k in h5ad_data.var.keys()])

    print("OBSM KEYS", [k for k in h5ad_data.obsm.keys()])
    for k in h5ad_data.obsm.keys():
        print('  OBSM', k, 'IS A', type(h5ad_data.obsm[k]))

    print("VARM KEYS", [k for k in h5ad_data.varm.keys()])
    for k in h5ad_data.varm.keys():
        print('  VARM', k, 'IS A', type(h5ad_data.varm[k]))

    print("OBSP KEYS", [k for k in h5ad_data.obsp.keys()])
    for k in h5ad_data.obsp.keys():
        print('  OBSP', k, 'IS A', type(h5ad_data.obsp[k]))

    print("VARP KEYS", [k for k in h5ad_data.varp.keys()])
    for k in h5ad_data.varp.keys():
        print('  VARP', k, type(h5ad_data.varp[k]))

    # Defer unstructured data for now:
    # print("UNS  KEYS", [k for k in h5ad_data.uns.keys()])
    # print()
    # show_uns_types(h5ad_data.uns)

# ----------------------------------------------------------------
def show_uns_types(uns, depth=0):
    """
    Recursive helper function for describe_ann_file, given that `uns` data
    can be arbitrarily nested.
    """
    leader = "  " * depth
    for k in uns.keys():
        v = uns[k]
        if isinstance(v, np.ndarray):
            print(leader, 'UNS', k, "IS A", type(uns[k]))
        elif isinstance(v, dict) or isinstance(v, ad.compat._overloaded_dict.OverloadedDict):
            print(leader, 'UNS', k, "IS A", type(uns[k]))
            show_uns_types(v, depth+1)
        else:
            print(leader, 'UNS', k, "IS A", type(uns[k]), "which is unrecognized")
