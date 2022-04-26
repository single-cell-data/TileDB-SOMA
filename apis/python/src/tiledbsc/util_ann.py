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
    anndata = ad.read_h5ad(input_path)
    anndata.var_names_make_unique()

    print()
    print("================================================================ {input_path}")
    print("ANNDATA SUMMARY:")
    print(anndata)

    print("X IS A   ", type(anndata.X))
    print("  X SHAPE  ", anndata.X.shape)
    print("  OBS  LEN ", len(anndata.obs))
    print("  VAR  LEN ", len(anndata.var))

    try: # not all groups have raw X
        print("RAW X IS A   ", type(anndata.raw.X))
        print("  X SHAPE  ", anndata.raw.X.shape)
        print("  OBS  LEN ", len(anndata.raw.X.obs_names))
        print("  VAR  LEN ", len(anndata.raw.X.var_names))
    except:
        pass

    print('OBS IS A', type(anndata.obs))
    print("  OBS  KEYS", [k for k in anndata.obs.keys()])
    print('VAR IS A', type(anndata.var))
    print("  VAR  KEYS", [k for k in anndata.var.keys()])

    print("OBSM KEYS", [k for k in anndata.obsm.keys()])
    for k in anndata.obsm.keys():
        print('  OBSM', k, 'IS A', type(anndata.obsm[k]))

    print("VARM KEYS", [k for k in anndata.varm.keys()])
    for k in anndata.varm.keys():
        print('  VARM', k, 'IS A', type(anndata.varm[k]))

    print("OBSP KEYS", [k for k in anndata.obsp.keys()])
    for k in anndata.obsp.keys():
        print('  OBSP', k, 'IS A', type(anndata.obsp[k]))

    print("VARP KEYS", [k for k in anndata.varp.keys()])
    for k in anndata.varp.keys():
        print('  VARP', k, type(anndata.varp[k]))

    # Defer unstructured data for now:
    # print("UNS  KEYS", [k for k in anndata.uns.keys()])
    # print()
    # show_uns_types(anndata.uns)

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
