import sys, os
import math
import anndata as ad
import pandas as pd
import numpy as np
import scipy
import tiledb

# ----------------------------------------------------------------
def find_csr_chunk_size(mat: scipy.sparse._csr.csr_matrix, start_row_index: int, capacity: int):
    """
    Given a CSR matrix and a start row index, returns the number of rows with cumulative nnz
    targeted to be around the capacity argument. Context is chunked-COO ingest of larger CSR
    matrices: if mat is say 8000x9000 but sparse, maybe we'll read rows 0:45 as one chunk and
    convert that to COO and ingest, then maybe rows 46:78 as a second chunk and convert that to COO
    and ingest, and so on.
    :param mat: The input CSR matrix.
    :param start_row_index: the row index at which to start a chunk.
    :param capacity: TileDB array-schema capacity parameter.
    """
    chunk_size = 1
    sum_nnz = 0
    for row_index in range(start_row_index, mat.shape[0]):
        sum_nnz += mat[row_index].nnz
        if sum_nnz > capacity:
            break
        chunk_size += 1

    return chunk_size

# ----------------------------------------------------------------
def describe_ann_file(input_path: str, types_only=False):
    """
    This is an anndata-describer that goes a bit beyond what h5ls does for us.
    In particular, it shows us that for one HDF5 file we have anndata.X being of type numpy.ndarray
    while for another HDF5 file we have anndata.X being of type 'scipy.sparse.csr.csr_matrix'.  This is
    crucial information for building I/O logic that accepts a diversity of anndata HDF5 files.
    """
    anndata = ad.read_h5ad(input_path)
    anndata.var_names_make_unique()

    print()
    print(f"================================================================ {input_path}")
    print("ANNDATA FILE TYPES:")

    namewidth = 30

    X = anndata.X
    print("%-*s %s" % (namewidth, "X/data", type(X)))
    m,n = X.shape
    print("%-*s (%d, %d)" % (namewidth, "X/data shape", m, n))
    if isinstance(X, scipy.sparse._csr.csr_matrix) or isinstance(X, scipy.sparse._csc.csc_matrix):
        density = X.nnz / (m*n)
        print("%-*s %.4f" % (namewidth, "X/data density", density))

    has_raw = False
    try: # not all groups have raw X
        X = anndata.raw.X
        has_raw = True
    except:
        pass

    if has_raw:
        X = anndata.raw.X
        print("%-*s %s" % (namewidth, "X/raw", type(X)))
        m,n = X.shape
        print("%-*s (%d, %d)" % (namewidth, "X/raw shape", m, n))
        if isinstance(X, scipy.sparse._csr.csr_matrix) or isinstance(X, scipy.sparse._csc.csc_matrix):
            density = X.nnz / (m*n)
            print("%-*s %.4f" % (namewidth, "X/raw density", density))

    print("%-*s %s" % (namewidth, "obs", type(anndata.obs)))
    print("%-*s %s" % (namewidth, "var", type(anndata.var)))

    for k in anndata.obsm.keys():
        print("%-*s %s" % (namewidth, "obsm/"+k, type(anndata.obsm[k])))

    for k in anndata.varm.keys():
        print("%-*s %s" % (namewidth, "varm/"+k, type(anndata.varm[k])))

    for k in anndata.obsp.keys():
        print("%-*s %s" % (namewidth, "obsp/"+k, type(anndata.obsp[k])))

    for k in anndata.varp.keys():
        print("%-*s %s" % (namewidth, "varp/"+k, type(anndata.varp[k])))

    if types_only:
        return

    print()
    print("ANNDATA SUMMARY:")
    print(anndata)

    print("X SHAPE  ", anndata.X.shape)
    print("OBS  LEN ", len(anndata.obs))
    print("VAR  LEN ", len(anndata.var))

    print('OBS IS A', type(anndata.obs))
    for name in anndata.obs.keys():
        print("  ", name, anndata.obs[name].dtype)
    print('VAR IS A', type(anndata.var))
    for name in anndata.var.keys():
        print("  ", name, anndata.var[name].dtype)

    try: # not all groups have raw X
        print("RAW X SHAPE  ", anndata.raw.X.shape)
        print(" RAW OBS  LEN ", len(anndata.raw.X.obs_names))
        print("RAW VAR  LEN ", len(anndata.raw.X.var_names))
    except:
        pass

    print("OBS  KEYS", list(anndata.obs.keys()))
    print("VAR  KEYS", list(anndata.var.keys()))

    print("OBSM KEYS", list(anndata.obsm.keys()))
    print("VARM KEYS", list(anndata.varm.keys()))
    print("OBSP KEYS", list(anndata.obsp.keys()))
    print("VARP KEYS", list(anndata.varp.keys()))

    # Defer unstructured data for now:
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
