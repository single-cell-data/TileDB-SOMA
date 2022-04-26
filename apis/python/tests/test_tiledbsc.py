import anndata
import tiledb
import tiledbsc

import pytest
import tempfile
import os

def test_import_anndata(request):
    # Make sure this works regardless of from what directory level the `python -m pytest ...` is invoked
    ourdir = request.fspath.dirname

    # Set up anndata input path and tiledb-group output path
    input_path = os.path.join(ourdir, '..', 'anndata', 'pbmc3k_processed.h5ad')
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    scdataset = tiledbsc.SCGroup(output_path, verbose=True)
    scdataset.from_h5ad(input_path)

    # Structure:
    #   X/data
    #   X/raw
    #   obs
    #   var
    #   obsm/X_pca
    #   obsm/X_tsne
    #   obsm/X_umap
    #   obsm/X_draw_graph_fr
    #   varm/PCs
    #   obsp/distances
    #   obsp/connectivities

    # Check X/data

    # Note: intentionally avoiding syntax `with tiledb.open(...) as A:`
    # to make this easier to interact with in the Python interpreter.
    A = tiledb.open(os.path.join(output_path, 'X', 'data'))
    df = A[:]
    keys = [k for k in df.keys()]
    assert keys == ['data', 'obs_id', 'var_id']
    assert A.ndim == 2
    A.close()

    # Check obs
    A = tiledb.open(os.path.join(output_path, 'obs'))
    df = A[:]
    keys = [k for k in df.keys()]
    assert keys == ['n_genes', 'percent_mito', 'n_counts', 'louvain', 'index']

    # Check var
    A = tiledb.open(os.path.join(output_path, 'var'))
    df = A[:]
    keys = [k for k in df.keys()]
    assert keys == ['n_cells', 'index']
    A.close()

    # Check some annotation matrices
    # Note: pbmc3k_processed doesn't have varp.

    A = tiledb.open(os.path.join(output_path, 'obsm', 'X_pca'))
    assert A.shape == (2638, 50)
    A.close()

    A = tiledb.open(os.path.join(output_path, 'varm', 'PCs'))
    assert A.shape == (1838, 50)
    A.close()

    A = tiledb.open(os.path.join(output_path, 'obsp', 'connectivities'))
    assert A.shape == (2638, 2638)
    A.close()

    tempdir.cleanup()
