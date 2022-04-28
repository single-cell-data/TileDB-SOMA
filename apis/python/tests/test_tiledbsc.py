import anndata
import tiledb
import tiledbsc

import pytest
import tempfile
import os

@pytest.fixture
def h5ad_file(request):
    # Make sure this works regardless of from what directory level the `python -m pytest ...` is invoked
    ourdir = request.fspath.dirname
    input_path = os.path.join(ourdir, '..', 'anndata', 'pbmc3k_processed.h5ad')
    return input_path

@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)

def test_import_anndata(h5ad_file):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    scdataset = tiledbsc.SCGroup(output_path, verbose=True)
    scdataset.from_h5ad(h5ad_file)

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
    with tiledb.open(os.path.join(output_path, 'X', 'data')) as A:
        df = A[:]
        keys = list(df.keys())
        assert keys == ['data', 'obs_id', 'var_id']
        assert A.ndim == 2

    # Check obs
    with tiledb.open(os.path.join(output_path, 'obs')) as A:
        df = A[:]
        keys = list(df.keys())
        assert keys == ['n_genes', 'percent_mito', 'n_counts', 'louvain', 'index']

    # Check var
    with tiledb.open(os.path.join(output_path, 'var')) as A:
        df = A[:]
        keys = list(df.keys())
        assert keys == ['n_cells', 'index']

    # Check some annotation matrices
    # Note: pbmc3k_processed doesn't have varp.

    with tiledb.open(os.path.join(output_path, 'obsm', 'X_pca')) as A:
        assert A.shape == (2638, 50)

    with tiledb.open(os.path.join(output_path, 'varm', 'PCs')) as A:
        assert A.shape == (1838, 50)

    with tiledb.open(os.path.join(output_path, 'obsp', 'connectivities')) as A:
        assert A.shape == (2638, 2638)

    tempdir.cleanup()
