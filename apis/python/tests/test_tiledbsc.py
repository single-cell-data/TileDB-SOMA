import anndata
import tiledb
import tiledbsc

import pytest
import tempfile
import os
from pathlib import Path

HERE = Path(__file__).parent

@pytest.fixture
def h5ad_file(request):
    input_path = HERE.parent / "anndata/pbmc3k_processed.h5ad"
    return input_path

@pytest.fixture
def adata(h5ad_file):
    return anndata.read_h5ad(h5ad_file)

def test_import_anndata(adata):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    orig = adata

    # Ingest
    scdataset = tiledbsc.SCGroup(output_path, verbose=True)
    scdataset.from_anndata(orig)

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
        df = A.df[:]
        assert df.columns.to_list() == orig.obs_keys()

    # Check var
    with tiledb.open(os.path.join(output_path, 'var')) as A:
        df = A.df[:]
        assert df.columns.to_list() == orig.var_keys()

    # Check some annotation matrices
    # Note: pbmc3k_processed doesn't have varp.
    for key in orig.obsm_keys():
        with tiledb.open(os.path.join(output_path, 'obsm', key)) as A:
            assert A.shape == orig.obsm[key].shape

    for key in orig.varm_keys():
        with tiledb.open(os.path.join(output_path, 'varm', key)) as A:
            assert A.shape == orig.varm[key].shape

    for key in list(orig.obsp.keys()):
        with tiledb.open(os.path.join(output_path, 'obsp', key)) as A:
            assert A.shape == orig.obsp[key].shape

    tempdir.cleanup()
