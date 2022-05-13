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
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
    return input_path

def test_soma_group_indexing(h5ad_file):
    """
    Verify basic group-member access at the tiledbsc-py level.
    """

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    soma = tiledbsc.SOMA(output_path, verbose=False)
    soma.from_h5ad(h5ad_file)

    # Structure:
    #   X/data
    #   obs
    #   var
    #   obsm/X_pca
    #   obsm/X_tsne
    #   obsm/X_umap
    #   obsm/X_draw_graph_fr
    #   varm/PCs
    #   obsp/distances
    #   obsp/connectivities
    #   raw/X/data
    #   raw/var
    #   raw/varm/PCs

    assert set(soma.get_member_names()) == set(['uns', 'varm', 'X', 'raw', 'obsp', 'varp', 'var', 'obsm', 'obs'])
    assert set(soma.X.get_member_names()) == set(['data'])

    assert set(soma.obsm.get_member_names()) == set(['X_pca', 'X_tsne'])
    assert isinstance(soma.obsm['X_pca'], tiledbsc.AnnotationMatrix)
    assert soma.obsm['nonesuch'] is None

    assert set(soma.varm.get_member_names()) == set(['PCs'])
    assert isinstance(soma.varm['PCs'], tiledbsc.AnnotationMatrix)
    assert soma.varm['nonesuch'] is None

    assert set(soma.obsp.get_member_names()) == set(['distances'])
    assert isinstance(soma.obsp['distances'], tiledbsc.AnnotationPairwiseMatrix)
    assert soma.varp['nonesuch'] is None

    assert set(soma.uns.get_member_names()) == set(['neighbors'])
    assert isinstance(soma.uns['neighbors'], tiledbsc.UnsGroup)
    assert set(soma.uns['neighbors'].get_member_names()) == set(['params'])
    assert isinstance(soma.uns['neighbors']['params'], tiledbsc.UnsGroup)
    assert set(soma.uns['neighbors']['params'].get_member_names()) == set(['method'])
    assert isinstance(soma.uns['neighbors']['params']['method'], tiledbsc.UnsArray)
    assert soma.uns['nonesuch'] is None
