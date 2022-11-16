import os
import tempfile
from pathlib import Path

import anndata
import pytest
import tiledb

import tiledbsoma
import tiledbsoma.io

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    # pbmc-small is faster for automated unit-test / CI runs.
    # input_path = HERE.parent / "anndata/pbmc3k_processed.h5ad"
    input_path = HERE.parent / "anndata/pbmc-small.h5ad"
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
    soma = tiledbsoma.SOMA(output_path)
    tiledbsoma.io.from_anndata(soma, orig)

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

    with tiledb.Group(output_path) as G:
        assert G.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "SOMA"

    # Check X/data (dense)
    with tiledb.open(os.path.join(output_path, "X", "data")) as A:
        df = A[:]
        keys = list(df.keys())
        assert keys == ["value", "obs_id", "var_id"]
        assert A.ndim == 2
        assert A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AssayMatrix"
    # Convenience accessors
    assert soma.X["data"].shape() == soma.X.data.shape()

    # Check X/raw (sparse)
    with tiledb.open(os.path.join(output_path, "raw", "X", "data")) as A:
        df = A.df[:]
        assert df.columns.to_list() == ["obs_id", "var_id", "value"]
        # verify sparsity of raw data
        assert df.shape[0] == orig.raw.X.nnz
        assert A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AssayMatrix"

    # Check obs
    with tiledb.open(os.path.join(output_path, "obs")) as A:
        df = A.df[:]
        assert df.columns.to_list() == orig.obs_keys()
        assert (
            A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
            == "AnnotationDataFrame"
        )
    assert sorted(soma.obs.ids()) == sorted(list(orig.obs_names))
    # Convenience accessors
    assert soma.obs_keys() == soma.obs_names
    assert soma.obs_names == soma.obs.ids()
    assert soma.n_obs == len(soma.obs.ids())

    # Check var
    with tiledb.open(os.path.join(output_path, "var")) as A:
        df = A.df[:]
        assert df.columns.to_list() == orig.var_keys()
        assert (
            A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
            == "AnnotationDataFrame"
        )
    assert sorted(soma.var.ids()) == sorted(list(orig.var_names))
    # Convenience accessors
    assert soma.var_keys() == soma.var_names
    assert soma.var_names == soma.var.ids()
    assert soma.n_var == len(soma.var.ids())

    # Check some annotation matrices
    # Note: pbmc3k_processed doesn't have varp.
    assert sorted(soma.obsm.keys()) == sorted(orig.obsm.keys())
    for key in orig.obsm_keys():
        with tiledb.open(os.path.join(output_path, "obsm", key)) as A:
            df = A.df[:]
            assert df.shape[0] == orig.obsm[key].shape[0]
            assert soma.obsm[key].shape() == orig.obsm[key].shape
            assert (
                A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
                == "AnnotationMatrix"
            )
    # Convenience accessors: soma.obsm.X_pca <-> soma.obsm['X_pca']
    for key in soma.obsm.keys():
        assert getattr(soma.obsm, key).shape() == soma.obsm[key].shape()

    assert sorted(soma.varm.keys()) == sorted(orig.varm.keys())
    for key in orig.varm_keys():
        with tiledb.open(os.path.join(output_path, "varm", key)) as A:
            df = A.df[:]
            assert df.shape[0] == orig.varm[key].shape[0]
            assert soma.varm[key].shape() == orig.varm[key].shape
            assert (
                A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
                == "AnnotationMatrix"
            )
    # Convenience accessors:
    for key in soma.varm.keys():
        assert getattr(soma.varm, key).shape() == soma.varm[key].shape()

    assert sorted(soma.obsp.keys()) == sorted(orig.obsp.keys())
    for key in list(orig.obsp.keys()):
        with tiledb.open(os.path.join(output_path, "obsp", key)) as A:
            df = A.df[:]
            assert df.columns.to_list() == ["obs_id_i", "obs_id_j", "value"]
            assert df.shape[0] == orig.obsp[key].nnz
            # https://github.com/single-cell-data/TileDB-SingleCell/issues/125
            # At present (without that PR's suggested enhancement) the best we
            # can get is the NNZ x attrs shape -- note that there are two
            # dims and one attr so the shape is nnz x 1.
            shape = soma.obsp[key].df().shape
            assert shape[0] == orig.obsp[key].nnz
            assert shape[1] == 1
            assert (
                A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AssayMatrix"
            )
    # Convenience accessors:
    for key in soma.obsp.keys():
        assert getattr(soma.obsp, key).shape() == soma.obsp[key].shape()

    tempdir.cleanup()


def test_export_anndata(adata):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    orig = adata

    # Ingest
    soma = tiledbsoma.SOMA(output_path)
    tiledbsoma.io.from_anndata(soma, orig)

    readback = tiledbsoma.io.to_anndata(soma)
    print("================================================================")
    print("READBACK")
    print(readback)
    print("================================================================")

    assert readback.obs.shape == orig.obs.shape
    assert readback.var.shape == orig.var.shape
    assert readback.X.shape == orig.X.shape

    for key in orig.obsm.keys():
        assert readback.obsm[key].shape == orig.obsm[key].shape
    for key in orig.varm.keys():
        assert readback.varm[key].shape == orig.varm[key].shape
    for key in orig.obsp.keys():
        assert readback.obsp[key].shape == orig.obsp[key].shape
    for key in orig.varp.keys():
        assert readback.varp[key].shape == orig.varp[key].shape


@pytest.mark.parametrize("X_capacity", [1000, 10000, 100000])
def test_X_capacity(adata, X_capacity):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    soma_options = tiledbsoma.SOMAOptions(X_capacity=X_capacity)
    soma = tiledbsoma.SOMA(output_path, soma_options=soma_options)
    tiledbsoma.io.from_anndata(soma, adata)

    with soma.X["data"]._open() as X:
        assert X.schema.capacity == X_capacity

    tempdir.cleanup()
