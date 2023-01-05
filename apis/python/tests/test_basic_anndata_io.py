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


@pytest.mark.parametrize(
    "ingest_modes",
    [
        ["write"],  # Standard ingest: normal use-case
        ["schema_only"],  # User only creates schema
        ["write", "schema_only"],  # User creates schema, then writes data
        [
            "write",
            "schema_only",
            "schema_only",
        ],  # User creates schema, then writes and re-writes the same data
        ["schema_only", "schema_only"],  # User writes and re-writes the same data
        [
            "write",
            "resume",
        ],  # User writes data, then a subsequent write creates nothing new
        ["resume"],  # "Resume" after no write at all does write new data
    ],
)
def test_import_anndata(adata, ingest_modes):

    for ingest_mode in ingest_modes:

        # Set up anndata input path and tiledb-group output path
        tempdir = tempfile.TemporaryDirectory()
        output_path = tempdir.name

        orig = adata

        # Ingest
        soma = tiledbsoma.SOMA(output_path)
        tiledbsoma.io.from_anndata(soma, orig, ingest_mode=ingest_mode)

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
            assert (
                A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AssayMatrix"
            )
        # Convenience accessors
        assert soma.X["data"].shape() == soma.X.data.shape()

        # Check X/raw (sparse)
        with tiledb.open(os.path.join(output_path, "raw", "X", "data")) as A:
            df = A.df[:]
            assert df.columns.to_list() == ["obs_id", "var_id", "value"]
            # verify sparsity of raw data
            if ingest_mode == "schema_only":
                assert df.shape[0] == 0
            else:
                assert df.shape[0] == orig.raw.X.nnz
            assert (
                A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AssayMatrix"
            )

        # Check obs
        with tiledb.open(os.path.join(output_path, "obs")) as A:
            df = A.df[:]
            assert df.columns.to_list() == orig.obs_keys()
            assert (
                A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
                == "AnnotationDataFrame"
            )
        if ingest_mode == "schema_only":
            assert sorted(soma.obs.ids()) == []
        else:
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
        if ingest_mode == "schema_only":
            assert sorted(soma.var.ids()) == []
        else:
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
                if ingest_mode == "schema_only":
                    assert df.shape[0] == 0
                    assert soma.obsm[key].shape() == (0, orig.obsm[key].shape[1])
                else:
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
                if ingest_mode == "schema_only":
                    assert df.shape[0] == 0
                    assert soma.varm[key].shape() == (0, orig.varm[key].shape[1])
                else:
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
                if ingest_mode == "schema_only":
                    assert df.shape[0] == 0
                else:
                    assert df.shape[0] == orig.obsp[key].nnz
                # https://github.com/single-cell-data/TileDB-SingleCell/issues/125
                # At present (without that PR's suggested enhancement) the best we
                # can get is the NNZ x attrs shape -- note that there are two
                # dims and one attr so the shape is nnz x 1.
                shape = soma.obsp[key].df().shape
                if ingest_mode == "schema_only":
                    assert shape[0] == 0
                else:
                    assert shape[0] == orig.obsp[key].nnz
                assert shape[1] == 1
                assert (
                    A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
                    == "AssayMatrix"
                )
        # Convenience accessors:
        for key in soma.obsp.keys():
            assert getattr(soma.obsp, key).shape() == soma.obsp[key].shape()

    tempdir.cleanup()


def _get_fragment_count(array_uri):
    return len(tiledb.fragment.FragmentInfoList(array_uri=array_uri))


def test_resume_mode(adata):
    """
    Makes sure resume-mode ingest after successful ingest of the same input data does not write
    anything new
    """

    tempdir1 = tempfile.TemporaryDirectory()
    tempdir2 = tempfile.TemporaryDirectory()
    output_path1 = tempdir1.name
    output_path2 = tempdir2.name

    soma1 = tiledbsoma.SOMA(output_path1)
    soma2 = tiledbsoma.SOMA(output_path2)
    tiledbsoma.io.from_anndata(soma1, adata, ingest_mode="write")
    tiledbsoma.io.from_anndata(soma2, adata, ingest_mode="write")
    tiledbsoma.io.from_anndata(soma2, adata, ingest_mode="resume")

    assert _get_fragment_count(soma1.obs.uri) == _get_fragment_count(soma2.obs.uri)
    assert _get_fragment_count(soma1.var.uri) == _get_fragment_count(soma2.var.uri)
    assert _get_fragment_count(soma1.X["data"].uri) == _get_fragment_count(
        soma2.X["data"].uri
    )

    assert _get_fragment_count(soma1.raw.var.uri) == _get_fragment_count(
        soma2.raw.var.uri
    )
    assert _get_fragment_count(soma1.raw.X["data"].uri) == _get_fragment_count(
        soma2.raw.X["data"].uri
    )

    assert _get_fragment_count(soma1.obsm["X_pca"].uri) == _get_fragment_count(
        soma2.obsm["X_pca"].uri
    )
    assert _get_fragment_count(soma1.obsm["X_tsne"].uri) == _get_fragment_count(
        soma2.obsm["X_tsne"].uri
    )

    assert _get_fragment_count(soma1.obsp["distances"].uri) == _get_fragment_count(
        soma2.obsp["distances"].uri
    )

    assert _get_fragment_count(soma1.varm["PCs"].uri) == _get_fragment_count(
        soma2.varm["PCs"].uri
    )

    assert _get_fragment_count(
        soma1.uns["neighbors"]["params"]["method"].uri
    ) == _get_fragment_count(soma2.uns["neighbors"]["params"]["method"].uri)

    tempdir1.cleanup()
    tempdir2.cleanup()


def test_export_anndata(adata):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    orig = adata

    # Ingest
    soma = tiledbsoma.SOMA(output_path)
    tiledbsoma.io.from_anndata(soma, orig)

    readback = tiledbsoma.io.to_anndata(soma)

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


@pytest.mark.parametrize("df_capacity", [1000, 10000, 100000])
def test_df_capacity(adata, df_capacity):

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    # Ingest
    soma_options = tiledbsoma.SOMAOptions(df_capacity=df_capacity)
    soma = tiledbsoma.SOMA(output_path, soma_options=soma_options)
    tiledbsoma.io.from_anndata(soma, adata)

    with soma.obs._open() as D:
        assert D.schema.capacity == df_capacity
    with soma.var._open() as D:
        assert D.schema.capacity == df_capacity

    tempdir.cleanup()
