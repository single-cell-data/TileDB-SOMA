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
    exp = tiledbsoma.SOMAExperiment(output_path)
    tiledbsoma.io.from_anndata(exp, orig, "mRNA")

    # Structure:
    # pbmc-small SOMAExperiment:
    #   obs SOMADataFrame (80,)
    #   ms Collection:
    #     mRNA SOMAMeasurement:
    #       X Collection:
    #         data SOMASparseNdArray (80, 20)
    #       obsp Collection:
    #         distances SOMASparseNdArray (80, 80)
    #       var SOMADataFrame (20,)
    #       obsm Collection:
    #         X_tsne SOMADenseNdArray (80, 2)
    #         X_pca SOMADenseNdArray (80, 19)
    #       varm Collection:
    #         PCs SOMADenseNdArray (20, 19)
    #     raw SOMAMeasurement:
    #       var SOMADataFrame (230,)
    #       X Collection:
    #         data SOMASparseNdArray (80, 230)

    with tiledb.Group(output_path) as G:
        assert G.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "SOMAExperiment"

    # Check obs
    df = exp.obs.read_as_pandas_all()
    assert sorted(df.columns.to_list()) == sorted(
        orig.obs_keys() + ["soma_rowid", "soma_joinid", "obs_id"]
    )
    assert (
        exp.obs.metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
        == "SOMADataFrame"
    )
    assert sorted(df["obs_id"]) == sorted(list(orig.obs_names))
    # Convenience accessor
    assert sorted(exp.obs.keys()) == sorted(
        list(orig.obs.keys()) + ["soma_rowid", "soma_joinid", "obs_id"]
    )

    # Check X/data (dense)
    #    with tiledb.open(os.path.join(output_path, "X", "data")) as A:
    #        df = A[:]
    #        keys = list(df.keys())
    #        assert keys == ["value", "obs_id", "var_id"]
    #        assert A.ndim == 2
    #        assert A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AssayMatrix"
    # Convenience accessors
    #    assert exp.X["data"].shape() == exp.X.data.shape()

    #    # Check X/raw (sparse)
    #    with tiledb.open(os.path.join(output_path, "raw", "X", "data")) as A:
    #        df = A.df[:]
    #        assert df.columns.to_list() == ["obs_id", "var_id", "value"]
    #        # verify sparsity of raw data
    #        assert df.shape[0] == orig.raw.X.nnz
    #        assert A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AssayMatrix"

    # TODO: PORT FROM V0 TO V1

    #    # Check var
    #    with tiledb.open(os.path.join(output_path, "var")) as A:
    #        df = A.df[:]
    #        assert df.columns.to_list() == orig.var_keys()
    #        assert (
    #            A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AnnotationDataFrame"
    #        )
    #    assert sorted(exp.var.ids()) == sorted(list(orig.var_names))
    #    # Convenience accessors
    #    assert exp.var_keys() == exp.var_names
    #    assert exp.var_names == exp.var.ids()
    #    assert exp.n_var == len(exp.var.ids())
    #
    #    # Check some annotation matrices
    #    # Note: pbmc-small doesn't have varp.
    #    assert sorted(exp.obsm.keys()) == sorted(orig.obsm.keys())
    #    for key in orig.obsm_keys():
    #        with tiledb.open(os.path.join(output_path, "obsm", key)) as A:
    #            df = A.df[:]
    #            assert df.shape[0] == orig.obsm[key].shape[0]
    #            assert exp.obsm[key].shape() == orig.obsm[key].shape
    #            assert (
    #                A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
    #                == "AnnotationMatrix"
    #            )
    #    # Convenience accessors: exp.obsm.X_pca <-> exp.obsm['X_pca']
    #    for key in exp.obsm.keys():
    #        assert getattr(exp.obsm, key).shape() == exp.obsm[key].shape()
    #
    #    assert sorted(exp.varm.keys()) == sorted(orig.varm.keys())
    #    for key in orig.varm_keys():
    #        with tiledb.open(os.path.join(output_path, "varm", key)) as A:
    #            df = A.df[:]
    #            assert df.shape[0] == orig.varm[key].shape[0]
    #            assert exp.varm[key].shape() == orig.varm[key].shape
    #            assert (
    #                A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
    #                == "AnnotationMatrix"
    #            )
    #    # Convenience accessors:
    #    for key in exp.varm.keys():
    #        assert getattr(exp.varm, key).shape() == exp.varm[key].shape()
    #
    #    assert sorted(exp.obsp.keys()) == sorted(orig.obsp.keys())
    #    for key in list(orig.obsp.keys()):
    #        with tiledb.open(os.path.join(output_path, "obsp", key)) as A:
    #            df = A.df[:]
    #            assert df.columns.to_list() == ["obs_id_i", "obs_id_j", "value"]
    #            assert df.shape[0] == orig.obsp[key].nnz
    #            # https://github.com/single-cell-data/TileDB-SOMA/issues/125
    #            # At present (without that PR's suggested enhancement) the best we
    #            # can get is the NNZ x attrs shape -- note that there are two
    #            # dims and one attr so the shape is nnz x 1.
    #            shape = exp.obsp[key].df().shape
    #            assert shape[0] == orig.obsp[key].nnz
    #            assert shape[1] == 1
    #            assert A.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "AssayMatrix"
    #    # Convenience accessors:
    #    for key in exp.obsp.keys():
    #        assert getattr(exp.obsp, key).shape() == exp.obsp[key].shape()

    tempdir.cleanup()


# def test_export_anndata(adata):
#
#    # Set up anndata input path and tiledb-group output path
#    tempdir = tempfile.TemporaryDirectory()
#    output_path = tempdir.name
#
#    orig = adata
#
#    # Ingest
#    exp = tiledbsoma.SOMAExperiment(output_path)
#    tiledbsoma.io.from_anndata(exp, orig)
#
#    readback = tiledbsoma.io.to_anndata(exp)
#
#    assert readback.obs.shape == orig.obs.shape
#    assert readback.var.shape == orig.var.shape
#    assert readback.X.shape == orig.X.shape
#
#    for key in orig.obsm.keys():
#        assert readback.obsm[key].shape == orig.obsm[key].shape
#    for key in orig.varm.keys():
#        assert readback.varm[key].shape == orig.varm[key].shape
#    for key in orig.obsp.keys():
#        assert readback.obsp[key].shape == orig.obsp[key].shape
#    for key in orig.varp.keys():
#        assert readback.varp[key].shape == orig.varp[key].shape
