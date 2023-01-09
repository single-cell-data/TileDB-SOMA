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
    exp = tiledbsoma.Experiment(output_path)
    tiledbsoma.io.from_anndata(exp, orig, "RNA")

    # Structure:
    # pbmc-small Experiment:
    #   obs DataFrame (80,)
    #   ms Collection:
    #     RNA Measurement:
    #       X Collection:
    #         data SparseNDArray (80, 20)
    #       obsp Collection:
    #         distances SparseNDArray (80, 80)
    #       var DataFrame (20,)
    #       obsm Collection:
    #         X_tsne DenseNDArray (80, 2)
    #         X_pca DenseNDArray (80, 19)
    #       varm Collection:
    #         PCs DenseNDArray (20, 19)
    #     raw Measurement:
    #       var DataFrame (230,)
    #       X Collection:
    #         data SparseNDArray (80, 230)

    with tiledb.Group(output_path) as G:
        assert G.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "SOMAExperiment"

    # Check obs
    obs = exp.obs.read().concat().to_pandas()
    assert sorted(obs.columns.to_list()) == sorted(
        orig.obs_keys() + ["soma_joinid", "obs_id"]
    )
    assert (
        exp.obs.metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
        == "SOMADataFrame"
    )
    assert sorted(obs["obs_id"]) == sorted(list(orig.obs_names))
    # Convenience accessor
    assert sorted(exp.obs.keys()) == sorted(
        list(orig.obs.keys()) + ["soma_joinid", "obs_id"]
    )

    # Check var
    var = exp.ms["RNA"].var.read().concat().to_pandas()
    assert sorted(var.columns.to_list()) == sorted(
        orig.var_keys() + ["soma_joinid", "var_id"]
    )
    assert (
        exp.ms["RNA"].var.metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
        == "SOMADataFrame"
    )
    assert sorted(var["var_id"]) == sorted(list(orig.var_names))
    # Convenience accessor
    assert sorted(exp.ms["RNA"].var.keys()) == sorted(
        list(orig.var.keys()) + ["soma_joinid", "var_id"]
    )

    # Check X/data (dense)
    X = exp.ms["RNA"].X["data"].read(coords=(slice(None), slice(None)))
    assert X.shape == orig.X.shape
    assert (
        exp.ms["RNA"]
        .X["data"]
        .metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
        == "SOMADenseNDArray"
    )

    # Check raw/X/data (sparse)
    X = next(exp.ms["raw"].X["data"].read(coords=(slice(None), slice(None))).coos())
    assert X.shape == orig.raw.X.shape
    assert (
        exp.ms["raw"]
        .X["data"]
        .metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
        == "SOMASparseNDArray"
    )

    # Check some annotation matrices
    # Note: pbmc-small doesn't have varp.

    obsm = exp.ms["RNA"].obsm
    assert sorted(obsm.keys()) == sorted(orig.obsm.keys())
    for key in list(orig.obsm.keys()):
        matrix = obsm[key].read(coords=(slice(None), slice(None)))
        assert matrix.shape == orig.obsm[key].shape
        assert (
            obsm[key].metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
            == "SOMADenseNDArray"
        )

    varm = exp.ms["RNA"].varm
    assert sorted(varm.keys()) == sorted(orig.varm.keys())
    for key in list(orig.varm.keys()):
        matrix = varm[key].read(coords=(slice(None), slice(None)))
        assert matrix.shape == orig.varm[key].shape
        assert (
            varm[key].metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
            == "SOMADenseNDArray"
        )

    obsp = exp.ms["RNA"].obsp
    assert sorted(obsp.keys()) == sorted(orig.obsp.keys())
    for key in list(orig.obsp.keys()):
        matrix = next(obsp[key].read(coords=(slice(None), slice(None))).coos())
        assert matrix.shape == orig.obsp[key].shape
        assert (
            obsp[key].metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
            == "SOMASparseNDArray"
        )

    tempdir.cleanup()


def test_export_anndata(adata):
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    orig = adata

    exp = tiledbsoma.Experiment(output_path)
    tiledbsoma.io.from_anndata(exp, orig, measurement_name="RNA")

    readback = tiledbsoma.io.to_anndata(exp, measurement_name="RNA")

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
