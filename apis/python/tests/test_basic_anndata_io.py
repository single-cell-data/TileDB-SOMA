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
        # Standard ingest: normal use-case:
        ["write"],
        # Schema only:
        ["schema_only"],
        # Schema only, then populate:
        ["schema_only", "resume"],
        # User writes data, then a subsequent write creates nothing new:
        ["write", "resume"],
        # "Resume" after no write at all does write new data:
        ["resume"],
    ],
)
def test_import_anndata(adata, ingest_modes):

    have_ingested = False

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    orig = adata

    for ingest_mode in ingest_modes:

        exp = tiledbsoma.Experiment(output_path)
        tiledbsoma.io.from_anndata(exp, orig, "RNA", ingest_mode=ingest_mode)
        if ingest_mode != "schema_only":
            have_ingested = True

        with tiledb.Group(output_path) as G:
            assert (
                G.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
                == "SOMAExperiment"
            )

        # Check obs
        obs = exp.obs.read_as_pandas_all()
        assert sorted(obs.columns.to_list()) == sorted(
            orig.obs_keys() + ["soma_joinid", "obs_id"]
        )
        assert (
            exp.obs.metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
            == "SOMADataFrame"
        )
        if have_ingested:
            assert sorted(obs["obs_id"]) == sorted(list(orig.obs_names))
        else:
            assert sorted(obs["obs_id"]) == []
        # Convenience accessor
        assert sorted(exp.obs.keys()) == sorted(
            list(orig.obs.keys()) + ["soma_joinid", "obs_id"]
        )

        # Check var
        var = exp.ms["RNA"].var.read_as_pandas_all()
        assert sorted(var.columns.to_list()) == sorted(
            orig.var_keys() + ["soma_joinid", "var_id"]
        )
        assert (
            exp.ms["RNA"].var.metadata.get(
                tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY
            )
            == "SOMADataFrame"
        )
        if have_ingested:
            assert sorted(var["var_id"]) == sorted(list(orig.var_names))
        else:
            assert sorted(var["var_id"]) == []
        # Convenience accessor
        assert sorted(exp.ms["RNA"].var.keys()) == sorted(
            list(orig.var.keys()) + ["soma_joinid", "var_id"]
        )

        # Check X/data (dense)
        if have_ingested:
            X = exp.ms["RNA"].X["data"].read_tensor(coords=(slice(None), slice(None)))
            assert X.shape == orig.X.shape
            assert (
                exp.ms["RNA"]
                .X["data"]
                .metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
                == "SOMADenseNDArray"
            )
        else:
            with pytest.raises(ValueError):
                X = (
                    exp.ms["RNA"]
                    .X["data"]
                    .read_tensor(coords=(slice(None), slice(None)))
                )

        # Check raw/X/data (sparse)
        if have_ingested:
            X = next(
                exp.ms["raw"]
                .X["data"]
                .read_sparse_tensor(coords=(slice(None), slice(None)), format="coo")
            )
            assert X.shape == orig.raw.X.shape
            assert (
                exp.ms["raw"]
                .X["data"]
                .metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
                == "SOMASparseNDArray"
            )
        else:
            with pytest.raises(StopIteration):
                X = next(
                    exp.ms["raw"]
                    .X["data"]
                    .read_sparse_tensor(coords=(slice(None), slice(None)), format="coo")
                )

        obsm = exp.ms["RNA"].obsm
        assert sorted(obsm.keys()) == sorted(orig.obsm.keys())
        for key in list(orig.obsm.keys()):
            if have_ingested:
                matrix = obsm[key].read_tensor(coords=(slice(None), slice(None)))
                assert matrix.shape == orig.obsm[key].shape
                assert (
                    obsm[key].metadata.get(
                        tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY
                    )
                    == "SOMADenseNDArray"
                )
            else:
                with pytest.raises(ValueError):
                    matrix = obsm[key].read_tensor(coords=(slice(None), slice(None)))

        varm = exp.ms["RNA"].varm
        assert sorted(varm.keys()) == sorted(orig.varm.keys())
        for key in list(orig.varm.keys()):
            if have_ingested:
                matrix = varm[key].read_tensor(coords=(slice(None), slice(None)))
                assert matrix.shape == orig.varm[key].shape
                assert (
                    varm[key].metadata.get(
                        tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY
                    )
                    == "SOMADenseNDArray"
                )
            else:
                with pytest.raises(ValueError):
                    matrix = varm[key].read_tensor(coords=(slice(None), slice(None)))

        obsp = exp.ms["RNA"].obsp
        assert sorted(obsp.keys()) == sorted(orig.obsp.keys())
        for key in list(orig.obsp.keys()):
            if have_ingested:
                matrix = next(
                    obsp[key].read_sparse_tensor(
                        coords=(slice(None), slice(None)), format="coo"
                    )
                )
                assert matrix.shape == orig.obsp[key].shape
                assert (
                    obsp[key].metadata.get(
                        tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY
                    )
                    == "SOMASparseNDArray"
                )
            else:
                with pytest.raises(StopIteration):
                    matrix = next(
                        obsp[key].read_sparse_tensor(
                            coords=(slice(None), slice(None)), format="coo"
                        )
                    )

        # pbmc-small has no varp

        tempdir.cleanup()


def _get_fragment_count(array_uri):
    return len(tiledb.fragment.FragmentInfoList(array_uri=array_uri))


@pytest.mark.parametrize(
    "resume_mode_h5ad_file",
    [
        HERE.parent / "anndata/pbmc-small-x-dense.h5ad",
        HERE.parent / "anndata/pbmc-small-x-csr.h5ad",
        HERE.parent / "anndata/pbmc-small-x-csc.h5ad",
    ],
)
def test_resume_mode(adata, resume_mode_h5ad_file):
    """
    Makes sure resume-mode ingest after successful ingest of the same input data does not write
    anything new
    """

    tempdir1 = tempfile.TemporaryDirectory()
    output_path1 = tempdir1.name
    exp1 = tiledbsoma.Experiment(output_path1)
    tiledbsoma.io.from_h5ad(exp1, resume_mode_h5ad_file, "RNA", ingest_mode="write")

    tempdir2 = tempfile.TemporaryDirectory()
    output_path2 = tempdir2.name
    exp2 = tiledbsoma.Experiment(output_path2)
    tiledbsoma.io.from_h5ad(exp2, resume_mode_h5ad_file, "RNA", ingest_mode="write")
    tiledbsoma.io.from_h5ad(exp2, resume_mode_h5ad_file, "RNA", ingest_mode="resume")

    assert _get_fragment_count(exp1.obs.uri) == _get_fragment_count(exp2.obs.uri)
    assert _get_fragment_count(exp1.ms["RNA"].var.uri) == _get_fragment_count(
        exp2.ms["RNA"].var.uri
    )
    assert _get_fragment_count(exp1.ms["RNA"].X["data"].uri) == _get_fragment_count(
        exp2.ms["RNA"].X["data"].uri
    )

    meas1 = exp1.ms["RNA"]
    meas2 = exp2.ms["RNA"]

    if "obsm" in meas1:
        for key in meas1.obsm.keys():
            assert _get_fragment_count(meas1.obsm[key].uri) == _get_fragment_count(
                meas2.obsm[key].uri
            )
    if "varm" in meas1:
        for key in meas1.varm.keys():
            assert _get_fragment_count(meas1.obsm[key].uri) == _get_fragment_count(
                meas2.obsm[key].uri
            )
    if "obsp" in meas1:
        for key in meas1.obsp.keys():
            assert _get_fragment_count(meas1.obsp[key].uri) == _get_fragment_count(
                meas2.obsp[key].uri
            )
    if "varp" in meas1:
        for key in meas1.varp.keys():
            assert _get_fragment_count(meas1.varm[key].uri) == _get_fragment_count(
                meas2.varm[key].uri
            )

    tempdir1.cleanup()
    tempdir2.cleanup()


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
