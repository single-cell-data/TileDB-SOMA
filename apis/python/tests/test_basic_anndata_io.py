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
        [
            "write",
            "resume",
        ],
        # "Resume" after no write at all does write new data
        ["resume"],
    ],
)
def test_import_anndata(adata, ingest_modes):

    have_ingested = False

    # Set up anndata input path and tiledb-group output path
    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    for ingest_mode in ingest_modes:

        orig = adata

        # Ingest
        exp = tiledbsoma.Experiment(output_path)
        tiledbsoma.io.from_anndata(exp, orig, "mRNA", ingest_mode=ingest_mode)

        if ingest_mode != "schema_only":
            have_ingested = True

        with tiledb.Group(output_path) as G:
            assert (
                G.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY]
                == "SOMAExperiment"
            )

        # Check obs
        df = exp.obs.read_as_pandas_all()
        assert sorted(df.columns.to_list()) == sorted(
            orig.obs_keys() + ["soma_joinid", "obs_id"]
        )
        assert (
            exp.obs.metadata.get(tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY)
            == "SOMADataFrame"
        )
        # Convenience accessor
        assert sorted(exp.obs.keys()) == sorted(
            list(orig.obs.keys()) + ["soma_joinid", "obs_id"]
        )
        if have_ingested:
            assert sorted(df["obs_id"]) == sorted(list(orig.obs_names))
        else:
            assert sorted(df["obs_id"]) == []

        # XXX MORE FROM MAIN-OLD PERHAPS AS UNDERDIFF

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

    # XXX TO DO
    #    assert _get_fragment_count(exp1.raw.var.uri) == _get_fragment_count(
    #        exp2.raw.var.uri
    #    )
    #    assert _get_fragment_count(exp1.raw.X["data"].uri) == _get_fragment_count(
    #        exp2.raw.X["data"].uri
    #    )
    #
    #    assert _get_fragment_count(exp1.obsm["X_pca"].uri) == _get_fragment_count(
    #        exp2.obsm["X_pca"].uri
    #    )
    #    assert _get_fragment_count(exp1.obsm["X_tsne"].uri) == _get_fragment_count(
    #        exp2.obsm["X_tsne"].uri
    #    )
    #
    #    assert _get_fragment_count(exp1.obsp["distances"].uri) == _get_fragment_count(
    #        exp2.obsp["distances"].uri
    #    )
    #
    #    assert _get_fragment_count(exp1.varm["PCs"].uri) == _get_fragment_count(
    #        exp2.varm["PCs"].uri
    #    )

    tempdir1.cleanup()
    tempdir2.cleanup()
