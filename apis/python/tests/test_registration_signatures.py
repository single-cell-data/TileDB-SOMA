import tempfile
from pathlib import Path

import anndata as ad
import pytest

import tiledbsoma.io
import tiledbsoma.io.registration.signatures as signatures

HERE = Path(__file__).parent


@pytest.fixture
def canned_h5ad_file(request):
    input_path = HERE.parent / "testdata/pbmc-small.h5ad"
    return input_path


@pytest.fixture
def canned_anndata(canned_h5ad_file):
    return ad.read_h5ad(canned_h5ad_file)


def test_signature_serdes(canned_h5ad_file, canned_anndata):
    sig = signatures.Signature.fromH5AD(canned_h5ad_file.as_posix())
    text1 = sig.toJSON()
    assert "obs_schema" in text1
    assert "var_schema" in text1
    assert sig == signatures.Signature.fromJSON(text1)

    sig = signatures.Signature.fromAnnData(canned_anndata)
    text2 = sig.toJSON()
    assert sig == signatures.Signature.fromJSON(text2)

    assert text1 == text2

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    uri = tiledbsoma.io.from_anndata(output_path, canned_anndata, "RNA")
    sig = signatures.Signature.fromSOMAExperiment(uri)
    text3 = sig.toJSON()
    assert sig == signatures.Signature.fromJSON(text3)

    assert text1 == text3


def test_compatible(canned_anndata):

    # Check that zero inputs result in zero incompatibility
    ok, msg = signatures.Signature.compatible({})
    assert ok
    assert msg is None

    sig1 = signatures.Signature.fromAnnData(canned_anndata)

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    uri = tiledbsoma.io.from_anndata(output_path, canned_anndata, "RNA")
    sig2 = signatures.Signature.fromSOMAExperiment(uri)

    # Check that single inputs result in zero incompatibility
    ok, msg = signatures.Signature.compatible({"anndata": sig1})
    assert ok
    assert msg is None

    ok, msg = signatures.Signature.compatible({"experiment": sig2})
    assert ok
    assert msg is None

    # Check that AnnData/H5AD is compatible with itself; likewise with SOMA Experiment
    ok, msg = signatures.Signature.compatible({"anndata": sig1, "same": sig1})
    assert ok
    assert msg is None

    ok, msg = signatures.Signature.compatible({"experiment": sig2, "same": sig2})
    assert ok
    assert msg is None

    # Check compatibility of identical AnnData / SOMA experiment.
    ok, msg = signatures.Signature.compatible({"anndata": sig1, "experiment": sig2})
    assert ok
    assert msg is None

    # Check incompatibility of modified AnnData
    adata3 = canned_anndata
    del adata3.obs["groups"]
    sig3 = signatures.Signature.fromAnnData(adata3)
    ok, msg = signatures.Signature.compatible({"orig": sig1, "anndata3": sig3})
    assert not ok
    assert msg is not None
