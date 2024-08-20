import tempfile

import pytest

import tiledbsoma.io
import tiledbsoma.io._registration.signatures as signatures

from ._util import assert_adata_equal


def test_signature_serdes(conftest_pbmc_small_h5ad_path, conftest_pbmc_small):
    sig = signatures.Signature.from_h5ad(conftest_pbmc_small_h5ad_path.as_posix())
    text1 = sig.to_json()
    assert "obs_schema" in text1
    assert "var_schema" in text1
    assert sig == signatures.Signature.from_json(text1)

    original = conftest_pbmc_small.copy()
    sig = signatures.Signature.from_anndata(conftest_pbmc_small)
    assert_adata_equal(original, conftest_pbmc_small)

    text2 = sig.to_json()
    assert sig == signatures.Signature.from_json(text2)

    assert text1 == text2

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    uri = tiledbsoma.io.from_anndata(output_path, conftest_pbmc_small, "RNA")
    assert_adata_equal(original, conftest_pbmc_small)

    sig = signatures.Signature.from_soma_experiment(uri)
    text3 = sig.to_json()
    assert sig == signatures.Signature.from_json(text3)

    assert text1 == text3


# Magical conftest.py fixture
def test_compatible(conftest_pbmc_small):
    # Check that zero inputs result in zero incompatibility
    signatures.Signature.check_compatible({})

    original = conftest_pbmc_small.copy()
    sig1 = signatures.Signature.from_anndata(conftest_pbmc_small)
    assert_adata_equal(original, conftest_pbmc_small)

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name
    uri = tiledbsoma.io.from_anndata(output_path, conftest_pbmc_small, "RNA")
    assert_adata_equal(original, conftest_pbmc_small)
    sig2 = signatures.Signature.from_soma_experiment(uri)

    # Check that single inputs result in zero incompatibility
    signatures.Signature.check_compatible({"anndata": sig1})  # no throw
    signatures.Signature.check_compatible({"experiment": sig2})  # no throw

    # Check that AnnData/H5AD is compatible with itself; likewise with SOMA Experiment
    signatures.Signature.check_compatible({"anndata": sig1, "same": sig1})  # no throw
    signatures.Signature.check_compatible(
        {"experiment": sig2, "same": sig2}
    )  # no throw

    # Check compatibility of identical AnnData / SOMA experiment.
    signatures.Signature.check_compatible(
        {"anndata": sig1, "experiment": sig2}
    )  # no throw

    # Check incompatibility of modified AnnData
    adata3 = conftest_pbmc_small
    del adata3.obs["groups"]

    original = adata3.copy()
    sig3 = signatures.Signature.from_anndata(adata3)
    assert_adata_equal(original, adata3)

    with pytest.raises(ValueError):
        signatures.Signature.check_compatible({"orig": sig1, "anndata3": sig3})
