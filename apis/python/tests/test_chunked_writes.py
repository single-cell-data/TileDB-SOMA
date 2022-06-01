import anndata
import tiledb
import tiledbsc
import tiledbsc.io

import pytest
import tempfile
import os
from pathlib import Path
import math

HERE = Path(__file__).parent


def relerr(a, b):
    m = max(math.fabs(a), math.fabs(b))
    if m == 0.0:
        return 0
    return math.fabs(a - b) / m


def relerr_ok(a, b):
    return relerr(a, b) < 1e-6


def test_chunked_writes(tmp_path):
    """
    This ingests datasets having X stored in dense, CSR, and CSC formats, and
    verifies the chunked-ingest logic is the same.
    """

    # We have three extra copies of pbmc-small.h5ad checked into source control.
    # See the README.md file there for how they were prepared -- essentially
    #   ann = anndata.read_h5ad('pbmc-small.h5ad')
    #   ann.X = scipy.sparse.csc_matrix(ann.X)
    #   ann.write_h5ad('pbmc-small-x-csc.h5ad')
    # along with nulling out obsm/uns/etc which are not used in this test (to make the files a
    # little smaller).
    ann_dense_path = HERE.parent / "anndata/pbmc-small-x-dense.h5ad"
    ann_csr_path = HERE.parent / "anndata/pbmc-small-x-csr.h5ad"
    ann_csc_path = HERE.parent / "anndata/pbmc-small-x-csc.h5ad"

    soma_unchunked_path = (tmp_path / "soma_unchunked").as_posix()
    soma_dense_path = (tmp_path / "soma_dense").as_posix()
    soma_csr_path = (tmp_path / "soma_csr").as_posix()
    soma_csc_path = (tmp_path / "soma_csc").as_posix()

    # This chunk-size value is specifically chosen to ensure that the curated datafiles on this test
    # will get their X ingested in multiple chunks -- this is important in order to test the
    # chunked-ingest logic for all three variants: dense, CSR, and CSC.
    unchunked_sopt = tiledbsc.SOMAOptions(write_X_chunked=False)
    chunked_sopt = tiledbsc.SOMAOptions(goal_chunk_nnz=250)

    # Do the ingests from anndata .h5ad files to TileDB SOMAs.
    soma_unchunked = tiledbsc.SOMA(soma_unchunked_path, soma_options=unchunked_sopt)
    soma_dense = tiledbsc.SOMA(soma_dense_path, soma_options=chunked_sopt)
    soma_csr = tiledbsc.SOMA(soma_csr_path, soma_options=chunked_sopt)
    soma_csc = tiledbsc.SOMA(soma_csc_path, soma_options=chunked_sopt)

    tiledbsc.io.from_h5ad(soma_unchunked, ann_dense_path)
    tiledbsc.io.from_h5ad(soma_dense, ann_dense_path)
    tiledbsc.io.from_h5ad(soma_csr, ann_csr_path)
    tiledbsc.io.from_h5ad(soma_csc, ann_csc_path)

    # Read the X arrays into memory as pandas dataframe objects.
    xdf_unchunked = soma_unchunked.X["data"].df()
    xdf_dense = soma_dense.X["data"].df()
    xdf_csr = soma_csr.X["data"].df()
    xdf_csc = soma_csc.X["data"].df()

    # IJV/COO triples are correct, but are not necessarily in the same order for CSC.
    # Sort them to facilitate integer-indexed value comparisons below.
    xdf_unchunked.sort_index(inplace=True)
    xdf_dense.sort_index(inplace=True)
    xdf_csr.sort_index(inplace=True)
    xdf_csc.sort_index(inplace=True)

    # Verify that the X data that made it through the chunking algorithms is the same.

    # A first check is the shape of the dataframes. Note that once densified
    # these are (2638, 1838) but as dataframes they're (2638*1838, 3): there
    # is an 'obs_id' column, a 'var_id' column, and a 'value' column.
    assert xdf_unchunked.shape == xdf_dense.shape
    assert xdf_unchunked.shape == xdf_csr.shape
    assert xdf_unchunked.shape == xdf_csc.shape

    # A second check is that the indices are the same.
    assert list(xdf_unchunked.index) == list(xdf_dense.index)
    assert list(xdf_unchunked.index) == list(xdf_csr.index)
    assert list(xdf_unchunked.index) == list(xdf_csc.index)

    # A third check is that the values are the same.
    # Check X values pointwise for all ('obs_id', 'var_id') pairs.
    for i in range(xdf_unchunked.shape[0]):
        assert relerr_ok(xdf_unchunked.iloc[i].value, xdf_dense.iloc[i].value)
        assert relerr_ok(xdf_unchunked.iloc[i].value, xdf_csr.iloc[i].value)
        assert relerr_ok(xdf_unchunked.iloc[i].value, xdf_csc.iloc[i].value)
