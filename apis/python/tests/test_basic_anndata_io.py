import inspect
import pathlib
import tempfile
from pathlib import Path
from typing import Optional

import anndata
import numpy as np
import pandas as pd
import pytest
import somacore
import tiledb

import tiledbsoma
import tiledbsoma.io
from tiledbsoma import _constants, _factory

HERE = Path(__file__).parent


@pytest.fixture
def h5ad_file(request):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
    # pbmc-small is faster for automated unit-test / CI runs.
    input_path = HERE.parent / "testdata/pbmc-small.h5ad"
    # input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_extended(request):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
    # This has more component arrays in it
    input_path = HERE.parent / "testdata/pbmc3k_processed.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_with_obsm_holes(request):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
    # This has zeroes in an obsm matrix so nnz is not num_rows * num_cols
    input_path = HERE.parent / "testdata/pbmc3k-with-obsm-zero.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_uns_string_arrays(request):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
    # This has uns["louvain_colors"] with dtype.char == "U".
    # It also has uns["more_colors"] in the form '[[...]]', as often occurs in the wild.
    input_path = HERE.parent / "testdata/pbmc3k.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_categorical_int_nan(request):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
    # This has obs["categ_int_nan"] as a categorical int but with math.nan as a
    # "not-in-the-category" indicator. Such H5AD files do arise in the wild.
    #
    # Reference:
    #   import anndata as ad
    #   import pandas  as pd
    #   import math
    #   adata = adata.read_h5ad("whatever.h5ad")
    #   s = pd.Series(list(range(80)), dtype="category")
    #   s[0] = math.nan
    #   adata.obs["categ_int_nan"] = s
    #   adata.write_h5ad("categorical_int_nan.h5ad")
    input_path = HERE.parent / "testdata/categorical_int_nan.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_X_empty(request):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
    """adata.X is a zero-cell sparse matrix"""
    input_path = HERE.parent / "testdata/x-empty.h5ad"
    return input_path


@pytest.fixture
def h5ad_file_X_none(request):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
    """
    adata.X has Python value None if read in non-backed mode; if read in backed
    mode, adata.X is not present as an attribute of adata.
    """
    input_path = HERE.parent / "testdata/x-none.h5ad"
    return input_path


@pytest.fixture
def adata(h5ad_file):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
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
@pytest.mark.parametrize(
    "X_kind",
    [tiledbsoma.SparseNDArray, tiledbsoma.DenseNDArray],
)
def test_import_anndata(adata, ingest_modes, X_kind):
    print(f"\nENTER {inspect.currentframe().f_code.co_name}\n")
    adata = adata.copy()

    have_ingested = False

    tempdir = tempfile.TemporaryDirectory()
    output_path = tempdir.name

    adata.layers["plus1"] = adata.X + 1
    orig = adata.copy()

    metakey = _constants.SOMA_OBJECT_TYPE_METADATA_KEY  # keystroke-saver
    all2d = (slice(None), slice(None))  # keystroke-saver

    for ingest_mode in ingest_modes:
        uri = tiledbsoma.io.from_anndata(
            output_path,
            orig,
            "RNA",
            ingest_mode=ingest_mode,
            X_kind=X_kind,
        )
        if ingest_mode != "schema_only":
            have_ingested = True

        exp = tiledbsoma.Experiment.open(uri)

        assert exp.metadata[metakey] == "SOMAExperiment"

        # Check obs
        obs = exp.obs.read().concat().to_pandas()
