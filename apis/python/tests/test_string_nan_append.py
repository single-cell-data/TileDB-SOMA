from contextlib import nullcontext

import pytest
import scanpy as sc
import numpy as np
import tiledbsoma
import tiledbsoma.io
import tempfile


@pytest.fixture
def adata_simple(adata):
    adata.obsm = None
    adata.varm = None
    adata.obsp = None
    adata.varp = None
    adata.uns = dict()
    return adata


def write_ads(adat1, adat2, exc=None):
    # Initial ingest
    SOMA_URI = tempfile.mkdtemp(prefix="soma-exp-")
    tiledbsoma.io.from_anndata(
        experiment_uri=SOMA_URI,
        anndata=adat1,
        measurement_name="RNA"
    )

    # Register the second anndata object
    rd = tiledbsoma.io.register_anndatas(
        experiment_uri=SOMA_URI,
        adatas=[adat2],
        measurement_name="RNA",
        obs_field_name="obs_id",
        var_field_name="var_id"
    )

    ctx = pytest.raises(exc) if exc else nullcontext()
    with ctx:
        # Append the second anndata object
        tiledbsoma.io.from_anndata(
            experiment_uri=SOMA_URI,
            anndata=adat2,
            measurement_name="RNA",
            registration_mapping=rd
        )


@pytest.mark.parametrize("rename_obs_index", [False, True])
@pytest.mark.parametrize("dtype", ["float64", "string"])
def test_string_nan_append_small(adata_simple, rename_obs_index, dtype):
    adat1 = adata_simple

    # Add empty column to obs
    adat1.obs['batch_id'] = np.nan
    adat1.obs['batch_id'] = adat1.obs['batch_id'].astype(dtype)

    # Create a copy of the anndata object
    adat2 = adat1.copy()
    if rename_obs_index:
        adat2.obs.index = adat1.obs.index + "-2"

    write_ads(adat1, adat2)


@pytest.mark.parametrize("rename_obs_index", [False, True])
def test_string_nan_append_3k(rename_obs_index):
    adat1 = sc.datasets.pbmc3k()

    # Add empty column to obs
    adat1.obs['batch_id'] = np.nan
    adat1.obs['batch_id'] = adat1.obs['batch_id'].astype('string')

    # Create a copy of the anndata object
    adat2 = adat1.copy()
    if rename_obs_index:
        adat2.obs.index = adat1.obs.index.str.replace("1", "2")

    write_ads(adat1, adat2)


def test_empty_append(adata_simple):
    adat1 = adata_simple
    adat2 = adat1[0:0].copy()
    write_ads(adat1, adat2, exc=NotImplementedError)
