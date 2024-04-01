import scanpy as sc
import numpy as np
import tiledbsoma
import tiledbsoma.io
import tempfile


def test_string_nan_append():
    # Load example anndata
    adat1 = sc.datasets.pbmc3k()

    # Add empty column to obs
    adat1.obs['batch_id'] = np.nan
    adat1.obs['batch_id'] = adat1.obs['batch_id'].astype('string')

    # Create a copy of the anndata object with new obs ids
    adat2 = adat1.copy()

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

    # Append the second anndata object
    tiledbsoma.io.from_anndata(
        experiment_uri=SOMA_URI,
        anndata=adat2,
        measurement_name="RNA",
        registration_mapping=rd
    )
