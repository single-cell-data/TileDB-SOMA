"""
Miscellaneous tests
"""

from __future__ import annotations

import anndata as ad
import pyarrow as pa

import tiledbsoma as soma
import tiledbsoma.io


def test_open(carrara_group_path: str, carrara_context: soma.SOMAContext) -> None:
    """Test that open works on direct paths to child objects."""

    soma.Collection.create(carrara_group_path, context=carrara_context).close()
    with soma.open(carrara_group_path, mode="w", context=carrara_context) as C:
        C.add_new_sparse_ndarray("snda1", type=pa.int8(), shape=(10, 11))
        C.add_new_collection("c1")

    with soma.open(carrara_group_path, mode="w", context=carrara_context) as C:
        C["c1"].add_new_collection("c1.1")

    with soma.open(carrara_group_path, mode="w", context=carrara_context) as C:
        C["c1"]["c1.1"].add_new_dense_ndarray("dnda1", type=pa.int64(), shape=(100, 10))

    assert soma.open(carrara_group_path).soma_type == "SOMACollection"
    assert soma.open(f"{carrara_group_path}/snda1").soma_type == "SOMASparseNDArray"
    assert soma.open(f"{carrara_group_path}/c1").soma_type == "SOMACollection"
    assert soma.open(f"{carrara_group_path}/c1/c1.1").soma_type == "SOMACollection"
    assert soma.open(f"{carrara_group_path}/c1/c1.1/dnda1").soma_type == "SOMADenseNDArray"


def test_experiment_axis_query(
    small_pbmc: ad.AnnData, carrara_group_path: str, carrara_context: soma.SOMAContext
) -> None:
    """Very basic test that we can run an ExperimentAxisQuery against an Experiment."""
    soma.io.from_anndata(carrara_group_path, small_pbmc, measurement_name="RNA", context=carrara_context)

    with soma.open(carrara_group_path, context=carrara_context) as exp:
        obs_filter = "n_genes>700 and is_b_cell == True"
        with exp.axis_query(obs_query=soma.AxisQuery(value_filter=obs_filter), measurement_name="RNA") as query:
            assert query.n_obs == len(small_pbmc.obs.query(obs_filter))
            assert query.n_vars == len(small_pbmc.var)
            assert (
                (
                    small_pbmc.obs.query(obs_filter)
                    == query.obs().concat().to_pandas().drop(columns=["soma_joinid"]).set_index("obs_id")
                )
                .all()
                .all()
            )

            exported_adata = query.to_anndata(X_name="data")
            assert exported_adata.n_obs == query.n_obs
            assert exported_adata.n_vars == query.n_vars
