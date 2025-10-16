"""
Miscellaneous tests
"""

from __future__ import annotations

import pyarrow as pa

import tiledbsoma as soma


def test_open(carrara_group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    """Test that open works on direct paths to child objects."""

    soma.Collection.create(carrara_group_path, context=carrara_context).close()
    with soma.open(carrara_group_path, mode="w") as C:
        C.add_new_sparse_ndarray("snda1", type=pa.int8(), shape=(10, 11))
        C.add_new_collection("c1")

    with soma.open(carrara_group_path, mode="w") as C:
        C["c1"].add_new_collection("c1.1")

    with soma.open(carrara_group_path, mode="w") as C:
        C["c1"]["c1.1"].add_new_dense_ndarray("dnda1", type=pa.int64(), shape=(100, 10))

    assert soma.open(carrara_group_path).soma_type == "SOMACollection"
    assert soma.open(f"{carrara_group_path}/snda1").soma_type == "SOMASparseNDArray"
    assert soma.open(f"{carrara_group_path}/c1").soma_type == "SOMACollection"
    assert soma.open(f"{carrara_group_path}/c1/c1.1").soma_type == "SOMACollection"
    assert soma.open(f"{carrara_group_path}/c1/c1.1/dnda1").soma_type == "SOMADenseNDArray"
