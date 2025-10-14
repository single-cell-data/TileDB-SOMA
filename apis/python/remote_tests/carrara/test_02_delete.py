"""
Delete tests with Carrara URIs

TODO:
- delete member from collection, with and without member object deletion
- move a member to a new group
etc
"""

from __future__ import annotations

import pyarrow as pa
import pytest

import tiledbsoma as soma
import tiledb


@pytest.mark.xfail(reason="Delete semantics are undefined/surprising")
@pytest.mark.carrara
def test_collection_delete(group_path: str, carrara_context: soma.SOMATileDBContext) -> None:
    soma.Collection.create(group_path, context=carrara_context).close()
    with soma.Collection.open(group_path, mode="w", context=carrara_context) as C:
        C.add_new_sparse_ndarray("array", type=pa.int8(), shape=(10,))
        C.add_new_collection("collection")

    with soma.open(group_path, context=carrara_context) as C:
        assert set(C) == {"array", "collection"}

    with soma.open(group_path, mode="d") as C:
        del C["array"]

    with soma.open(group_path, context=carrara_context) as C:
        assert set(C) == {"collection"}

    tiledb.Array.delete_array(f"{group_path}/array", ctx=carrara_context.tiledb_ctx)
