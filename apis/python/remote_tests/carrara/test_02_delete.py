"""
Delete tests with Carrara URIs.

Carrara semantics differ from S3/Posix/et al in the following:
- removing a member from a group will deregister (but not delete), the underlying TileDB object. This makes the child object
  in-accessible via tiledb URIs, requiring direct storage system access to delete/cleanup (if desired).
- ...

"""

from __future__ import annotations

import pyarrow as pa
import pytest

import tiledbsoma as soma
import tiledb

from ._util import get_asset_info


@pytest.mark.carrara
def test_collection_delete_member_by_name(carrara_group_path: str, carrara_context: soma.SOMAContext) -> None:
    """
    This tests the behavior when an object is removed from a group.
    1. create parent collection
    2. add members
    3. remove one of the members from the parent
    4. verify that
    """
    soma.Collection.create(carrara_group_path, context=carrara_context).close()
    with soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C:
        C.add_new_sparse_ndarray("array", type=pa.int8(), shape=(10,))
        C.add_new_collection("collection")

    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert set(C) == {"array", "collection"}

    with soma.open(carrara_group_path, mode="d") as C:
        del C["array"]

    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert set(C) == {"collection"}

    assert not soma.SparseNDArray.exists(f"{carrara_group_path}/array", context=carrara_context)


@pytest.mark.carrara
def test_collection_delete_member_object(carrara_group_path: str, carrara_context: soma.SOMAContext) -> None:
    """
    This tests the behavior when an object is deleted - does it correctly clean up parent group?.
    """
    import tiledb.client
    soma.Collection.create(carrara_group_path, context=carrara_context).close()
    with soma.Collection.open(carrara_group_path, mode="w", context=carrara_context) as C:
        C.add_new_sparse_ndarray("array", type=pa.int8(), shape=(10,))
        C.add_new_collection("collection")

    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert set(C) == {"array", "collection"}

        # stash asset info for all paths in use
        array_asset_info = get_asset_info(C["array"].uri)

    # Delete the array object, which is a member of the parent group
    tiledb.client.assets.delete_asset(f"{carrara_group_path}/array", delete_storage=True)

    # verify that the parent group no longer contains the object, and that
    # the other objects still exist.
    with soma.open(carrara_group_path, context=carrara_context) as C:
        assert set(C) == {"collection"}

    # verify that the object does not exist on the storage layer
    assert not soma.SparseNDArray.exists(array_asset_info["uri"], context=carrara_context)
