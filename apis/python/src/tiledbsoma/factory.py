"""
This module exists to avoid what would otherwise be cyclic-module-import issues within
Collection.
"""

from typing import Any, Optional, Union

import tiledb

from .soma_collection import Collection, CollectionBase
from .soma_dataframe import DataFrame
from .soma_dense_nd_array import DenseNdArray
from .soma_experiment import Experiment
from .soma_indexed_dataframe import IndexedDataFrame
from .soma_measurement import Measurement
from .soma_sparse_nd_array import SparseNdArray
from .util import SOMA_OBJECT_TYPE_METADATA_KEY

SOMAObjectTypes = Union[
    Experiment,
    Measurement,
    Collection,
    DataFrame,
    IndexedDataFrame,
    DenseNdArray,
    SparseNdArray,
]


def _construct_member(
    member_uri: str,
    parent: CollectionBase[Any],
    ctx: Optional[tiledb.Ctx] = None,
    object_type: Optional[str] = None,
) -> Optional[SOMAObjectTypes]:
    """
    Given a name/uri from a Collection, create a SOMA object matching the type
    of the underlying object. In other words, if the name/uri points to a DataFrame,
    instantiate a DataFrame pointing at the underlying array.

    Returns None if the URI does not point at a TileDB object, or if the TileDB
    object is not recognized as a SOMA object.

    Solely for the use of ``Collection``. In fact this would/should be a method of the ``Collection`` class,
    but there are cyclic-module-import issues.  This allows us to examine storage metadata and invoke the appropriate
    per-type constructor when reading SOMA groups/arrays from storage.  See also ``_set_object_type_metadata`` and
    ``_get_object_type_metadata`` within ``TileDBObject``.
    """

    # Get the class name from TileDB storage. At the TileDB level there are just "arrays" and
    # "groups", with separate metadata-getters.
    if object_type is None:
        object_type = tiledb.object_type(member_uri, ctx=ctx)

    # auto-detect class name from metadata
    try:
        if object_type == "array":
            with tiledb.open(member_uri, ctx=ctx) as A:
                class_name = A.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
        elif object_type == "group":
            with tiledb.Group(member_uri, mode="r", ctx=ctx) as G:
                class_name = G.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
        else:
            return None

    except tiledb.TileDBError:
        return None

    assert class_name is not None

    # Now invoke the appropriate per-class constructor.
    if class_name == "Experiment":
        assert object_type is None or object_type == "group"
        return Experiment(uri=member_uri, parent=parent, ctx=ctx)
    elif class_name == "Measurement":
        assert object_type is None or object_type == "group"
        return Measurement(uri=member_uri, parent=parent, ctx=ctx)
    elif class_name == "Collection":
        assert object_type is None or object_type == "group"
        return Collection(uri=member_uri, parent=parent, ctx=ctx)
    elif class_name == "DataFrame":
        assert object_type is None or object_type == "array"
        return DataFrame(uri=member_uri, parent=parent, ctx=ctx)
    elif class_name == "IndexedDataFrame":
        assert object_type is None or object_type == "array"
        return IndexedDataFrame(uri=member_uri, parent=parent, ctx=ctx)
    elif class_name == "DenseNdArray":
        assert object_type is None or object_type == "array"
        return DenseNdArray(uri=member_uri, parent=parent, ctx=ctx)
    elif class_name == "SparseNdArray":
        assert object_type is None or object_type == "array"
        return SparseNdArray(uri=member_uri, parent=parent, ctx=ctx)
    else:
        # TODO: SOMA needs better exception typing
        raise Exception(f'class name "{class_name}" unrecognized')
