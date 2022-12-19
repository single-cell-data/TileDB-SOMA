"""
This module exists to avoid what would otherwise be cyclic-module-import issues within
Collection.
"""

from typing import Any, Optional, Union

import tiledb

from .collection import Collection, CollectionBase
from .dataframe import DataFrame
from .dense_nd_array import DenseNDArray
from .exception import SOMAError
from .experiment import Experiment
from .measurement import Measurement
from .sparse_nd_array import SparseNDArray
from .util import SOMA_OBJECT_TYPE_METADATA_KEY, SPEC_NAMES_TO_CLASS_NAMES

ObjectTypes = Union[
    Experiment,
    Measurement,
    Collection,
    DataFrame,
    DenseNDArray,
    SparseNDArray,
]


def _construct_member(
    member_uri: str,
    parent: CollectionBase[Any],
    ctx: Optional[tiledb.Ctx] = None,
    object_type: Optional[str] = None,
) -> Optional[ObjectTypes]:
    """
    Given a name/uri from a Collection, create a SOMA object matching the type
    of the underlying object. In other words, if the name/uri points to an DataFrame,
    instantiate an DataFrame pointing at the underlying array.

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
                spec_name = A.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
        elif object_type == "group":
            with tiledb.Group(member_uri, mode="r", ctx=ctx) as G:
                spec_name = G.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
        else:
            return None

    except tiledb.TileDBError:
        return None

    assert spec_name is not None
    if spec_name not in SPEC_NAMES_TO_CLASS_NAMES:
        raise SOMAError(f'name "{spec_name}" unrecognized')
    class_name = SPEC_NAMES_TO_CLASS_NAMES[spec_name]

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
    elif class_name in ["DenseNDArray", "DenseNdArray"]:
        assert object_type is None or object_type == "array"
        return DenseNDArray(uri=member_uri, parent=parent, ctx=ctx)
    elif class_name in ["SparseNDArray", "SparseNdArray"]:
        assert object_type is None or object_type == "array"
        return SparseNDArray(uri=member_uri, parent=parent, ctx=ctx)
    else:
        raise SOMAError(
            f'internal coding error: class name "{class_name}" unrecognized'
        )
