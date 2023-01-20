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
from .options import SOMATileDBContext
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
    context: Optional[SOMATileDBContext] = None,
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

    context = context or SOMATileDBContext()

    # Get the class name from TileDB storage. At the TileDB level there are just "arrays" and
    # "groups", with separate metadata-getters.
    if object_type is None:
        object_type = tiledb.object_type(member_uri, ctx=context.tiledb_ctx)

    # auto-detect class name from metadata
    try:
        if object_type == "array":
            with tiledb.open(member_uri, ctx=context.tiledb_ctx) as A:
                spec_name = A.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
        elif object_type == "group":
            with tiledb.Group(member_uri, mode="r", ctx=context.tiledb_ctx) as G:
                spec_name = G.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
        else:
            return None

    except tiledb.TileDBError:
        return None

    if spec_name is None:
        raise SOMAError("internal error: spec_name was not found")
    if spec_name not in SPEC_NAMES_TO_CLASS_NAMES:
        raise SOMAError(f'name "{spec_name}" unrecognized')
    class_name = SPEC_NAMES_TO_CLASS_NAMES[spec_name]

    # Now invoke the appropriate per-class constructor.
    if class_name == "Experiment":
        _check_object_type(object_type, "group")
        return Experiment(uri=member_uri, parent=parent)
    elif class_name == "Measurement":
        _check_object_type(object_type, "group")
        return Measurement(uri=member_uri, parent=parent)
    elif class_name == "Collection":
        _check_object_type(object_type, "group")
        return Collection(uri=member_uri, parent=parent)
    elif class_name == "DataFrame":
        _check_object_type(object_type, "array")
        return DataFrame(uri=member_uri, parent=parent)
    elif class_name in ["DenseNDArray", "DenseNdArray"]:
        _check_object_type(object_type, "array")
        return DenseNDArray(uri=member_uri, parent=parent)
    elif class_name in ["SparseNDArray", "SparseNdArray"]:
        _check_object_type(object_type, "array")
        return SparseNDArray(uri=member_uri, parent=parent)
    else:
        raise SOMAError(f'internal error: class name "{class_name}" unrecognized')


def _check_object_type(
    actual_object_type: Union[None, str], expected_object_type: str
) -> None:
    """Helper function for `_construct_member`"""
    if actual_object_type is not None and actual_object_type != expected_object_type:
        raise SOMAError(
            f'internal error: expected object_type None or "{expected_object_type}"; got "{actual_object_type}"'
        )
