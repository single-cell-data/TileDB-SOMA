"""
This module exists to avoid what would otherwise be cyclic-module-import issues within
Collection.
"""

from typing import Mapping, Optional, Type, Union

import tiledb

from .collection import Collection, CollectionBase
from .dataframe import DataFrame
from .dense_nd_array import DenseNDArray
from .exception import SOMAError
from .experiment import Experiment
from .measurement import Measurement
from .options import SOMATileDBContext
from .sparse_nd_array import SparseNDArray
from .tiledb_array import TileDBArray
from .tiledb_object import TileDBObject
from .util import (
    SOMA_ENCODING_VERSION,
    SOMA_ENCODING_VERSION_METADATA_KEY,
    SOMA_OBJECT_TYPE_METADATA_KEY,
)

SPEC_NAME_TO_CLASS: Mapping[str, Type[TileDBObject]] = {
    # Maps languge-independent-spec names to Python-implementation class
    "SOMAExperiment": Experiment,
    "SOMAMeasurement": Measurement,
    "SOMACollection": Collection,
    "SOMADataFrame": DataFrame,
    "SOMADenseNDArray": DenseNDArray,
    "SOMADenseNdArray": DenseNDArray,
    "SOMASparseNDArray": SparseNDArray,
    "SOMASparseNdArray": SparseNDArray,
}


def _construct_member(
    member_uri: str,
    context: SOMATileDBContext,
    object_type: Type[Union[tiledb.Array, tiledb.Group]],
) -> Optional[TileDBObject]:
    """
    Given a name/uri from a Collection, create a SOMA object matching the type
    of the underlying object. In other words, if the name/uri points to an DataFrame,
    instantiate an DataFrame pointing at the underlying array.

    Returns None if the URI does not point at a TileDB object, or if the TileDB
    object is not recognized as a SOMA object.

    Solely for the use of ``Collection``. In fact this would/should be a method of the
    ``Collection`` class, but there are cyclic-module-import issues.  This allows us to
    examine storage metadata and invoke the appropriate per-type constructor when reading
    SOMA groups/arrays from storage.  See also ``_set_object_type_metadata`` and
    ``_get_object_type_metadata`` within ``TileDBObject``.
    """
    # auto-detect class name from metadata
    if object_type is tiledb.Array:
        tdb_open = tiledb.open
    elif object_type is tiledb.Group:
        tdb_open = tiledb.Group
    else:
        return None

    try:
        with tdb_open(member_uri, ctx=context.tiledb_ctx) as o:
            spec_name = o.meta.get(SOMA_OBJECT_TYPE_METADATA_KEY, None)
            encoding_version = o.meta.get(SOMA_ENCODING_VERSION_METADATA_KEY)
    except tiledb.TileDBError:
        return None

    if spec_name is None:
        raise SOMAError("internal error: spec_name was not found")
    if encoding_version is None:
        raise SOMAError("internal error: encoding_version not found")
    if encoding_version != SOMA_ENCODING_VERSION:
        raise ValueError("Unsupported SOMA object encoding version")

    # Now invoke the appropriate per-class constructor.
    cls = SPEC_NAME_TO_CLASS.get(spec_name)
    if cls is None:
        raise SOMAError(f'name "{spec_name}" unrecognized')

    if issubclass(cls, CollectionBase) and object_type is not tiledb.Group:
        raise SOMAError(
            f'internal error: expected "group" object_type; got "{object_type}"'
        )

    if issubclass(cls, TileDBArray) and object_type is not tiledb.Array:
        raise SOMAError(
            f'internal error: expected "array" object_type; got "{object_type}"'
        )

    return cls(uri=member_uri, context=context)
