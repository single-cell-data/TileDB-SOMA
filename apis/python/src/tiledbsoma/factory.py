"""
This module exists to avoid what would otherwise be cyclic-module-import issues within
SOMACollection.
"""

from typing import Optional, Union

import tiledb

from .soma_collection import SOMACollection as SOMACollection
from .soma_dataframe import SOMADataFrame as SOMADataFrame
from .soma_dense_nd_array import SOMADenseNdArray as SOMADenseNdArray
from .soma_experiment import SOMAExperiment as SOMAExperiment
from .soma_indexed_dataframe import SOMAIndexedDataFrame as SOMAIndexedDataFrame
from .soma_measurement import SOMAMeasurement as SOMAMeasurement
from .soma_sparse_nd_array import SOMASparseNdArray as SOMASparseNdArray
from .util import SOMA_OBJECT_TYPE_METADATA_KEY

MemberType = Union[
    SOMAExperiment,
    SOMAMeasurement,
    SOMACollection,
    SOMADataFrame,
    SOMAIndexedDataFrame,
    SOMADenseNdArray,
    SOMASparseNdArray,
]


def _construct_member(
    member_name: str,
    member_uri: str,
    parent: SOMACollection,
    ctx: Optional[tiledb.Ctx] = None,
) -> Optional[MemberType]:
    """
    Given a name/uri from a SOMACollection, create a SOMA object matching the type
    of the underlying object. In other words, if the name/uri points to a SOMADataFrame,
    instantiate a SOMADataFrame pointing at the underlying array.

    Returns None if the URI does not point at a TileDB object.

    Solely for the use of ``SOMACollection``. In fact this would/should be a method of the ``SOMACollection`` class,
    but there are cyclic-module-import issues.  This allows us to examine storage metadata and invoke the appropriate
    per-type constructor when reading SOMA groups/arrays from storage.  See also ``_set_object_type_metadata`` and
    ``_get_object_type_metadata`` within ``TileDBObject``.
    """

    # Get the class name from TileDB storage. At the TileDB level there are just "arrays" and
    # "groups", with separate metadata-getters.
    class_name = None
    object_type = tiledb.object_type(member_uri, ctx=ctx)
    if object_type is None:
        return None
    elif object_type == "array":
        with tiledb.open(member_uri, ctx=ctx) as A:
            class_name = A.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
    elif object_type == "group":
        with tiledb.Group(member_uri, mode="r", ctx=ctx) as G:
            class_name = G.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
    else:
        raise Exception(f"object type {object_type} unrecognized")
    assert class_name is not None

    # Now invoke the appropriate per-class constructor.
    if class_name == "SOMAExperiment":
        return SOMAExperiment(uri=member_uri, name=member_name, parent=parent, ctx=ctx)
    elif class_name == "SOMAMeasurement":
        return SOMAMeasurement(uri=member_uri, name=member_name, parent=parent, ctx=ctx)
    elif class_name == "SOMACollection":
        return SOMACollection(uri=member_uri, name=member_name, parent=parent, ctx=ctx)
    elif class_name == "SOMADataFrame":
        return SOMADataFrame(uri=member_uri, name=member_name, parent=parent, ctx=ctx)
    elif class_name == "SOMADenseNdArray":
        return SOMADenseNdArray(
            uri=member_uri, name=member_name, parent=parent, ctx=ctx
        )
    elif class_name == "SOMASparseNdArray":
        return SOMASparseNdArray(
            uri=member_uri, name=member_name, parent=parent, ctx=ctx
        )
    else:
        raise Exception(f'class name "{class_name}" unrecognized')
