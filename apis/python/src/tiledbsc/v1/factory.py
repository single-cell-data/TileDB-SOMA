from typing import Union

import tiledb

from .soma_collection import SOMACollection as SOMACollection
from .soma_dataframe import SOMADataFrame as SOMADataFrame
from .soma_dense_nd_array import SOMADenseNdArray as SOMADenseNdArray
from .soma_experiment import SOMAExperiment as SOMAExperiment
from .soma_measurement import SOMAMeasurement as SOMAMeasurement
from .soma_sparse_nd_array import SOMASparseNdArray as SOMASparseNdArray
from .tiledb_group import TileDBGroup
from .util import SOMA_OBJECT_TYPE_METADATA_KEY

MemberType = Union[
    SOMAExperiment,
    SOMAMeasurement,
    SOMACollection,
    SOMADataFrame,
    SOMADenseNdArray,
    SOMASparseNdArray,
]


# TODO: temp class name
def _construct_member(
    member_uri: str, temp_class_name: str, parent: TileDBGroup
) -> MemberType:
    """
    TODO: COMMENT
    """
    # TODO: xref to TileDBObject _set_object_type_metadata/_get_object_type_metadata.
    # and/or, put some of this there as a class/static method.

    # sketch:
    # get class name from meta -- with due respect for:
    # * is-array vs is-group
    # * cloud-ops-count minimization

    class_name = None
    object_type = tiledb.object_type(member_uri)
    if object_type is None:
        raise Exception(f"URI {member_uri} not found")
    elif object_type == "array":
        with tiledb.open(member_uri) as A:  # TODO: CTX
            class_name = A.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
    elif object_type == "group":
        with tiledb.Group(member_uri, mode="r") as G:  # TODO: CTX
            class_name = G.meta[SOMA_OBJECT_TYPE_METADATA_KEY]
    else:
        raise Exception(f"object type {object_type} unrecognized")
    assert class_name is not None

    if class_name == "SOMAExperiment":
        return SOMAExperiment(uri=member_uri, parent=parent)
    elif class_name == "SOMAMeasurement":
        return SOMAMeasurement(uri=member_uri, parent=parent)
    elif class_name == "SOMACollection":
        return SOMACollection(uri=member_uri, parent=parent)
    elif class_name == "SOMADataFrame":
        return SOMADataFrame(uri=member_uri, parent=parent)
    elif class_name == "SOMADenseNdArray":
        return SOMADenseNdArray(uri=member_uri, parent=parent)
    elif class_name == "SOMASparseNdArray":
        return SOMASparseNdArray(uri=member_uri, parent=parent)
    else:
        raise Exception(f'class name "{class_name}" unrecognized')
