# from typing import List, Optional
from typing import Optional

#
# from .logging import log_io
# from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup

# from typing import Optional, Sequence, Set, Tuple

# import numpy as np
# import pandas as pd
# import tiledb

# import tiledbsc.v1.util as util


# from .types import Ids


class SOMACollection(TileDBGroup):
    """
    Contains a key-value mapping where the keys are string names and the values are any SOMA-defined
    foundational or composed type, including SOMACollection, SOMADataFrame, SOMADenseNdArray,
    SOMASparseNdArray or SOMAExperiment.
    """

    # create(uri)
    # Create a SOMACollection named with the URI.

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """

        super().__init__(uri=uri, name=name, parent=parent)


# delete(uri)
# Delete the SOMACollection specified with the URI.

#    # TODO: static/class method?
#    #    def delete(uri: str) -> None
#    #        """
#    #        Delete the SOMADataFrame specified with the URI.
#    #        """

# exists(uri) -> bool
# Return true if object exists and is a SOMACollection.

#    # TODO: static/class method?
#    #    def exists(uri: str) -> bool
#    #        """
#    #        Return true if object exists and is a SOMADataFrame.
#    #        """

# get metadata
# Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)

#    #    def get_metadata():
#    #        """
#    #        Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)
#    #        """

#    # get_type() is inherited from TileDBObject

# get(string key)
# Get the object associated with the key

# has(string key)
# Test for the existence of key in collection.

# set(string key, ValueType value)
# Set the key/value in the collection.

# del(string key)
# Remove the key/value from the collection.

# iterator
# Iterate over the collection.
