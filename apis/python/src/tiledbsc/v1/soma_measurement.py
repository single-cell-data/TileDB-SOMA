# from typing import List, Optional

# from .logging import log_io
from .soma_collection import SOMACollection

# from .tiledb_group import TileDBGroup

# from typing import Optional, Sequence, Set, Tuple

# import numpy as np
# import pandas as pd
# import tiledb

# import tiledbsc.v1.util as util


# from .types import Ids


class SOMAMeasurement(SOMACollection):
    """
    TBD
    """

    pass


# The `SOMAMeasurement` is a sub-element of a SOMAExperiment, and is otherwise a specialized SOMACollection with pre-defined fields:
#
# Field name
# Field type
# Field description
#
# `var`
# `SOMADataFrame`
# Primary annotations on the _variable_ axis, for variables in this measurement (i.e., annotates
# columns of `X`). The contents of the `__rowid` pseudo-column define the _variable_ index domain,
# aka `varid`. All variables for this measurement _must_ be defined in this dataframe.

# `X`
# `SOMACollection[string, SOMASparseNdArray]`
# A collection of sparse matrices, each containing measured feature values. Each matrix is indexed
# by `[obsid, varid]`

# `obsm`
# `SOMACollection[string, SOMADenseNdArray]`
# A collection of dense matrices containing annotations of each _obs_ row. Has the same shape as
# `obs`, and is indexed with `obsid`.

# `obsp`
# `SOMACollection[string, SOMASparseNdArray]`
# A collection of sparse matrices containing pairwise annotations of each _obs_ row. Indexed with
# `[obsid_1, obsid_2].`

# `varm`
# `SOMACollection[string, SOMADenseNdArray]`
# A collection of dense matrices containing annotations of each _var_ row. Has the same shape as
# `var`, and is indexed with `varid`

# `varp`
# `SOMACollection[string, SOMASparseNdArray]`
# A collection of sparse matrices containing pairwise annotations of each _var_ row. Indexed with
# `[varid_1, varid_2]`
#
# For the entire `SOMAExperiment`, the index domain for the elements within `obsp`, `obsm` and `X`
# (first dimension) are the values defined by the `obs` `SOMADataFrame` `__rowid` column. For each
# `SOMAMeasurement`, the index domain for `varp`, `varm` and `X` (second dimension) are the values
# defined by the `var` `SOMADataFrame` `__rowid` column in the same measurement. In other words, all
# predefined fields in the `SOMAMeasurement` share a common `obsid` and `varid` domain, which is
# defined by the contents of the respective columns in `obs` and `var` SOMADataFrames.
#
# As with other SOMACollections, the `SOMAExperiment` and `SOMAMeasurement` also have a `metadata`
# field, and may contain other user-defined elements. Keys in a `SOMAExperiment` and
# `SOMAMeasurement` beginning with the characters `_`, `.`, or `$` are reserved for ad hoc use, and
# will not be utilized by this specification. All other keys are reserved for future specifications.
#
# The following naming and indexing constraints are defined for the `SOMAExperiment` and `SOMAMeasurement`:
#
# Field name
# Field constraints
#
# `obs`, `var`
# Field type is a `SOMADataFrame`

# `obsp`, `varp`, `X`
# Field type is a `SOMACollection`, and each element in the collection has a value of type `SOMASparseNdArray`

# `obsm`, `varm`
# Field type is a `SOMACollection`, and each element in the collection has a value of type `SOMADenseNdArray`

# `obsm`, `obsp`, `varm`, `varp`
# Fields may be empty collections.

# `X` collection values
# All matrices must have the shape `(#obs, #var)`. The domain of the first dimension is the values
# of `obs.__rowid`, and the index domain of the second dimension is the values of `var.__rowid` in
# the containing `SOMAMeasurement`.

# `obsm` collection values
# All matrices must have the shape `(#obs, M)`, where `M` is user-defined. The domain of the first
# dimension is the values of `obs.__rowid`.

# `obsp` collection values
# All matrices must have the shape `(#obs, #obs)`. The domain of both dimensions is the values of
# `obs.__rowid`.

# `varm` collection values
# All matrices must have the shape `(#var, M)`, where `M` is user-defined. The domain of the first
# dimension is the values of `var.__rowid`.

# `varp` collection values
# All matrices must have the shape `(#var, #var)`. The domain of both dimensions is the values of
# `var.__rowid`.
