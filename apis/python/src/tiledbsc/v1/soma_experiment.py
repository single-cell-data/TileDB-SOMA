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


class SOMAExperiment(SOMACollection):
    """
    TBD
    """

    pass


# Field name
# Field type
# Field description
#
# `obs`
# `SOMADataFrame`
# Primary annotations on the _observation_ axis. The contents of the `__rowid` pseudo-column define
# the _observation_ index domain, aka `obsid`. All observations for the SOMAExperiment _must_ be
# defined in this dataframe.
#
# `ms`
# `SOMACollection[string, SOMAMeasurement]`
# A collection of named measurements.
