# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Common types used across SOMA query modules."""

from typing import Callable, Union

import numpy as np
import numpy.typing as npt
import pyarrow as pa
from typing_extensions import Protocol

IntegerArray = Union[npt.NDArray[np.int64], pa.IntegerArray]


class IndexLike(Protocol):
    """The basics of what we expect an Index to be.

    This is a basic description of the parts of the ``pandas.Index`` type
    that we use. It is intended as a rough guide so an implementor can know
    that they are probably passing the right "index" type into a function,
    not as a full specification of the types and behavior of ``get_indexer``.
    """

    def get_indexer(self, target: IntegerArray) -> npt.NDArray[np.intp]:
        """Something compatible with Pandas' Index.get_indexer method."""


IndexFactory = Callable[[IntegerArray], IndexLike]
"""Function that builds an index over the given ``IntegerArray``.

This interface is implemented by the callable ``pandas.Index``.
"""
