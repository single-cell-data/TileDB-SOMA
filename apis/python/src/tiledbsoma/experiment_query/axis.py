import collections.abc
from typing import Any, Optional

import attrs
import numpy as np
import pyarrow as pa

from ..types import SparseDataFrameCoordinates


@attrs.define(kw_only=True)
class AxisQuery:
    """
    A class which defines a single-axis dataframe query based upon coordinates and/or a
    value filter predicate [lifecycle: experimental].

    Per dimension, the AxisQuery can have value of:
    * None - all data
    * Coordinates - a set of coordinates on the axis dataframe index, expressed in any
      type or format supported by ``DataFrame.read()``.
    * A SOMA `value_filter` across columns in the axis dataframe, expressed as string
    * Or, a combination of coordinates and value filter.

    Parameters
    ----------
    coords : Optional[SparseDataFrameCoordinates]
        Query (slice) by dimension. Tuple must have lenth less than or equal to the number
        of dimensions, and be of a type supported by ``DataFrame``.
    value_filter : Optional[str]
        A string specifying a SOMA value_filter.

    Examples
    --------
    ```
        AxisQuery()  # all data
        AxisQuery(coords=(slice(1,10),))  # 1D, slice
        AxisQuery(coords=([0,1,2]))  # 1D, point indexing using array-like
        AxisQuery(coords=(slice(None), numpy.array([0,88,1001])))  # 2D
        AxisQuery(value_filter="tissue == 'lung'")
        AxisQuery(coords=(slice(1,None),), value_filter="tissue == 'lung'")
    ```
    """

    value_filter: Optional[str] = attrs.field(
        default=None,
        validator=attrs.validators.optional(attrs.validators.instance_of(str)),
    )
    coords: SparseDataFrameCoordinates = attrs.field(
        default=(slice(None),),
        validator=attrs.validators.instance_of(tuple),
    )

    @coords.validator
    def _validate_coords(self, _: Any, value: Any) -> None:
        """
        This should arguably be delegated to DataClass, but that would
        defer error reporting to the user. Proactive validation has
        UX benefit.
        """
        for c in value:
            if c is None:
                continue
            if isinstance(
                c,
                (
                    int,
                    slice,
                    collections.abc.Sequence,
                    pa.Array,
                    pa.ChunkedArray,
                    np.ndarray,
                ),
            ):
                continue

            raise TypeError("AxisQuery coordinate type is unsupported.")
