from dataclasses import dataclass
from typing import Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
import pyarrow as pa

AxisCoordinate = Union[slice, int, npt.ArrayLike]
AxisCoordinates = Tuple[AxisCoordinate, ...]
AxisValueFilter = str


@dataclass()
class AxisQuery:
    """
    A dataclass which defines a single-axis dataframe query based upon coordinates and/or a
    value filter predicate [lifecycle: experimental].

    Per dimension, the AxisQuery can have value of:
    * None - all data
    * Coordinates - a set of coordinates on the axis dataframe index (or soma_rowids if
      a dense dataframe), expressed as an int, a slice or a NumPy "array like" of ints.
    * A SOMA `value_filter` across columns in the axis dataframe, expressed as string
    * Or, a combination of coordinates and value filter.

    Parameters
    ----------
    coords : Optional[tuple]
        Query (slice) by dimension. Tuple must have lenth less than or equal to the number
        of dimensions. For each dimension, coordinate is an ``int``, ``slice`` or NumPy
        ``array-like`` of ints.
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

    value_filter: Optional[AxisValueFilter] = None
    coords: Optional[AxisCoordinates] = None

    def __post_init__(self) -> None:
        if self.value_filter is not None:
            if not isinstance(self.value_filter, str):
                raise TypeError("AxisQuery - value_filter must be a str type")

        if self.coords is None:
            self.coords = (slice(None),)
        else:
            if not isinstance(self.coords, tuple):
                raise TypeError(
                    "AxisQuery - coords must be tuple of int, slice or numpy.array_like"
                )
            coords = []
            for c in self.coords:
                if isinstance(c, int) or isinstance(c, slice):
                    coords.append(c)
                else:
                    coords.append(pa.array(np.array(c, dtype=np.int64)))
            self.coords = tuple(coords)
