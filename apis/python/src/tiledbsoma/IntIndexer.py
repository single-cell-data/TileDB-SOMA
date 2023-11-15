import numpy as np

from . import pytiledbsoma as clib
from ._experiment import _indexer_map_locations

"""External API"""


class IntIndexer:
    @staticmethod
    def map_locations(
        keys: np.typing.NDArray[np.int64],
    ) -> clib.IntIndexer:  # type ignore
        return _indexer_map_locations(keys)


"""For internal use"""
