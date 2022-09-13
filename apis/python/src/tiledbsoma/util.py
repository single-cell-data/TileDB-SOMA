import time
from typing import TypeVar

import numpy as np
import pandas as pd
import scipy.sparse as sp

T = TypeVar("T", np.ndarray, pd.Series, pd.DataFrame, sp.spmatrix)

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


def get_start_stamp() -> float:
    """
    Returns information about start time of an event. Nominally float seconds since the epoch,
    but articulated here as being compatible with the format_elapsed function.
    """
    return time.time()


def format_elapsed(start_stamp: float, message: str) -> str:
    """
    Returns the message along with an elapsed-time indicator, with end time relative to start
    start from `get_start_stamp`. Used for annotating elapsed time of a task.
    """
    return "%s TIME %.3f seconds" % (message, time.time() - start_stamp)


def is_local_path(path: str) -> bool:
    if path.startswith("file://"):
        return True
    if "://" in path:
        return False
    return True
