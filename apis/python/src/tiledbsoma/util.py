import pathlib
import time
import urllib.parse
from typing import Optional, Tuple, TypeVar

import numpy as np
import pandas as pd
import scipy.sparse as sp

T = TypeVar("T", np.ndarray, pd.Series, pd.DataFrame, sp.spmatrix)

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


def get_start_stamp() -> float:
    """
    Returns information about start time of an event. Nominally float seconds since the epoch, but articulated here as being compatible with the format_elapsed function.
    """
    return time.time()


def format_elapsed(start_stamp: float, message: str) -> str:
    """
    Returns the message along with an elapsed-time indicator, with end time relative to start start from ``get_start_stamp``. Used for annotating elapsed time of a task.
    """
    return "%s TIME %.3f seconds" % (message, time.time() - start_stamp)


def is_local_path(path: str) -> bool:
    if path.startswith("file://"):
        return True
    if "://" in path:
        return False
    return True


def make_relative_path(uri: str, relative_to: str) -> str:
    """
    Return a URI relative to another URI. If not possible, raise a ValueError.

    This function assumes that the URI scheme follows posix path conventions
    and only contains a scheme, netloc and path. It does not handle query params,
    fragments, etc.

    The default scheme, if one is not specified, is assumed to be `file`
    """
    p_uri = urllib.parse.urlparse(uri)
    p_relative_to = urllib.parse.urlparse(relative_to)

    uri_scheme = p_uri.scheme if p_uri.scheme != "" else "file"
    relative_to_scheme = p_relative_to.scheme if p_relative_to.scheme != "" else "file"
    if uri_scheme != relative_to_scheme:
        raise ValueError(
            "Unable to make relative path between URIs with different scheme"
        )

    relpath = pathlib.PurePath(p_uri.path).relative_to(p_relative_to.path).as_posix()
    return relpath


def uri_joinpath(base: str, path: str) -> str:
    """
    Join a path to a URI.

    Supports relative paths for `file` or unspecified schemes, assuming
    they are file system paths.  Assumes NO suport for relative paths
    otherwise.
    """
    p_base = urllib.parse.urlparse(base)
    parts = [*p_base]

    if len(path) == 0:
        return base

    if p_base.scheme == "" or p_base.scheme == "file":
        # if a file path, just use pathlib.
        parts[2] = pathlib.PurePath(p_base.path).joinpath(path).as_posix()
    else:
        if ".." in path:
            raise ValueError("Relative paths unsupported")
        if path.startswith("/"):
            # if absolute, just use the path
            parts[2] = path
        else:
            # join, being careful about extraneous path sep
            if parts[2].endswith("/"):
                parts[2] = parts[2] + path
            else:
                parts[2] = parts[2] + "/" + path

    return urllib.parse.urlunparse(parts)


def slice_to_range(ids: slice) -> Optional[Tuple[int, int]]:
    """
    For the interface between ``SOMADataFrame::read`` et al. (Python) and ``SOMAReader`` (C++).
    """
    if ids.start is None and ids.stop is None:
        return None
    if ids.start is None or ids.stop is None:
        # TODO: https://github.com/single-cell-data/TileDB-SOMA/issues/457
        raise ValueError("slice start and stop must be both specified, or neither")
    if ids.start > ids.stop:
        raise ValueError("slice start must be <= slice stop")
    if ids.step is not None and ids.step != 1:
        raise ValueError("slice step must be 1 or None")
    return (ids.start, ids.stop)  # XXX TEMP
