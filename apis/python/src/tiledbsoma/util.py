import pathlib
import time
import urllib.parse
from typing import Any, List, Optional, Tuple, Type, Union

import somacore
from somacore import options


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


def is_relative_uri(uri: str) -> bool:
    """Detects whether a provided URI is a child relative URI to the parent."""
    return "://" not in uri and not uri.startswith("/")


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


def slice_to_range(
    ids: slice, nonempty_domain: Tuple[int, int]
) -> Optional[Tuple[int, int]]:
    """
    For the interface between ``DataFrame::read`` et al. (Python) and ``SOMAReader`` (C++).
    """
    if ids.step is not None:
        raise ValueError("slice steps are not supported")
    if ids == slice(None):
        return None

    start = ids.start
    stop = ids.stop

    # TODO: with future C++ improvements, move half-slice logic to SOMAReader
    if start is None:
        start = nonempty_domain[0]
    if stop is None:
        stop = nonempty_domain[1]

    # Indexing beyond end of array should "trim" to end of array
    if stop > nonempty_domain[1]:
        stop = nonempty_domain[1]

    if start > stop:
        raise ValueError("slice start must be <= slice stop")
    return (start, stop)


def dense_indices_to_shape(
    coords: options.DenseNDCoords,
    array_shape: Tuple[int, ...],
    result_order: somacore.ResultOrder,
) -> Tuple[int, ...]:
    """
    Given a subarray index specified as a tuple of per-dimension slices or scalars
    (e.g., ``([:], 1, [1:2])``), and the shape of the array, return the shape of
    the subarray. Note that the number of coordinates may be less than or equal
    to the number of dimensions in the array.
    """
    if len(coords) > len(array_shape):
        raise ValueError(
            f"coordinate length ({len(coords)}) must be <="
            f" array dimension count ({len(array_shape)})"
        )

    shape: List[int] = []
    for i, extent in enumerate(array_shape):
        if i < len(coords):
            shape.append(dense_index_to_shape(coords[i], extent))
        else:
            shape.append(extent)

    if result_order is somacore.ResultOrder.ROW_MAJOR:
        return tuple(shape)

    return tuple(reversed(shape))


def dense_index_to_shape(
    coord: Union[int, slice],
    array_length: int,
) -> int:
    """
    Given a subarray per-dimension index specified as a slice or scalar (e.g, ``[:], 1, [1:2]``),
    and the shape of the array in that dimension, return the shape of the subarray in
    that dimension.

    Note that Python slice semantics are right-endpoint-exclusive whereas SOMA slice semantics are
    doubly inclusive.
    """
    if type(coord) is int:
        return 1

    if type(coord) is slice:
        start, stop, step = coord.indices(array_length)
        if step != 1:
            raise ValueError("stepped slice ranges are not supported")
        # This is correct for doubly-inclusive slices which SOMA uses.
        stop = min(stop, array_length - 1)
        return stop - start + 1

    raise TypeError("coordinates must be tuple of int or slice")


def check_type(
    name: str,
    actual_value: Any,
    expected_types: Tuple[Type[Any], ...],
) -> None:
    """
    Verifies the type of an argument, or produces a useful error message.
    """
    if not isinstance(actual_value, expected_types):
        if len(expected_types) == 1:
            raise TypeError(
                f"expected {name} argument to be of type {expected_types[0]}; got {type(actual_value)}"
            )
        raise TypeError(
            f"expected {name} argument to be one of {expected_types!r}; got {type(actual_value)}"
        )


def check_unpartitioned(partitions: Optional[options.ReadPartitions]) -> None:
    """Ensures that we're not being asked for a partitioned read.

    Because we currently don't support partitioned reads, we should reject all
    reads that request partitions to avoid giving the user duplicate data across
    sharded tasks.
    """
    if not partitions or partitions == options.IOfN(0, 1):
        return
    raise ValueError("Paritioned reads are not currently supported")
