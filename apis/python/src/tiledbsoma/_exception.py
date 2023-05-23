# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Exceptions.
"""

import tiledb


class SOMAError(Exception):
    """Base error type for SOMA-specific exceptions.

    Lifecycle: maturing
    """

    pass


class DoesNotExistError(SOMAError):
    """Raised when attempting to open a non-existent or inaccessible SOMA object.

    Lifecycle: maturing
    """

    pass


def is_does_not_exist_error(e: tiledb.TileDBError) -> bool:
    """Given a TileDBError, return true if it indicates the object does not exist

    Lifecycle: maturing

    Example:
        try:
            with tiledb.open(uri):
                ...
        except tiledb.TileDBError as e:
            if is_does_not_exist_error(e):
                ...
            raise e
    """
    stre = str(e)
    # Local-disk/S3 does-not-exist exceptions say 'Group does not exist'; TileDB Cloud
    # does-not-exist exceptions are worded less clearly.
    if (
        "does not exist" in stre
        or "Unrecognized array" in stre
        or "HTTP code 401" in stre
        or "HTTP code 404" in stre
    ):
        return True

    return False


def is_duplicate_group_key_error(e: tiledb.TileDBError) -> bool:
    """Given a TileDBError, return try if it indicates a duplicate member
    add request in a tiledb.Group.

    Lifecycle: maturing
    """
    stre = str(e)
    if "member already exists in group" in stre:
        return True

    return False
