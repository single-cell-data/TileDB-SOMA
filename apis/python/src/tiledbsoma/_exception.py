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


class AlreadyExistsError(SOMAError):
    """Raised when attempting to create an already existing SOMA object.

    Lifecycle: experimental
    """

    pass


def is_already_exists_error(e: tiledb.TileDBError) -> bool:
    """Given a TileDBError, return true if it indicates the object already exists

    Lifecycle: experimental

    Example:
        try:
            tiledb.Array.create(uri, schema, ctx=ctx)
                ...
        except tiledb.TileDBError as e:
            if is_already_exists_error(e):
                ...
            raise e
    """
    stre = str(e)
    # Local-disk, S3, and TileDB Cloud exceptions all have the substring
    # "already exists". Here we lower-case the exception message just
    # in case someone ever uppercases it on the other end.
    if "already exists" in stre.lower():
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
