# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Exceptions.
"""

from typing import Union


class SOMAError(Exception):
    """Base error type for SOMA-specific exceptions.

    Lifecycle: Maturing.
    """

    pass


class DoesNotExistError(SOMAError):
    """Raised when attempting to open a non-existent or inaccessible SOMA object.

    Lifecycle: Maturing.
    """

    pass


def is_does_not_exist_error(e: Union[RuntimeError, SOMAError]) -> bool:
    """Given a RuntimeError or SOMAError, return true if it indicates the object
    does not exist

    Lifecycle: Maturing.

    Example:
        try:
            with tiledbsoma.open(uri):
                ...
        except (RuntimeError, SOMAError) as e:
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
        or "[SOMAObject::open] " in stre
    ):
        return True

    return False


class AlreadyExistsError(SOMAError):
    """Raised when attempting to create an already existing SOMA object.

    Lifecycle: Maturing.
    """

    pass


def is_already_exists_error(e: SOMAError) -> bool:
    """Given a SOMAError, return true if it indicates the object already exists

    Lifecycle: Maturing.

    Example:
        try:
            clib.SOMASparseNDArray.create(...)
                ...
        except SOMAError as e:
            if is_already_exists_error(e):
                ...
            raise e
    """
    stre = str(e)
    # Local-disk, S3, and TileDB Cloud exceptions all have the substring
    # "already exists". Here we lower-case the exception message just
    # in case someone ever uppercases it on the other end.
    return "already exists" in stre.lower()


class NotCreateableError(SOMAError):
    """Raised when attempting to create an already existing SOMA object.

    Lifecycle: Maturing
    """

    pass


def is_duplicate_group_key_error(e: SOMAError) -> bool:
    """Given a TileDBError, return try if it indicates a duplicate member
    add request in a tiledb.Group.

    Lifecycle: Maturing
    """
    return "member already exists in group" in str(e)


def is_domain_setting_error(e: SOMAError) -> bool:
    """Given a SOMAError, return whether it attempted to create
    the ArraySchema but the passed in domain was invalid.

    Lifecycle: Maturing
    """
    return "Cannot set domain" in str(e)


def map_exception_for_create(e: SOMAError, uri: str) -> Exception:
    if is_already_exists_error(e):
        return AlreadyExistsError(f"{uri!r} already exists")
    if is_domain_setting_error(e):
        return ValueError(e)
    return e
