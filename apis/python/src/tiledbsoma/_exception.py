# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""Exceptions.
"""


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


def is_does_not_exist_error(e: RuntimeError) -> bool:
    """Given a RuntimeError, return true if it indicates the object does not exist

    Lifecycle: Maturing.

    Example:
        try:
            with tiledbsoma.open(uri):
                ...
        except RuntimeError as e:
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


def is_not_createable_error(e: SOMAError) -> bool:
    """Given a TileDBError, return true if it indicates the object cannot be created

    Lifecycle: Maturing

    Example:
        try:
            tiledb.Array.create(uri, schema, ctx=ctx)
                ...
        except tiledb.TileDBError as e:
            if is_not_createable_error(e):
                ...
            raise e
    """
    stre = str(e)
    # Context:
    # * A recurring paradigm in tiledbsoma.io is open for write (if exists) else create --
    #   or, equivalently, create (if doesn't already exist), else open for write
    # * A priori either seems fine
    # * There are performance implications for trying the create first: when an
    #   object _does_ already exist we get that quickly.
    # * Therefore it's more performant to try-create-catch-open-for-write
    # * However we have the following semantics for cloud URIs:
    #   o For writes: must be "creation URIs" of the form "tiledb://namespace/s3://bucket/some/path"
    #   o For read:   can be "creation URIs" _or_ non-creation URIs of the form
    #     "tiledb://namespace/groupname" or "tiledb://namespace/uuid"
    # * Put together: when we try-create-catch-open-for-write, _and when_ the URI provided
    #   is a non-creation URI, we need to catch that fact and treat it as a non-error.
    stre = stre.lower()
    if "storage backend local not supported" in stre:
        return True
    if "storage backend not supported: local" in stre:
        return True
    return False


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
    if is_not_createable_error(e):
        return NotCreateableError(f"{uri!r} cannot be created")
    if is_domain_setting_error(e):
        return ValueError(e)
    return e
