# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Logging configuration helpers."""

from __future__ import annotations

import logging

logger = logging.getLogger("tiledbsoma")


def warning() -> None:
    """Sets :mod:`tiledbsoma.logging` to a WARNING level.
    Use :func:`tiledbsoma.logging.info` in notebooks to suppress
    progress indicators for data ingestion.

    Lifecycle:
        Maturing.
    """
    _set_level(logging.WARNING)


def info() -> None:
    """Sets :mod:`tiledbsoma.logging` to an INFO level.
    Use :func:`tiledbsoma.logging.info` in notebooks to see
    progress indicators for data ingestion.

    Lifecycle:
        Maturing.
    """
    _set_level(logging.INFO)


def debug() -> None:
    """Sets :mod:`tiledbsoma.logging` to a DEBUG level.
    Use :func:`tiledbsoma.logging.debug` in notebooks to see
    more detailed progress indicators for data ingestion.

    Lifecycle:
        Maturing.
    """
    _set_level(logging.DEBUG)


def _set_level(level: int) -> None:
    logger.setLevel(level)
    # Without this check, if someone does tiledbsoma.logging.info twice, or
    # tiledbsoma.logging.info then tiledbsoma.logging.debug, etc., then log messages will appear
    # twice.
    if not logger.hasHandlers():
        logger.addHandler(logging.StreamHandler())


def log_io_same(message: str) -> None:
    """Log message to both INFO and DEBUG."""
    log_io(message, message)


def log_io(info_message: str | None, debug_message: str) -> None:
    """Data-ingestion timeframes range widely.
    Some folks won't want details for smaller uploads; some will want details for larger ones.
    For I/O and for I/O only, it's helpful to print a short message at INFO level,
    or a different, longer message at/beyond DEBUG level.

    Lifecycle:
        Maturing.
    """
    if logger.level == logging.INFO:
        if info_message is not None:
            logger.info(info_message)
    elif logger.level <= logging.DEBUG and debug_message is not None:
        logger.debug(debug_message)
