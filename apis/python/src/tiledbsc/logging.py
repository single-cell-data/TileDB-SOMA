import logging
from typing import Optional

logger = logging.getLogger("tiledbsc")


def warning() -> None:
    """
    Sets tiledbsc logging to a WARNING level. Use `tiledbsc.logging.info()` in notebooks to suppress
    progress indicators for data ingestion.
    """
    _set_level(logging.WARNING)


def info() -> None:
    """
    Sets tiledbsc logging to an INFO level. Use `tiledbsc.logging.info()` in notebooks to see
    progress indicators for data ingestion.
    """
    _set_level(logging.INFO)


def debug() -> None:
    """
    Sets tiledbsc logging to an DEBUG level. Use `tiledbsc.logging.debug()` in notebooks to see more
    detailed progress indicators for data ingestion.
    """
    _set_level(logging.DEBUG)


def _set_level(level: int) -> None:
    logger.setLevel(level)
    # Without this check, if someone does `tiledbsc.logging.info()` twice, or
    # `tiledbsc.logging.info()` then `tiledbsc.logging.debug()`, etc., then log messages will appear
    # twice.
    if not logger.hasHandlers():
        logger.addHandler(logging.StreamHandler())


def log_io(info_message: Optional[str], debug_message: str) -> None:
    """
    Data-ingesti timeframes range widely, from seconds to the better part of an hour.  Some folks
    won't want details in the former; some will want details in the latter.  For I/O and for I/O
    only, it's helpfulto print a short message at INFO level, or a different, longer message
    at/beyond DEBUG level.
    """
    if logger.level == logging.INFO:
        if info_message is not None:
            logger.info(info_message)
    elif logger.level <= logging.DEBUG:
        if debug_message is not None:
            logger.debug(debug_message)
