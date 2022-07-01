import logging

logger = logging.getLogger("tiledbsc")


def warning():
    """
    Sets tiledbsc logging to a WARNING level. Use `tiledbsc.logging.info()` in notebooks to suppress
    progress indicators for data ingestion.
    """
    _set_level(logging.WARNING)


def info():
    """
    Sets tiledbsc logging to an INFO level. Use `tiledbsc.logging.info()` in notebooks to see
    progress indicators for data ingestion.
    """
    _set_level(logging.INFO)


def debug():
    """
    Sets tiledbsc logging to an DEBUG level. Use `tiledbsc.logging.debug()` in notebooks to see more
    detailed progress indicators for data ingestion.
    """
    _set_level(logging.INFO)


def _set_level(level: int) -> None:
    logger.setLevel(level)
    # Without this check, if someone does `tiledbsc.logging.info()` twice, or
    # `tiledbsc.logging.info()` then `tiledbsc.logging.debug()`, etc., then log messages will appear
    # twice.
    if not logger.hasHandlers():
        logger.addHandler(logging.StreamHandler())
