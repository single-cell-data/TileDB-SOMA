import logging

logger = logging.getLogger(__name__)  # Nominally __name__ is 'tiledbsc'


def info():
    """
    Sets tiledbsc logging to an INFO level. Use `tiledbsc.logging.info()` in notebooks to see
    progress indicators for data ingestion.
    """
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())


def warning():
    """
    Sets tiledbsc logging to a WARNING level. Use `tiledbsc.logging.info()` in notebooks to suppress
    progress indicators for data ingestion.
    """
    logger.setLevel(logging.WARNING)
    logger.addHandler(logging.StreamHandler())
