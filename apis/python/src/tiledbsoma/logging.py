# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

import logging
import sys
from typing import Optional


def warning() -> None:
    """Sets ``tiledbsoma.logging`` to a WARNING level.
    Use ``tiledbsoma.logging.info()`` in notebooks to suppress
    progress indicators for data ingestion.

    Lifecycle:
        Experimental.
    """
    _set_level(logging.WARNING)


def info() -> None:
    """Sets ``tiledbsoma.logging`` to an INFO level.
    Use ``tiledbsoma.logging.info()`` in notebooks to see
    progress indicators for data ingestion.

    Lifecycle:
        Experimental.
    """
    _set_level(logging.INFO)


def debug() -> None:
    """Sets ``tiledbsoma.logging`` to an DEBUG level.
    Use ``tiledbsoma.logging.debug()`` in notebooks to see
    more detailed progress indicators for data ingestion.

    Lifecycle:
        Experimental.
    """
    _set_level(logging.DEBUG)


def _set_level(level: int) -> None:
    logger.setLevel(level)
    # Without this check, if someone does ``tiledbsoma.logging.info()`` twice, or
    # ``tiledbsoma.logging.info()`` then ``tiledbsoma.logging.debug()``, etc., then log messages will appear
    # twice.
    if not logger.hasHandlers():
        logger.addHandler(logging.StreamHandler(stream=sys.stdout))


def log_io(info_message: Optional[str], debug_message: str) -> None:
    """Data-ingestion timeframes range widely.
    Some folks won't want details for smaller uploads; some will want details for larger ones.
    For I/O and for I/O only, it's helpful to print a short message at INFO level,
    or a different, longer message at/beyond DEBUG level.

    Lifecycle:
        Experimental.
    """
    if logger.level == logging.INFO:
        if info_message is not None:
            logger.info(info_message)
    elif logger.level <= logging.DEBUG:
        if debug_message is not None:
            logger.debug(debug_message)


def _in_notebook() -> bool:
    try:
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


logger = logging.getLogger("tiledbsoma")
# For interactive notebook use it's _crucial_ that data ingests (which can take several minutes)
# must make a progress indicator to show something is happening -- by default and with zero user
# intervention.
#
# Additional note: without this first _set_level, default logging will go to stderr. We want it to
# stdout.
if _in_notebook():
    _set_level(logging.INFO)
else:
    _set_level(logging.WARN)
