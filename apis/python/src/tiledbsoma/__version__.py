# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Version information for tiledbsoma package."""

# This file is generated/updated during the build process
# The actual version comes from ../../version.py

import sys
from pathlib import Path

# Try to import from the version.py in the root
try:
    # Add the apis/python directory to the path
    _root = Path(__file__).parent.parent.parent
    if str(_root) not in sys.path:
        sys.path.insert(0, str(_root))

    from version import get_version

    __version__ = get_version()
except (ImportError, FileNotFoundError, AttributeError):
    # Fallback version if we can't get it from git
    __version__ = "0.0.0+unknown"
