# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Version information for tiledbsoma package."""

# Version is automatically managed by setuptools_scm

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("tiledbsoma")
except PackageNotFoundError:
    # Package is not installed, use fallback
    __version__ = "0.0.0+unknown"
