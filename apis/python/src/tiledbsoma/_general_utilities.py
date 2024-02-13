# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""General utility functions.
"""

import platform
import sys

import tiledb
import importlib.metadata

from .pytiledbsoma import version as libtiledbsoma_version


def get_SOMA_version() -> str:
    """Returns semver-compatible version of the supported SOMA API.

    Lifecycle: maturing
    """
    return "0.2.0-dev"


def get_implementation() -> str:
    """Returns the implementation name, e.g., "python-tiledb".

    Lifecycle: maturing
    """
    return "python-tiledb"


def get_implementation_version() -> str:
    """Returns the package implementation version as a semver.

    Lifecycle: maturing
    """
    try:
        return importlib_metadata.version("tiledbsoma")
    except importlib_metadata.PackageNotFoundError
        return "unknown"


def get_storage_engine() -> str:
    """Returns underlying storage engine name, e.g., "tiledb".

    Lifecycle: maturing
    """
    return "tiledb"


def show_package_versions() -> None:
    """Nominal use is for bug reports, so issue filers and issue fixers can be on
    the same page.

    Lifecycle: maturing
    """
    print("tiledbsoma.__version__       ", get_implementation_version())
    print("TileDB-Py tiledb.version()   ", tiledb.version())
    print(
        "TileDB core version          ",
        ".".join(str(ijk) for ijk in list(tiledb.libtiledb.version())),
    )
    print("libtiledbsoma version()      ", libtiledbsoma_version())
    print("python version               ", ".".join(str(v) for v in sys.version_info))
    u = platform.uname()
    print("OS version                   ", u.system, u.release)
