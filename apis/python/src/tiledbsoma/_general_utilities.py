# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""General utility functions.
"""
import importlib.metadata
import platform
import sys
import warnings
from re import fullmatch

from .pytiledbsoma import version as libtiledbsoma_core_version_str


def get_SOMA_version() -> str:
    """Returns semver-compatible version of the supported SOMA API.

    Lifecycle: Maturing
    """
    return "0.2.0-dev"


def get_implementation() -> str:
    """Returns the implementation name, e.g., "python-tiledb".

    Lifecycle: Maturing.
    """
    return "python-tiledb"


def get_implementation_version() -> str:
    """Returns the package implementation version as a semver.

    Lifecycle: Maturing.
    """

    try:
        return importlib.metadata.version("tiledbsoma")
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def assert_version_before(major: int, minor: int) -> None:
    version_string = get_implementation_version()
    if version_string == "unknown":
        warnings.warn(
            "`assert_version_before` could not retrieve the current "
            "implementation version"
        )
        return

    version = version_string.split(".")
    assert (int(version[0]), int(version[1])) < (major, minor)


def get_storage_engine() -> str:
    """Returns underlying storage engine name, e.g., "tiledb".

    Lifecycle: Maturing.
    """
    return "tiledb"


def get_libtiledbsoma_core_version() -> str:
    """Returns the version of libtiledb ("core") used by libtiledbsoma.

    Lifecycle: Maturing.
    """
    v = libtiledbsoma_core_version_str()
    m = fullmatch(r"libtiledb=(\d+\.\d+\.\d+)", v)
    if m is None:
        raise ValueError(f"Unexpected libtiledbsoma_core_version: {v}")
    return m.group(1)


# Set this env var to "err" to print an error to stderr when TileDB-Py's and libtiledbsoma's core
# versions mismatch (by default, an AssertionError is raised).
TILEDB_CORE_MISMATCHED_VERSIONS_ERROR_LEVEL_VAR = (
    "TILEDB_CORE_MISMATCHED_VERSIONS_ERROR_LEVEL"
)


def show_package_versions() -> None:
    """Nominal use is for bug reports, so issue filers and issue fixers can be on
    the same page.

    Lifecycle: Maturing.
    """
    u = platform.uname()
    # fmt: off
    print("tiledbsoma.__version__             ", get_implementation_version())
    print("TileDB core version (libtiledbsoma)", get_libtiledbsoma_core_version())
    print("python version                     ", ".".join(str(v) for v in sys.version_info))
    print("OS version                         ", u.system, u.release)
    # fmt: on
