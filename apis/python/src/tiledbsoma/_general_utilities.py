# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""General utility functions."""

import importlib.metadata
import platform
import sys
from re import fullmatch

from .pytiledbsoma import expected_tiledb_version, tiledb_version
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


def show_package_versions() -> None:
    """Nominal use is for bug reports, so issue filers and issue fixers can be on
    the same page.

    Lifecycle: Maturing.
    """
    u = platform.uname()
    # fmt: off
    print("tiledbsoma.__version__             ", get_implementation_version())  # noqa: T201
    print("TileDB core version (libtiledbsoma)", get_libtiledbsoma_core_version())  # noqa: T201
    print("python version                     ", ".".join(str(v) for v in sys.version_info))  # noqa: T201
    print("OS version                         ", u.system, u.release)  # noqa: T201
    # fmt: on


def _verify_expected_tiledb_version() -> None:
    expected = expected_tiledb_version()
    found = tiledb_version()
    if found != expected:
        raise RuntimeError(
            f"TileDB version mismatch - expected version {expected}, but found {found}. This should not occur, and"
            " is likely the result of a corrupted package installation. Recommend uninstalling/reinstalling the"
            " tiledbsoma package. Alternatively, if you are using a Python virtual environment (e.g., conda)"
            " remove and reinstall the Python virtual environment.",
        )
