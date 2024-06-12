# Copyright (c) 2021-2023 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2023 TileDB, Inc.
#
# Licensed under the MIT License.

"""General utility functions.
"""
import os
import platform
import sys
import warnings
from re import fullmatch

import tiledb

from .pytiledbsoma import version as libtiledbsoma_core_version_str


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
    if sys.version_info < (3, 8, 0):
        from pkg_resources import DistributionNotFound, get_distribution

        try:
            return get_distribution("tiledbsoma").version
        except DistributionNotFound:
            return "unknown"
    else:
        import importlib.metadata

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

    Lifecycle: maturing
    """
    return "tiledb"


def get_tiledb_py_core_version() -> str:
    """Returns the version of libtiledb ("core") used by the `tiledb` Python library.

    Lifecycle: maturing
    """
    return ".".join(str(ijk) for ijk in list(tiledb.libtiledb.version()))


def get_libtiledbsoma_core_version() -> str:
    """Returns the version of libtiledb ("core") used by libtiledbsoma.

    Lifecycle: maturing
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


def verify_core_versions() -> None:
    """Verify that the versions of libtiledb used by the `tiledb` Python library and
    libtiledbsoma are the same.

    See discussion on https://github.com/single-cell-data/TileDB-SOMA/issues/1837; this
    will be unnecessary when libtiledbsoma is the only "path to core" (cf.
    https://github.com/single-cell-data/TileDB-SOMA/issues/1632).

    Lifecycle: maturing
    """
    tiledb_py_core_version = get_tiledb_py_core_version()
    libtiledbsoma_core_version = get_libtiledbsoma_core_version()
    if tiledb_py_core_version != libtiledbsoma_core_version:
        msg = "libtiledb versions used by tiledb and libtiledbsoma differ: %s != %s" % (
            tiledb_py_core_version,
            libtiledbsoma_core_version,
        )
        if os.environ.get(TILEDB_CORE_MISMATCHED_VERSIONS_ERROR_LEVEL_VAR) == "err":
            print(msg, file=sys.stderr)
            print(
                f"Continuing, since ${TILEDB_CORE_MISMATCHED_VERSIONS_ERROR_LEVEL_VAR} is set, but it is highly recommended you fix the core version mismatch, as undefined behavior and segfaults can result.",
                file=sys.stderr,
            )
        else:
            raise AssertionError(
                f"libtiledb versions used by tiledb and libtiledbsoma differ: {tiledb_py_core_version} != {libtiledbsoma_core_version}"
            )


def show_package_versions() -> None:
    """Nominal use is for bug reports, so issue filers and issue fixers can be on
    the same page.

    Lifecycle: maturing
    """
    u = platform.uname()
    # fmt: off
    print("tiledbsoma.__version__             ", get_implementation_version())
    print("TileDB-Py version                  ", ".".join(str(v) for v in tiledb.version()))
    print("TileDB core version (tiledb)       ", get_tiledb_py_core_version())
    print("TileDB core version (libtiledbsoma)", get_libtiledbsoma_core_version())
    print("python version                     ", ".".join(str(v) for v in sys.version_info))
    print("OS version                         ", u.system, u.release)
    # fmt: on
    verify_core_versions()
