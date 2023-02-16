import sys

import tiledb
from pkg_resources import DistributionNotFound, get_distribution


def get_SOMA_version() -> str:
    """
    Return semver-compatible version of the supported SOMA API.

    [lifecycle: experimental]
    """
    return "0.2.0-dev"


def get_implementation() -> str:
    """
    Return the implementation name, e.g., "python-tiledb".

    [lifecycle: experimental]
    """
    return "python-tiledb"


def get_implementation_version() -> str:
    """
    Return the package implementation version as a semver

    [lifecycle: experimental]
    """
    try:
        return get_distribution("tiledbsoma").version
    except DistributionNotFound:
        return "unknown"


def get_storage_engine() -> str:
    """
    Return underlying storage engine name, e.g., "tiledb".

    [lifecycle: experimental]
    """
    return "tiledb"


def show_package_versions() -> None:
    """
    Nominal use is for bug reports, so issue filers and issue fixers can be on the same page.
    """
    print("tiledbsoma.__version__   ", get_implementation_version())
    print("tiledb.__version__       ", tiledb.__version__)
    print(
        "core version             ",
        ".".join(str(ijk) for ijk in list(tiledb.libtiledb.version())),
    )
    print("python version           ", ".".join(str(v) for v in sys.version_info))
