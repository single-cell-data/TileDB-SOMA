from pkg_resources import DistributionNotFound, get_distribution


def get_SOMA_version() -> str:
    """
    Return semver-compatible version of the supported SOMA API.
    """
    return "0.2.0-dev"


def get_implementation() -> str:
    """
    Return the implementation name, e.g., "python-tiledb".
    """
    return "python-tiledb"


def get_implementation_version() -> str:
    """
    Return the package implementation version as a semver
    """
    try:
        return get_distribution("tiledbsoma").version
    except DistributionNotFound:
        return "unknown"


def get_storage_engine() -> str:
    """
    Return underlying storage engine name, e.g., "tiledb".
    """
    return "tiledb"
