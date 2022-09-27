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
    from ._version import __version__

    return str(__version__)


def get_storage_engine() -> str:
    """
    Return underlying storage engine name, e.g., "tiledb".
    """
    return "tiledb"
