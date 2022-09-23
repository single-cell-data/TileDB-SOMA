import tiledbsoma


def test_general_utilities():
    assert tiledbsoma.get_storage_engine() == "tiledb"
    assert tiledbsoma.get_implementation() == "python-tiledb"


def test_versions_api():
    assert tiledbsoma.get_SOMA_version() == "0.2.0-dev"
    assert tiledbsoma.get_implementation_version() == tiledbsoma.__version__
