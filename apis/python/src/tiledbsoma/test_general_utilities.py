import tiledbsoma


def test_general_utilities() -> None:
    assert tiledbsoma.get_storage_engine() == "tiledb"
    assert tiledbsoma.get_implementation() == "python-tiledb"
