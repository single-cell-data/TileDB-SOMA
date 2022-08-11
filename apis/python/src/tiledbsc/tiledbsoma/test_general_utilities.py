import tiledbsc.tiledbsoma


def test_general_utilities() -> None:
    assert tiledbsc.tiledbsoma.get_storage_engine() == "tiledb"
    assert tiledbsc.tiledbsoma.get_implementation() == "python-tiledb"
