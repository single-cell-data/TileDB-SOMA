import tiledbsc.v1


def test_general_utilities():
    assert tiledbsc.v1.get_storage_engine() == "tiledb"
    assert tiledbsc.v1.get_implementation() == "python-tiledb"
