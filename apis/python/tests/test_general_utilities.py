from unittest import mock

import pytest

import tiledbsoma
import tiledbsoma.pytiledbsoma as clib


def test_general_utilities():
    assert tiledbsoma.get_storage_engine() == "tiledb"
    assert tiledbsoma.get_implementation() == "python-tiledb"


def test_versions_api():
    assert tiledbsoma.get_SOMA_version() == "0.2.0-dev"
    assert tiledbsoma.get_implementation_version() == tiledbsoma.__version__


def test_expected_version_api():
    assert clib.tiledb_version() == clib.expected_tiledb_version()


def test_verify_expected_tiledb_version() -> None:
    from tiledbsoma._general_utilities import _verify_expected_tiledb_version

    with mock.patch("tiledbsoma._general_utilities.expected_tiledb_version") as mock_expected_tiledb_version:
        mock_expected_tiledb_version.return_value = (1, 2, 3)

        with pytest.raises(RuntimeError):
            _verify_expected_tiledb_version()

        mock_expected_tiledb_version.assert_called_once()
