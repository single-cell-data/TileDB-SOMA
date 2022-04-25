import os
import sys
import pathlib
import tempfile
import pytest

import tiledb

from .. import anndata_to_tiledb

_file_path = os.path.dirname(__file__)

def test_anndata_import(tmp_path):
    #array_tmp = tempfile.mkdtemp()
    array_tmp = tmp_path
    data_path = pathlib.Path(_file_path) / "../../data/10x-pbmc-multiome-v1.0/subset_100_100.h5ad"

    exporter = anndata_to_tiledb.AnnDataExporter.from_h5ad(str(data_path))
    exporter.write_tiledb(array_tmp)

    array_path = pathlib.Path(array_tmp)
    X_path = array_path / "X"
    var_path = array_path / "var"
    obs_path = array_path / "obs"
    assert tiledb.object_type(str(X_path)) == "array"
    assert tiledb.object_type(str(obs_path)) == "array"
    assert tiledb.object_type(str(var_path)) == "array"