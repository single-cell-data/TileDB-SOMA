import anndata
import tiledb
import tiledbsc
import tiledbsc.io

import pytest
import tempfile
import os
from pathlib import Path

HERE = Path(__file__).parent


def test_import_anndata(tmp_path):

    ann1 = HERE.parent / "anndata/pbmc-small.h5ad"
    ann2 = HERE.parent / "anndata/pbmc-small-x-csr.h5ad"

    soco_dir = tmp_path.as_posix()
    soma1_dir = (tmp_path / "soma1").as_posix()
    soma2_dir = (tmp_path / "soma2").as_posix()

    soma1 = tiledbsc.SOMA(soma1_dir, name="soma1")
    tiledbsc.io.from_h5ad(soma1, ann1)

    soma2 = tiledbsc.SOMA(soma2_dir, name="soma2")
    tiledbsc.io.from_h5ad(soma2, ann2)

    soco = tiledbsc.SOMACollection(soco_dir)

    with tiledb.Group(soma1_dir) as G:
        assert G.meta[tiledbsc.util_tiledb.SOMA_OBJECT_TYPE_METADATA_KEY] == "SOMA"

    soco._create()
    assert len(soco._get_member_names()) == 0

    soco.add(soma1)
    assert len(soco._get_member_names()) == 1
    soco.add(soma2)
    assert len(soco._get_member_names()) == 2

    soco.remove(soma1)
    assert len(soco._get_member_names()) == 1
    soco.remove(soma2)
    assert len(soco._get_member_names()) == 0
