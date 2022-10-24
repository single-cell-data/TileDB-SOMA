from pathlib import Path

import tiledb

import tiledbsoma
import tiledbsoma.io
import tiledbsoma.util

HERE = Path(__file__).parent


def test_import_anndata(tmp_path):

    ann1 = HERE.parent / "anndata/pbmc-small.h5ad"
    ann2 = HERE.parent / "anndata/pbmc-small-x-csr.h5ad"
    ann3 = HERE.parent / "anndata/pbmc-small-x-csc.h5ad"

    soco_dir = tmp_path.as_posix()
    soma1_dir = (tmp_path / "soma1").as_posix()
    soma2_dir = (tmp_path / "soma2").as_posix()
    soma3_dir = (tmp_path / "soma3").as_posix()

    soma1 = tiledbsoma.SOMA(soma1_dir, name="soma1")
    tiledbsoma.io.from_h5ad(soma1, ann1)

    soma2 = tiledbsoma.SOMA(soma2_dir, name="soma2")
    tiledbsoma.io.from_h5ad(soma2, ann2)

    soma3 = tiledbsoma.SOMA(soma3_dir, name="soma3")
    tiledbsoma.io.from_h5ad(soma3, ann3)

    soco = tiledbsoma.SOMACollection(soco_dir)

    with tiledb.Group(soma1_dir) as G:
        assert G.meta[tiledbsoma.util.SOMA_OBJECT_TYPE_METADATA_KEY] == "SOMA"

    soco.create_unless_exists()
    assert len(soco.get_member_names()) == 0

    soco.add(soma1)
    assert len(soco.get_member_names()) == 1
    soco.add(soma2)
    assert len(soco.get_member_names()) == 2
    soco.add(soma3)
    assert len(soco.get_member_names()) == 3

    soco.remove(soma1)
    assert len(soco.get_member_names()) == 2
    del soco["soma2"]
    assert len(soco.get_member_names()) == 1
    del soco.soma3
    assert len(soco.get_member_names()) == 0

    assert tiledbsoma.util.is_soma(soma1.uri)
    assert tiledbsoma.util.is_soma(soma2.uri)
    assert not tiledbsoma.util.is_soma(soma1.obs.uri)
    assert not tiledbsoma.util.is_soma(soma2.var.uri)

    assert not tiledbsoma.util.is_soma_collection(soma2.var.uri)
    assert not tiledbsoma.util.is_soma_collection(soma3.uri)
    assert tiledbsoma.util.is_soma_collection(soco.uri)
