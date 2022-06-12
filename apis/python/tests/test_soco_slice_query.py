import anndata
import tiledb
import tiledbsc
import tiledbsc.io

import pytest
import tempfile
import os
from pathlib import Path

HERE = Path(__file__).parent


def test_soco_slice_query(tmp_path):
    # Populate the collection

    soco_dir = tmp_path.as_posix()

    soco = tiledbsc.SOMACollection(soco_dir)
    soco.create_unless_exists()

    for name, h5ad_path in [
        ("subset-soma-01", HERE.parent / "anndata/subset-soma-01.h5ad"),
        ("subset-soma-02", HERE.parent / "anndata/subset-soma-02.h5ad"),
        ("subset-soma-03", HERE.parent / "anndata/subset-soma-03.h5ad"),
        ("subset-soma-04", HERE.parent / "anndata/subset-soma-04.h5ad"),
    ]:
        soma_path = os.path.join(soco_dir, name)
        print("Writing", soma_path, "...")
        soma = tiledbsc.SOMA(soma_path, verbose=False)
        tiledbsc.io.from_h5ad(soma, h5ad_path)
        soco.add(soma)
    print("Done.")

    # Do the slice query
    obs_attr_names = ["tissue"]
    obs_query_string = 'tissue == "blood"'
    var_attr_names = ["feature_name"]
    var_query_string = 'feature_name == "MT-CO3"'

    soma_slices = []
    for soma in soco:
        # E.g. querying for 'cell_type == "blood"' but this SOMA doesn'tiledbsc have a cell_type column in
        # its obs at all.
        if not soma.obs.has_attr_names(obs_attr_names):
            continue
        # E.g. querying for 'feature_name == "MT-CO3"' but this SOMA doesn'tiledbsc have a feature_name
        # column in its var at all.
        if not soma.var.has_attr_names(var_attr_names):
            continue

        soma_slice = soma.attribute_filter(obs_query_string, var_query_string)
        if soma_slice != None:
            print("...", soma.name, soma_slice.ann.X.shape)
            soma_slices.append(soma_slice)

    result_soma_slice = tiledbsc.SOMASlice.concat(soma_slices)
    assert result_soma_slice is not None

    ann = result_soma_slice.to_anndata()

    assert ann.obs.shape == (400, 17)
    assert ann.var.shape == (1, 0)
    assert ann.X.shape == (400, 1)
