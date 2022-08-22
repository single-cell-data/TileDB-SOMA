import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse

import tiledbsc
import tiledbsc.io

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
        soma = tiledbsc.SOMA(soma_path)
        tiledbsc.io.from_h5ad(soma, h5ad_path)
        soco.add(soma)

    # Do the slice query
    obs_attrs = ["tissue"]
    obs_query_string = 'tissue == "blood"'
    var_attrs = ["feature_name"]
    var_query_string = 'feature_name == "MT-CO3"'

    soma_slices = []
    for soma in soco:
        # E.g. querying for 'cell_type == "blood"' but this SOMA doesn't have a cell_type column in
        # its obs at all.
        if not soma.obs.has_attr_names(obs_attrs):
            continue
        # E.g. querying for 'feature_name == "MT-CO3"' but this SOMA doesn't have a feature_name
        # column in its var at all.
        if not soma.var.has_attr_names(var_attrs):
            continue

        soma_slice = soma.query(
            obs_query_string=obs_query_string, var_query_string=var_query_string
        )
        if soma_slice is not None:
            soma_slices.append(soma_slice)

    result_soma_slice = tiledbsc.SOMASlice.concat(soma_slices)
    assert result_soma_slice is not None

    ann = result_soma_slice.to_anndata()

    assert ann.obs.shape == (400, 17)
    assert ann.var.shape == (1, 3)
    assert ann.X.shape == (400, 1)


def test_soco_slice_query_nans(tmp_path):

    # Example of the problem:
    # ---- SOMA1      # ---- SOMA2      # ---- SLICE (WITH QUERY ALL) AND CONCAT
    # X:              # X:              # X:
    # [[1. 2. 3.]     # [[ 7.  8.]      # [[ 1.  2.  3.    ]
    #  [4. 5. 6.]]    #  [ 9. 10.]]     #  [ 4.  5.  6.    ]
    # obs:            # obs:            #  [ 7.          8.]
    #        oa ob    #        oa ob    #  [ 9.         10.]]
    # cell1  10  a    # cell3  60  f    # obs:
    # cell2  20  b    # cell4  70  g    #        oa ob
    # var:            # var:            # cell1  10  a
    #        va vb    #        va vb    # cell2  20  b
    # gene1  30  c    # gene1  80  h    # cell3  60  f
    # gene2  40  d    # gene4  90  i    # cell4  70  g
    # gene3  50  e                      # var:
    #          va   vb
    # gene1  30.0    c
    # gene2  40.0    d
    # gene3  50.0    e
    # gene4   NaN  NaN <---- make sure this doesn't happen

    # ----------------------------------------------------------------
    # Make AnnData object in memory
    X1 = np.asarray([[1, 2, 3], [4, 5, 6]], dtype=np.float32)
    X1 = scipy.sparse.csr_matrix(X1)
    obs1 = pd.DataFrame(
        index=np.asarray(["cell1", "cell2"]),
        data={
            "oa": np.asarray([10, 20]),
            "ob": np.asarray(["a", "b"]),
        },
    )
    var1 = pd.DataFrame(
        index=np.asarray(["gene1", "gene2", "gene3"]),
        data={
            "va": np.asarray([30, 40, 50]),
            "vb": np.asarray([True, False, True]),
            "vc": np.asarray(["c", "d", "e"]),
        },
    )
    adata1 = ad.AnnData(X=X1, obs=obs1, var=var1)

    X2 = np.asarray([[7, 8], [9, 10]], dtype=np.float32)
    X2 = scipy.sparse.csr_matrix(X2)
    obs2 = pd.DataFrame(
        index=np.asarray(["cell3", "cell4"]),
        data={
            "oa": np.asarray([60, 70]),
            "ob": np.asarray(["f", "g"]),
        },
    )
    var2 = pd.DataFrame(
        index=np.asarray(["gene1", "gene4"]),
        data={
            "va": np.asarray([80, 90]),
            "vb": np.asarray([True, False]),
            "vc": np.asarray(["h", "i"]),
        },
    )
    adata2 = ad.AnnData(X=X2, obs=obs2, var=var2)

    # ----------------------------------------------------------------
    # Write them into SOMAs
    soma1 = tiledbsc.SOMA((tmp_path / "soma1").as_posix())
    soma2 = tiledbsc.SOMA((tmp_path / "soma2").as_posix())

    tiledbsc.io.from_anndata(soma1, adata1)
    tiledbsc.io.from_anndata(soma2, adata2)

    # ----------------------------------------------------------------
    # Slice query
    somas = [soma1, soma2]
    (obs_attrs, var_attrs) = tiledbsc.SOMA.find_common_obs_and_var_keys(somas)
    soma_slices = tiledbsc.SOMA.queries(somas, obs_attrs=obs_attrs, var_attrs=var_attrs)

    assert len(soma_slices) == 2

    result_soma_slice = tiledbsc.SOMASlice.concat(soma_slices)
    assert result_soma_slice is not None

    adatac = result_soma_slice.to_anndata()

    # ----------------------------------------------------------------
    # Store the slice to SOMA
    somac = tiledbsc.SOMA((tmp_path / "somac").as_posix())
    tiledbsc.io.from_anndata(somac, adatac)

    # ----------------------------------------------------------------
    # The primary success of this test is that it runs to the end without an exception.
    # We would get exceptions like
    #   NotImplementedError: boolean inferred dtype not supported
    # if a bool column with NaNs were encountered first, or
    #   tiledb.cc.TileDBError: Failed to convert buffer for attribute: 'vc'
    # if a string column with NaNs were encountered first.

    assert somac.obs.shape() == (4, 2)
    assert somac.var.shape() == (4, 3)
    assert somac.X.data.shape() == (4, 4)
