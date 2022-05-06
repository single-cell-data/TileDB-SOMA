from anndata import AnnData
import tiledb
from tiledbsc import SOMA
import pandas as pd
import numpy as np
from scipy import sparse

import pytest

"""
Verify that the AnnData.uns persists correctly. Currently focused on a simple
subset used by cellxgene schema.

TODO: round-trip via to_anndata() when that is implemented.
TODO: enable test when `uns` supported (issue #54)

"""

# TODO: re-enable when #54 is resolved
@pytest.mark.skip(reason="Unimplemented: filed as issue #54")
def test_from_anndata_uns(tmp_path):
    """
    Test very simple `uns` is persisted.
    """

    adata = AnnData(
        X=np.ones(10, 10, dtype=np.float32),
        uns={
            "number": 1,
            "string": "a string",
            "list of string": list(str(i) for i in range(10)),
            "list_of_int": list(i*10 for i in range(10)),
            "list_of_float": list(i*1.25 for i in range(10)),
            "simple dict": {"A": 0, "B": "one"},
            "numpy_ndarray_1d": np.asarray([1.,2.,3.]),
            "numpy_ndarray_2d": np.asarray([[1.,2.,3.],[4.,5.,6.]]),
            "pandas_dataframe": pd.DataFrame(
                index=np.arange(n_obs).astype(str), data={"A": np.arange(10, dtype=np.int32)}
            ),

        },
    )
    SOMA(tmp_path.as_posix()).from_anndata(adata)

    # TODO: persistance for UNS is not yet defined. Update when it is.
    assert (tmp_path / "uns").exists()

    # TODO: add check that contents were saved correctly
