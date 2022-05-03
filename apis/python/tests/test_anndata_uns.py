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
            "list of strings": list(str(i) for i in range(10)),
            "simple dict": {"A": 0, "B": "one"},
        },
    )
    SOMA(tmp_path.as_posix()).from_anndata(adata)

    # TODO: persistance for UNS is not yet defined. Update when it is.
    assert (tmp_path / "uns").exists()

    # TODO: add check that contents were saved correctly
