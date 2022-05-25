import anndata
import tiledb
import tiledbsc

import pytest
import tempfile
import os
from pathlib import Path
import math

HERE = Path(__file__).parent


def relerr(actual, expect):
    return math.fabs(actual - expect) / expect


def relerr_ok(actual, expect):
    return relerr(actual, expect) < 1e-6


def test_chunked_writes(tmp_path):
    """
    This ingests datasets having X stored in dense, CSR, and CSC formats, and
    verifies the chunked-ingest logic is the same.
    """

    # We have three copies of pbmc3k_processed.h5ad checked into source control.
    # See the README.md file there for how they were prepared -- essentially
    #   ann = anndata.read_h5ad('pbmc3k_processed.h5ad')
    #   ann.X = scipy.sparse.csc_matrix(ann.X)
    #   ann.write_h5ad('pbmc3k-x-csc.h5ad')
    # along with nulling out obsm/uns/etc which are not used in this test (to make the files a
    # little smaller).
    ann_dense_path = HERE.parent / "anndata/pbmc3k-x-dense.h5ad"
    ann_csr_path = HERE.parent / "anndata/pbmc3k-x-csr.h5ad"
    ann_csc_path = HERE.parent / "anndata/pbmc3k-x-csc.h5ad"

    soma_dense_path = (tmp_path / "soma_dense").as_posix()
    soma_csr_path = (tmp_path / "soma_csr").as_posix()
    soma_csc_path = (tmp_path / "soma_csc").as_posix()

    # This value is specifically chosen to ensure that the curated datafiles
    # on this test will get their X ingested in multiple chunks -- this is
    # important in order to test the chunked-ingest logic for all three variants:
    # dense, CSR, and CSC.
    sopt = tiledbsc.SOMAOptions(goal_chunk_nnz=1000000)

    # Do the ingests from anndata .h5ad files to TileDB SOMAs.
    soma_dense = tiledbsc.SOMA(soma_dense_path, soma_options=sopt)
    soma_csr = tiledbsc.SOMA(soma_csr_path, soma_options=sopt)
    soma_csc = tiledbsc.SOMA(soma_csc_path, soma_options=sopt)

    soma_dense.from_h5ad(ann_dense_path)
    soma_csr.from_h5ad(ann_csr_path)
    soma_csc.from_h5ad(ann_csc_path)

    # Read the X arrays into memory as pandas dataframe objects.
    xdf_dense = soma_dense.X.data.df()
    xdf_csr = soma_csr.X.data.df()
    xdf_csc = soma_csc.X.data.df()

    # We want to verify that the X data that made it through the chunking algorithms
    # is the same.

    # A first check is the shape of the dataframes. Note that once densified
    # these are (2638, 1838) but as dataframes they're (2638*1838, 3): there
    # is an 'obs_id' column, a 'var_id' column, and a 'value' column.
    assert xdf_dense.shape == (4848644, 3)
    assert xdf_csr.shape == (4848644, 3)
    assert xdf_csc.shape == (4848644, 3)

    # Compare indices
    assert list(xdf_dense["obs_id"]) == list(xdf_csr["obs_id"])
    assert list(xdf_dense["obs_id"]) == list(xdf_csc["obs_id"])
    assert list(xdf_dense["var_id"]) == list(xdf_csr["var_id"])
    assert list(xdf_dense["var_id"]) == list(xdf_csc["var_id"])

    # This takes shockingly long -- for this reason, we instead compare
    # sum-of-abs, and sums over group-bys, which are quicker.
    # for i in range(4848644):
    #    assert(relerr_ok(xdf_dense['value'][i], xdf_csr['value'][i]))
    #    assert(relerr_ok(xdf_dense['value'][i], xdf_csc['value'][i]))

    # Next compute the sum of absolute values over the X arrays. These
    # won't be 100% identical since floating-point arithmetic isn't perfectly
    # commutative, but we can compare them within a tolerance.
    xdf_dense_abs_sum = xdf_dense["value"].abs().sum()
    xdf_csr_abs_sum = xdf_csr["value"].abs().sum()
    xdf_csc_abs_sum = xdf_csc["value"].abs().sum()

    assert relerr_ok(xdf_dense_abs_sum, 2042387.8)
    assert relerr_ok(xdf_csr_abs_sum, 2042387.8)
    assert relerr_ok(xdf_csc_abs_sum, 2042387.8)

    # Next we can compute sums grouped by obs_id, check shapes, and do some point-checks.
    #
    #                       value                        value                        value
    # obs_id                       obs_id                       obs_id
    # AAACATACAACCAC-1 -54.403599  AAACATACAACCAC-1 -54.403599  AAACATACAACCAC-1 -54.403599
    # AAACATTGAGCTAC-1 -74.252998  AAACATTGAGCTAC-1 -74.252998  AAACATTGAGCTAC-1 -74.252998
    # AAACATTGATCAGC-1  14.642952  AAACATTGATCAGC-1  14.642952  AAACATTGATCAGC-1  14.642952
    # AAACCGTGCTTCCG-1  25.272526  AAACCGTGCTTCCG-1  25.272526  AAACCGTGCTTCCG-1  25.272526
    # AAACCGTGTATGCG-1   6.548214  AAACCGTGTATGCG-1   6.548214  AAACCGTGTATGCG-1   6.548214
    # ...                     ...  ...                     ...  ...                     ...
    # TTTCGAACTCTCAT-1  90.523659  TTTCGAACTCTCAT-1  90.523659  TTTCGAACTCTCAT-1  90.523659
    # TTTCTACTGAGGCA-1 -11.782106  TTTCTACTGAGGCA-1 -11.782106  TTTCTACTGAGGCA-1 -11.782106
    # TTTCTACTTCCTCG-1 -64.215446  TTTCTACTTCCTCG-1 -64.215446  TTTCTACTTCCTCG-1 -64.215446
    # TTTGCATGAGAGGC-1  39.785194  TTTGCATGAGAGGC-1  39.785194  TTTGCATGAGAGGC-1  39.785194
    # TTTGCATGCCTCAC-1 -83.011658  TTTGCATGCCTCAC-1 -83.011658  TTTGCATGCCTCAC-1 -83.011658
    # [2638 rows x 1 columns]      [2638 rows x 1 columns]      [2638 rows x 1 columns]
    xdf_dense_obs_sums = xdf_dense.groupby("obs_id").sum()
    xdf_csr_obs_sums = xdf_csr.groupby("obs_id").sum()
    xdf_csc_obs_sums = xdf_csc.groupby("obs_id").sum()

    assert xdf_dense_obs_sums.shape == (2638, 1)
    assert xdf_csr_obs_sums.shape == (2638, 1)
    assert xdf_csc_obs_sums.shape == (2638, 1)

    assert relerr_ok(xdf_dense_obs_sums["value"]["TTTCGAACTCTCAT-1"], 90.523659)
    assert relerr_ok(xdf_csr_obs_sums["value"]["TTTCGAACTCTCAT-1"], 90.523659)
    assert relerr_ok(xdf_csc_obs_sums["value"]["TTTCGAACTCTCAT-1"], 90.523659)

    # Similarly, we can compute sums grouped by var_id, check shapes, and do some point-checks.
    #
    #                value                    value                    value
    # var_id                   var_id                   var_id
    # AAGAB  -6.655109e+00     AAGAB  -6.655109e+00     AAGAB  -6.655109e+00
    # AAR2   -3.985365e+00     AAR2   -3.985365e+00     AAR2   -3.985365e+00
    # AATF    8.046627e-07     AATF    8.046627e-07     AATF    8.046627e-07
    # ABCB1  -1.818954e+01     ABCB1  -1.818954e+01     ABCB1  -1.818954e+01
    # ABCC10 -1.747384e+01     ABCC10 -1.747384e+01     ABCC10 -1.747384e+01
    # ...              ...     ...              ...     ...              ...
    # ZRANB3 -4.315807e+01     ZRANB3 -4.315807e+01     ZRANB3 -4.315807e+01
    # ZSWIM6 -3.478199e+01     ZSWIM6 -3.478199e+01     ZSWIM6 -3.478199e+01
    # ZUFSP  -7.315041e+00     ZUFSP  -7.315041e+00     ZUFSP  -7.315041e+00
    # ZWINT  -4.822803e+01     ZWINT  -4.822803e+01     ZWINT  -4.822803e+01
    # ZYX    -6.118789e-07     ZYX    -6.118789e-07     ZYX    -6.118789e-07
    # [1838 rows x 1 columns]  [1838 rows x 1 columns]  [1838 rows x 1 columns]
    xdf_dense_var_sums = xdf_dense.groupby("var_id").sum()
    xdf_csr_var_sums = xdf_csr.groupby("var_id").sum()
    xdf_csc_var_sums = xdf_csc.groupby("var_id").sum()

    assert xdf_dense_var_sums.shape == (1838, 1)
    assert xdf_csr_var_sums.shape == (1838, 1)
    assert xdf_csc_var_sums.shape == (1838, 1)

    assert relerr_ok(xdf_dense_var_sums["value"]["ZRANB3"], -4.315807e01)
    assert relerr_ok(xdf_csr_var_sums["value"]["ZRANB3"], -4.315807e01)
    assert relerr_ok(xdf_csc_var_sums["value"]["ZRANB3"], -4.315807e01)
