import pandas as pd
import scipy.sparse as sp


def csr_from_tiledb_df(df: pd.DataFrame, num_rows: int, num_cols: int) -> sp.csr_matrix:
    """
    Given a tiledb dataframe, return a ``scipy.sparse.csr_matrx``.
    """
    return sp.csr_matrix(
        (df["soma_data"], (df["soma_dim_0"], df["soma_dim_1"])),
        shape=(num_rows, num_cols),
    )
