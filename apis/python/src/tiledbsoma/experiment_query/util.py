import pandas as pd
import pyarrow as pa


def X_as_series(tbl: pa.Table) -> pd.Series:
    """
    Convert SOMA 2D data from Arrow Table to Pandas Series [lifecycle: experimental].

    NOTE: this operation is not zero copy.
    """
    data = tbl["soma_data"].to_numpy()
    dim_0 = tbl["soma_dim_0"].to_numpy()
    dim_1 = tbl["soma_dim_1"].to_numpy()
    return pd.Series(
        data,
        pd.MultiIndex.from_arrays((dim_0, dim_1), names=("soma_dim_0", "soma_dim_1")),
        dtype=pd.SparseDtype(data.dtype, fill_value=0),
        name="soma_data",
    )
