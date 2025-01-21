# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from typing import cast

import pandas as pd
import pyarrow as pa

from ._types import NPNDArray, PDSeries


def X_as_series(tbl: pa.Table) -> PDSeries:
    """Convert the COO 2D matrix data returned from ``SparseNDArray.read_table()``
    into a Pandas Series, with coordindates as a Pandas MultiIndex.

    Note:
        This operation is not zero-copy.

    Args:
        tbl:
            A Table containing COO-formated sparse matrix data.

    Returns:
        pandas.Series:
            COO data, with coordinates in a MultiIndex.

    Example:
        >>> tbl = pa.Table.from_arrays(
        ...     [ np.array([0, 2]), np.array([1, 3]), np.array([1.1, 4.2], dtype=np.float32) ],
        ...     names=['soma_dim_0', 'soma_dim_1', 'soma_data']
        ... )
        >>> X_as_series(tbl)
        soma_dim_0  soma_dim_1
        0           1             1.1
        2           3             4.2
        Name: soma_data, dtype: Sparse[float32, 0]

    Lifecycle:
        Maturing.
    """

    data = tbl["soma_data"].to_numpy()
    dim_0 = tbl["soma_dim_0"].to_numpy()
    dim_1 = tbl["soma_dim_1"].to_numpy()
    return pd.Series(
        cast(NPNDArray, data),
        pd.MultiIndex.from_arrays((dim_0, dim_1), names=("soma_dim_0", "soma_dim_1")),
        dtype=pd.SparseDtype(data.dtype, fill_value=0),
        name="soma_data",
    )
