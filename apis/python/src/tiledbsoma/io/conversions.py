from typing import TypeVar, cast

import numpy as np
import pandas as pd
import pandas._typing as pdt
import scipy.sparse as sp
from pandas.api.types import infer_dtype, is_categorical_dtype

from ..funcs import typeguard_ignore
from ..types import NPNDArray, PDSeries

_DT = TypeVar("_DT", bound=pdt.Dtype)
_MT = TypeVar("_MT", NPNDArray, sp.spmatrix, PDSeries)
_str_to_type = {"boolean": bool, "string": str, "bytes": bytes}


def decategoricalize_obs_or_var(obs_or_var: pd.DataFrame) -> pd.DataFrame:
    """
    Performs a typecast into types that TileDB can persist.
    """
    if len(obs_or_var.columns) > 0:
        return pd.DataFrame.from_dict(
            {k: to_tiledb_supported_array_type(v) for k, v in obs_or_var.items()},
        )
    else:
        return obs_or_var


@typeguard_ignore
def _to_tiledb_supported_dtype(dtype: _DT) -> _DT:
    """A handful of types are cast into the TileDB type system."""
    # TileDB has no float16 -- cast up to float32
    return cast(_DT, np.dtype("float32")) if dtype == np.dtype("float16") else dtype


def to_tiledb_supported_array_type(x: _MT) -> _MT:
    """
    Converts datatypes unrepresentable by TileDB into datatypes it can represent.
    E.g., categorical strings -> string.

    See also [https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html](https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html).

    Preferentially converts to the underlying primitive type, as TileDB does not support
    most complex types. NOTE: this does not support ``datetime64`` conversion.

    Categoricals are a special case. If the underlying categorical type is a primitive,
    convert to that. If the array contains NA/NaN (i.e. not in the category, code == -1),
    raise error unless it is a float or string.
    """
    if isinstance(x, (np.ndarray, sp.spmatrix)) or not is_categorical_dtype(x):
        # mypy issues a spurious error here, but only when
        # _to_tiledb_supported_dtype is decorated with @typeguard_ignore???
        target_dtype = _to_tiledb_supported_dtype(x.dtype)  # type: ignore[arg-type]
        return x if target_dtype == x.dtype else x.astype(target_dtype)

    categories = x.cat.categories
    cat_dtype = categories.dtype
    if cat_dtype.kind in ("f", "u", "i"):
        if x.hasnans and cat_dtype.kind == "i":
            raise ValueError(
                "Categorical array contains NaN -- unable to convert to TileDB array."
            )
        # More mysterious spurious mypy errors.
        target_dtype = _to_tiledb_supported_dtype(cat_dtype)  # type: ignore[arg-type]
    else:
        # Into the weirdness. See if Pandas can help with edge cases.
        inferred = infer_dtype(categories)
        if x.hasnans and inferred in ("boolean", "bytes"):
            raise ValueError(
                "Categorical array contains NaN -- unable to convert to TileDB array."
            )
        target_dtype = np.dtype(  # type: ignore[assignment]
            _str_to_type.get(inferred, object)
        )

    return x.astype(target_dtype)


def csr_from_tiledb_df(df: pd.DataFrame, num_rows: int, num_cols: int) -> sp.csr_matrix:
    """
    Given a tiledb dataframe, return a ``scipy.sparse.csr_matrx``.
    """
    return sp.csr_matrix(
        (df["soma_data"], (df["soma_dim_0"], df["soma_dim_1"])),
        shape=(num_rows, num_cols),
    )
