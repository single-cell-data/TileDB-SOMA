import pandas as pd

from . import util_tiledb


def _decategoricalize_obs_or_var(obs_or_var: pd.DataFrame) -> pd.DataFrame:
    """
    Performs a typecast into types that TileDB can persist.
    """
    if len(obs_or_var.columns) > 0:
        return pd.DataFrame.from_dict(
            {
                k: util_tiledb.to_tiledb_supported_array_type(v)
                for k, v in obs_or_var.items()
            },
        )
    else:
        return obs_or_var
