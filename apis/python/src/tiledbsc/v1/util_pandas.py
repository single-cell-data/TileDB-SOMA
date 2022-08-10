import pandas as pd


def ascii_to_unicode_pandas_readback(df: pd.DataFrame) -> pd.DataFrame:
    """
    Implements the 'decode on read' part of our ASCII/Unicode logic
    """
    # TODO: COMMENT/LINK HEAVILY
    for k in df:
        dfk = df[k]
        if len(dfk) > 0 and type(dfk.iat[0]) == bytes:
            df[k] = dfk.map(lambda e: e.decode())
    return df
