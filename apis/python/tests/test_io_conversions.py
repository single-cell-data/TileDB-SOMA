import numpy as np
import pandas as pd
import pyarrow as pa
import pytest

from tiledbsoma.io.conversions import df_to_arrow_table


@pytest.mark.parametrize(
    ("input_df", "expected"),
    [
        [
            pd.DataFrame(
                data={
                    "id": np.array([1, 3, 5, 6], dtype=np.int64),
                    "alpha": np.arange(4, dtype=np.float32),
                },
                index=[1, 3, 5, 6],
            ),
            pa.Table.from_pydict(
                {
                    "id": pa.array([1, 3, 5, 6], type=pa.int64()),
                    "alpha": pa.array([0, 1, 2, 3], type=pa.float32()),
                }
            ),
        ],
    ],
)
def test_df_to_arrow_table(input_df: pd.DataFrame, expected: pa.Table):
    actual = df_to_arrow_table(input_df)
    assert actual == expected
