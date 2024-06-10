import numpy as np
import pandas as pd
import pytest

import tiledbsoma as soma

from .common import TestReadPythonWriteR


class TestDataframeWriteRReadPython(TestReadPythonWriteR):
    @pytest.fixture(scope="class")
    def R_dataframe(self):
        base_script = f"""
        library(tiledbsoma)
        library(arrow)

        df_schema <- schema(
            field("foo", int32()),
            field("bar", float64()),
            field("baz", string()),
            field("quux", bool())
        )

        sdf <- SOMADataFrameCreate("{self.uri}", df_schema, "foo")

        df <- data.frame(
            soma_joinid = bit64::as.integer64(c(1,2,3,4,5)),
            foo = as.integer(c(10, 20, 30, 40, 50)),
            bar = c(4.1, 5.2, 6.3, 7.4, 8.5),
            baz = c("apple", "ball", "cat", "dog", "egg"),
            quux = c(TRUE, FALSE, FALSE, TRUE, FALSE)
        )
        tbl <- arrow_table(df)
        sdf$write(tbl)
        """
        self.execute_R_script(base_script)

    def test_dataframe_length_matches(self, R_dataframe):
        with soma.open(self.uri) as sdf:
            df = sdf.read().concat().to_pandas()
            assert len(df) == 5

    def test_dataframe_columns_match(self, R_dataframe):
        with soma.open(self.uri) as sdf:
            df = sdf.read().concat().to_pandas()
            assert df["foo"].equals(pd.Series([10, 20, 30, 40, 50], dtype=np.int32))
            assert df["bar"].equals(
                pd.Series([4.1, 5.2, 6.3, 7.4, 8.5], dtype=np.float64)
            )
            assert df["baz"].equals(pd.Series(["apple", "ball", "cat", "dog", "egg"]))
            assert df["quux"].equals(
                pd.Series([True, False, False, True, False], dtype=bool)
            )
