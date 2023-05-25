import pyarrow as pa
import pytest

import tiledbsoma as soma

from .common import TestWritePythonReadR, embed_python_list_into_R_code


class TestDataframeWritePythonReadR(TestWritePythonReadR):
    @pytest.fixture(scope="class")
    def dataframe(self):
        asch = pa.schema(
            [
                ("foo", pa.int32()),
                ("bar", pa.float64()),
                ("baz", pa.large_string()),
                ("quux", pa.bool_()),
            ]
        )

        soma.DataFrame.create(self.uri, schema=asch, index_column_names=["foo"]).close()

        pydict = {}
        pydict["soma_joinid"] = [0, 1, 2, 3, 4]
        pydict["foo"] = [10, 20, 30, 40, 50]
        pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
        pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
        pydict["quux"] = [True, False, False, True, False]
        rb = pa.Table.from_pydict(pydict)
        with soma.DataFrame.open(self.uri, "w") as sdf:
            sdf.write(rb)

        df = rb.to_pandas()  # Alternatively, we could re-open the dataframe
        return df

    # Prepares an R script with the dependencies and loads the data.frame in `df`.
    # Future assertions can be done by appending to this script
    def base_R_script(self):
        return f"""
        library("tiledbsoma")
        soma_df <- SOMADataFrameOpen("{self.uri}")
        table = soma_df$read()$concat()
        df = as.data.frame(table)
        """

    # tests
    def test_dataframe_length_matches(self, dataframe):
        self.r_assert(f"stopifnot(length(df) == {len(dataframe)})")

    def test_dataframe_columns_match(self, dataframe):
        for key in dataframe.keys():
            col = dataframe[key].tolist()
            R_list = embed_python_list_into_R_code(col)
            self.r_assert(f"""stopifnot(all.equal(as.list(df)$"{key}", c({R_list})))""")
