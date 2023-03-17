import shutil
import subprocess

import pyarrow as pa

import tiledbsoma as soma

# Remove artifacts from previous runs, if they exist
uri = "test-dataframe.auto.soma"
shutil.rmtree(uri)

asch = pa.schema(
    [
        ("foo", pa.int32()),
        ("bar", pa.float64()),
        ("baz", pa.large_string()),
        ("quux", pa.bool_()),
    ]
)

soma.DataFrame.create(uri, schema=asch, index_column_names=["foo"]).close()

pydict = {}
pydict["soma_joinid"] = [0, 1, 2, 3, 4]
pydict["foo"] = [10, 20, 30, 40, 50]
pydict["bar"] = [4.1, 5.2, 6.3, 7.4, 8.5]
pydict["baz"] = ["apple", "ball", "cat", "dog", "egg"]
pydict["quux"] = [True, False, False, True, False]
rb = pa.Table.from_pydict(pydict)
with soma.DataFrame.open(uri, "w") as sdf:
    sdf.write(rb)

df = rb.to_pandas()  # Alternatively, we could re-open the dataframe

# Prepares an R script with the dependencies and loads the data.frame in `df`.
# Future assertions can be done by appending to this script
R_script_base = f"""
require("tiledbsoma")
require("testthat")
platform_config <- NULL
tiledbsoma_ctx <- NULL
soma_df <- SOMADataFrame$new("{uri}", platform_config, tiledbsoma_ctx)
table = soma_df$read()
df = as.data.frame(table)
"""


def r_assert(code: str):
    R_script = R_script_base + code
    with open("test-dataframe-read.R", "w") as f:
        f.write(R_script)
    subprocess.run(["Rscript", "test-dataframe-read.R"])


r_assert(f"expect_length(df, {len(df)})")


def to_R(x):
    if isinstance(x, bool):
        return str(x).upper()
    elif isinstance(x, str):
        return f'"{x}"'
    else:
        return str(x)


def create_R_list(xs):
    return ",".join([to_R(x) for x in xs])


for key in df.keys():
    col = df[key].tolist()
    R_list = create_R_list(col)

    r_assert(f"""expect_equal(as.list(df)$"{key}", c({R_list}))""")
