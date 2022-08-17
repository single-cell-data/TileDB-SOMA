#!/usr/bin/env python

import tiledbsc
import tiledbsc.util
import tiledb
import pandas as pd

print(
    f"tiledbsc={tiledbsc.__version__} tiledb={tiledb.__version__} libtiledb={tiledb.libtiledb.version()}"
)

uri = "s3://tiledb-singlecell-data/soco/soco3/Krasnow"
long_test = False
verbose = False
write_csv = False

# Setup config
cfg = tiledb.Config()
cfg["py.init_buffer_bytes"] = 4 * 1024**3

if verbose:
    cfg["config.logging_level"] = 5

# Open SOMA
s = tiledbsc.util.get_start_stamp()
soma = tiledbsc.SOMA(uri, ctx=tiledb.Ctx(cfg))
print(tiledbsc.util.format_elapsed(s, f"Open SOMA: URI = {uri}"))

# Query obs, var, and sliced X/data
s = tiledbsc.util.get_start_stamp()
cell_type = "alveolar macrophage" if long_test else "pulmonary ionocyte"
slice = soma.query(obs_query_string=f'cell_type=="{cell_type}"')
print(tiledbsc.util.format_elapsed(s, f"X/data query"))

# Print shape of result
print(slice)

# Save results to csv
if write_csv:
    df = slice.X["data"].reset_index()
    df.to_csv(f"py.csv", index=False)
