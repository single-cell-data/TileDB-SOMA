#!/usr/bin/env python

import libtiledbsc
import tiledbsc.util
import pandas as pd

print(libtiledbsc.version())

uri = "s3://tiledb-singlecell-data/soco/soco3/Krasnow"
long_test = False
verbose = False
write_csv = False

# Setup config
config = {}
config["soma.init_buffer_bytes"] = f"{4 * 1024**3}"

if verbose:
    config["config.logging_level"] = "5"
    libtiledbsc.config_logging("debug")

# Open SOMA
s = tiledbsc.util.get_start_stamp()
sq = libtiledbsc.SOMA(uri, config).query()
print(tiledbsc.util.format_elapsed(s, f"Open SOMA: URI = {uri}"))

# Query obs, var, and sliced X/data
s = tiledbsc.util.get_start_stamp()
cell_type = "alveolar macrophage" if long_test else "pulmonary ionocyte"
sq.set_obs_condition("cell_type", cell_type, libtiledbsc.Condition.EQ)
sq.select_obs_attrs(["cell_type"])
sq.select_var_attrs(["var_id"])

dfs = []
while chunk := sq.next_results():
    dfs.append(chunk["soma/X/data"].to_pandas())
df = pd.concat(dfs)
print(tiledbsc.util.format_elapsed(s, f"X/data query"))

# Print shape of result
print(f"X/data: {df.shape}")

# Save results to csv
if write_csv:
    # reorder to match python
    df = df[["obs_id", "var_id", "value"]]
    df.to_csv(f"lib.csv", index=False)
