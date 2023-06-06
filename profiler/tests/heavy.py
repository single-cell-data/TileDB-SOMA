import sys
from io import StringIO

import cellxgene_census
import tiledb

""" This is heavy read and process SOMA script sample developed for early testing of the profiler."""
print("Activating TileDB stats\n")

tiledb.stats_enable()
tiledb.stats_reset()

census = cellxgene_census.open_soma(census_version="2023-05-15")
exp = census["census_data"]["mus_musculus"]
# For humans, we can use this instead
# exp = census["census_data"]["homo_sapiens"]
# Grab the experiment containing human data, and the measurement therein with RNA
q = exp.axis_query(measurement_name="RNA")
print("Reading X\n")
itr = q.X("raw")

# '''
print("\n Iterations")
t = next(itr.tables())
it = 0
while t:
    print(f"Iteration {it} {t}")
    try:
        t = next(itr.tables())
    except StopIteration:
        break
tmp = sys.stdout
my_result = StringIO()
sys.stdout = my_result
tiledb.stats_dump()
tiledb.stats_disable()
sys.stdout = tmp

with open("tiledb_stats.txt", "w") as f:
    f.write(my_result.getvalue())
