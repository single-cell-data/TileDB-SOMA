This tool is profiler for [SOMA](https://github.com/single-cell-data/SOMA/tree/main) workloads collecting both performance and memory usage across runs and alos can be used for detecting performance or memory hot spots.

Here's an example for how to profile a SOMA script (here ```tests/objects.py```)

This example shows how to generate the profile data for a given SOMA run (here ```tests/objects.py```)
```shell
python main.py python tests/objects.py
```

The profiled data includes the following items (for the time being):

* process: The process and its parameters to be profile 
* custom_out: list of custom profilers to be stored
* date
* time

Time (perf) stats:

* rt: Real time
* ut: User time
* st: System time

Memory stats:

* max_set_size
* page_reclaims
* page_faults
* cycles_elapsed
* peak_memory
* tiledb_stats
* somacore_version
* tiledbsoma_version

Context data:

* uname: uname -a
* total_virtual_mem
* total_physical_mem
* swap_mem 
* cpu_count
* python_version

We also store tileDB stats.

This example shows how to use a profiler plot generator on the same generated profile:
```shell
python profiler_plot.py python tests/objects.py -m ut
```
Finally, this example shows how to use the profiler inside your code (assuming your file is in the same directory as the profiler):
```code
from profiler import data

db = data.FileBasedProfileDB()
db.extract('python tests/objects.py', 'peak_memory')
```
If you want TileDB stats to also be added to the profile, please surround your code as follows:
```code
tiledb.stats_enable()
tiledb.stats_reset()

...

tmp = sys.stdout
my_result = StringIO()
sys.stdout = my_result
tiledb.stats_dump()
tiledb.stats_disable()
sys.stdout = tmp

with open("tiledb_stats.txt", "w") as f:
    f.write(my_result.getvalue())

```
The profiler automatically reads and stores the content of ```tiledb_stats.txt``` file.