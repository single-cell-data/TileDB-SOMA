This tool is a profiler for [SOMA](https://github.com/single-cell-data/SOMA/tree/main) workloads collecting both performance and memory usage across runs and can also be used for detecting performance or memory hot spots.

Here's an example for how to profile a SOMA script (here `tests/objects.py`)

This example shows how to generate the profile data for a given SOMA run (here `tests/objects.py`)
```shell
python -m profiler python tests/objects.py
```

This will capture the following context information and metrics in a databased, stored under `profiling_runs`:

Context info:
* command
* timestamp
* stdout
* stderr
* tiledb_stats
* somacore_version
* tiledbsoma_version
* host_context
* custom_out

Metrics captured:
* user_time_sec
* system_time_sec
* pct_of_cpu
* elapsed_time_sec
* avg_shared_text_sz_kb
* avg_unshared_text_sz_kb
* avg_stack_sz_kb
* avg_total_sz_kb
* max_res_set_sz_kb
* avg_res_set_sz_kb
* major_page_faults
* minor_page_faults
* voluntary_context_switches
* involuntary_context_switches
* swaps
* file_system_inputs
* file_system_outputs
* socket_messages_sent
* socket_messages_received
* signals_delivered
* page_size_bytes
* exit_status


To report the metrics from multiple runs of the profiled script as JSON output:
```shell
python -m profiler.report -j "python tests/objects.py"
```

To report on a metric from multiple runs of the profiled script with a graphical plot:
```shell
python -m profiler.report  -m <metric_name> "python tests/objects.py"
```

The profiling data includes the following metrics:

 

To capture TileDB stats, you must instrument your command script as follows:

```code
from profiler import data

db = data.FileBasedProfileDB()
db.extract('python tests/objects.py', 'peak_memory')
```
If you want TileDB stats to also be added to the profile, please surround your code as follows:
```code
tiledb.stats_enable()
tiledb.stats_reset()

# your profiling code goes here
...

tiledb.stats_disable()

import json
with open("tiledb_stats.json", "w") as f:
    json.dump(tiledb.stats_dump(json=True), f)

```
The profiler automatically reads and stores the content of `tiledb_stats.txt` file.



