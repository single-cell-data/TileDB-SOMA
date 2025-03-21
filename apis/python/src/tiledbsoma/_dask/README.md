# `tiledbsoma._dask`
Utilities for exporting [`Experiment`]s and [`ExperimentAxisQuery`]s to [anndata], with [Dask Array] `X` matrices.

See also: [Tutorial: Scanpy + Dask-AnnData][tutorial]

Passing a `tiledbsoma._dask.util.DaskConfig` to `Experiment.to_anndata` or `ExperimentAxisQuery.to_anndata` will cause the resulting [anndata]'s `X` matrix to be a [Dask Array].

## Memory management

### Disable `distributed.worker.memory.{spill,pause,terminate}`
Dask spilling and worker pause/terminate logic mostly interferes with our usage. It doesn't take into account available swap space and, once workers begin to spill
or pause, HVG/PCA jobs below never recover.

```python
import dask

dask.config.set(
    {
        # "distributed.worker.memory.target": 0.7,
        "distributed.worker.memory.spill": False,
        "distributed.worker.memory.pause": False,
        "distributed.worker.memory.terminate": False,
    }
)
```

### `threads_per_worker=1`, `n_workers` based on available memory
It's recommended to use one worker thread per Dask process; performance (wall-clock and memory use) suffers when more than one Dask worker thread runs in the same Dask worker process.

In addition, each Dask worker process will use an amount of memory
correlated with the number of rows per chunk (`DaskConfig.chunk_size`). 10 GiB per worker process works well with 100k-row chunks:

```python
from os import cpu_count
from dask.distributed import Client, LocalCluster
from psutil import virtual_memory
available_mem = virtual_memory().available
worker_mem_allowance = 10 * 2**30  # 10GiB, good estimate for 100k-row Dask chunks

cluster = LocalCluster(
    n_workers=max(2, min(
        available_mem // worker_mem_allowance,
        cpu_count(),
    )),
    threads_per_worker=1,  # Performance is better with just one thread per worker process
)
client = Client(cluster)
```

[`Experiment`]: https://tiledbsoma.readthedocs.io/en/stable/python-tiledbsoma-experiment.html
[`ExperimentAxisQuery`]: https://tiledbsoma.readthedocs.io/en/stable/python-tiledbsoma-experimentaxisquery.html
[anndata]: https://anndata.readthedocs.io/en/stable/
[Dask Array]: https://docs.dask.org/en/stable/array.html
[Scanpy]: https://scanpy.readthedocs.io/en/stable/
[scanpy dask nb]: https://scanpy.readthedocs.io/en/stable/tutorials/experimental/dask.html
[tutorial]: ../../../notebooks/tutorial_scanpy_pca_dask.ipynb
