# `tiledbsoma._dask`
Utilities for exporting [`Experiment`]s and [`ExperimentAxisQuery`]s to [anndata], with [Dask Array] `X` matrices.

See also: [Tutorial: Scanpy + Dask-AnnData][tutorial]

Passing a `tiledbsoma._dask.util.DaskConfig` to `Experiment.to_anndata` or `ExperimentAxisQuery.to_anndata` will cause the resulting [anndata]'s `X` matrix to be a [Dask Array].

It's recommended to use one worker thread per Dask process; performance (wall-clock and memory use) suffers when more than one Dask worker thread runs in the same Dask worker process:

```python
from os import cpu_count
from dask.distributed import Client, LocalCluster

cluster = LocalCluster(
    n_workers=min(16, cpu_count()),
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
