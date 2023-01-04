import asyncio
import concurrent.futures
import contextvars
import functools
from typing import (
    TYPE_CHECKING,
    Any,
    AsyncIterator,
    Callable,
    ContextManager,
    Dict,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
    TypeVar,
    Union,
    cast,
)

import anndata
import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
from typing_extensions import Literal, ParamSpec, TypedDict

from ..dataframe import DataFrame as SOMADataFrame

if TYPE_CHECKING:
    from ..experiment import Experiment
from ..sparse_nd_array import SparseNDArray as SOMASparseNDArray
from .anndata import make_anndata
from .axis import AxisQuery, MatrixAxisQuery
from .eq_types import AxisColumnNames, ExperimentQueryReadArrowResult

AxisJoinIds = TypedDict(
    "AxisJoinIds",
    {
        "obs": pa.Array,
        "var": pa.Array,
    },
)


class ExperimentQuery(ContextManager["ExperimentQuery"]):
    """
    ExperimentQuery allows easy selection and extraction of data from a single soma.Measurement
    in a soma.Experiment [lifecycle: experimental].

    IMPORTANT: this class is not thread safe.

    IMPORTANT: this query class assumes it can store the full result of both axis dataframe
    queries in memory, and only provides incremental access to the underlying X NDArray. API
    features such as `n_obs` and `n_vars` codify this in the API.

    IMPORTANT: you must call `close()` on any instance of this class in order to release
    underlying resources. The ExperimentQuery is a context manager -- it is recommended
    that you use the following pattern to make this easy and safe:
    ```
        with ExperimentQuery(...) as query:
            ...
    ```
    """

    experiment: "Experiment"
    ms: str

    _query: MatrixAxisQuery
    _joinids: AxisJoinIds
    _indexer: "AxisIndexer"
    _default_threadpool: Optional[
        concurrent.futures.ThreadPoolExecutor
    ]  # lazily created, use property accessor

    def __init__(
        self,
        experiment: "Experiment",
        measurement_name: str,
        *,
        obs_query: Optional[AxisQuery] = None,
        var_query: Optional[AxisQuery] = None,
    ):
        if not experiment.exists():
            raise ValueError("Experiment does not exist")
        if measurement_name not in experiment.ms:
            raise ValueError("Measurement does not exist in the experiment")

        self.experiment = experiment
        self.ms = measurement_name

        self._query = {
            "obs": obs_query if obs_query is not None else AxisQuery(),
            "var": var_query if var_query is not None else AxisQuery(),
        }
        self._joinids = {
            "obs": None,
            "var": None,
        }
        self._indexer = AxisIndexer(self)

        self._default_threadpool = None

    @property
    def default_threadpool(self) -> concurrent.futures.ThreadPoolExecutor:
        if self._default_threadpool is None:
            # TODO: the user should be able to set their own threadpool, a la asyncio's
            # loop.set_default_executor().  This is important for managing the level of
            # concurrency, etc.
            self._default_threadpool = concurrent.futures.ThreadPoolExecutor()
        return self._default_threadpool

    def close(self) -> None:
        """
        Cleanup and close all resources. This must be called or the thread pool
        will not be release, etc.
        """
        if self._default_threadpool is not None:
            self._default_threadpool.shutdown()
            self._default_threadpool = None

    def __enter__(self) -> "ExperimentQuery":
        return self

    def __exit__(self, *excinfo: Any) -> None:
        self.close()

    def _read_axis_dataframe(
        self,
        axis: Literal["obs", "var"],
        axis_df: SOMADataFrame,
        *,
        column_names: Optional[Sequence[str]],
    ) -> pa.Table:
        """
        Read the specified axis. Will load and save the resulting soma_joinids for the
        axis, if they are not already known.
        """
        query = self._query[axis]
        need_joinids = self._joinids[axis] is None

        query_columns = column_names
        if (
            need_joinids
            and column_names is not None
            and "soma_joinid" not in column_names
        ):
            query_columns = ["soma_joinid"] + list(column_names)

        tbl = axis_df.read_all(
            ids=query.coords,
            value_filter=query.value_filter,
            column_names=query_columns,
        )

        if need_joinids:
            self._joinids[axis] = tbl.column("soma_joinid").combine_chunks()
        assert self._joinids[axis] is not None

        if column_names is not None:
            tbl = tbl.select(column_names)
        return tbl

    def _read_axis_joinids(
        self, axis: Literal["obs", "var"], axis_df: SOMADataFrame
    ) -> pa.Array:
        if self._joinids[axis] is None:
            self._read_axis_dataframe(axis, axis_df, column_names=["soma_joinid"])
        return self._joinids[axis]

    def obs(self, *, column_names: Optional[Sequence[str]] = None) -> pa.Table:
        """Return obs as an Arrow table."""
        return self._read_axis_dataframe(
            "obs", self.experiment.obs, column_names=column_names
        )

    def var(self, *, column_names: Optional[Sequence[str]] = None) -> pa.Table:
        """Return var as an Arrow table."""
        return self._read_axis_dataframe(
            "var", self.experiment.ms[self.ms].var, column_names=column_names
        )

    def obs_joinids(self) -> pa.Array:
        return self._read_axis_joinids("obs", self.experiment.obs)

    def var_joinids(self) -> pa.Array:
        return self._read_axis_joinids("var", self.experiment.ms[self.ms].var)

    @property
    def n_obs(self) -> int:
        return len(self.obs_joinids())

    @property
    def n_vars(self) -> int:
        return len(self.var_joinids())

    def _ensure_joinids_loaded(self) -> None:
        """Private. Ensure joinids for both axis are in-memory."""
        futures = []
        if not self._joinids["obs"]:
            futures.append(self.default_threadpool.submit(self.obs_joinids))
        if not self._joinids["var"]:
            futures.append(self.default_threadpool.submit(self.var_joinids))
        if futures:
            concurrent.futures.wait(futures)
        assert self._joinids["obs"] is not None
        assert self._joinids["var"] is not None

    def _fetchX(
        self, X: SOMASparseNDArray, prefetch: bool = False
    ) -> Iterator[pa.Table]:
        assert self._joinids["obs"] is not None
        assert self._joinids["var"] is not None

        obs_joinids = self._joinids["obs"]
        var_joinids = self._joinids["var"]

        if len(obs_joinids) == 0 or len(var_joinids) == 0:
            yield pa.Table.from_pylist([], schema=X.schema)
            return

        if not prefetch:
            # yield for clarity
            yield from X.read_table((obs_joinids, var_joinids))

        else:
            # prefetch
            fn = wrap_generator(X.read_table((obs_joinids, var_joinids)))
            _prefetch_future = self.default_threadpool.submit(fn)
            while True:
                value, done = _prefetch_future.result()
                if done:
                    return
                assert value is not None
                _prefetch_future = self.default_threadpool.submit(fn)
                yield value

    def X(self, layer: str, prefetch: bool = False) -> Iterator[pa.Table]:
        """
        Return an X layer as an iterator of Arrow Tables.
        """
        if not (layer and layer in self.experiment.ms[self.ms].X):
            raise ValueError("Must specify X layer")

        X = self.experiment.ms[self.ms].X[layer]
        if X.soma_type != "SOMASparseNDArray":
            raise NotImplementedError("Dense array unsupported")
        assert isinstance(X, SOMASparseNDArray)

        self._ensure_joinids_loaded()

        yield from self._fetchX(X, prefetch=prefetch)

    def read(
        self,
        X_name: str,
        *,
        use_position_indexing: bool = False,
        column_names: Optional[AxisColumnNames] = None,
        X_layers: Optional[List[str]] = None,
    ) -> ExperimentQueryReadArrowResult:
        """
        Read the _entire_ query result into Arrow Tables. Low-level routine
        intended to be the basis for exporting to other in-core formats, such
        as AnnData.
        """
        X_collection = self.experiment.ms[self.ms].X
        X_layers = [] if X_layers is None else X_layers
        all_X_names = [X_name] + X_layers
        all_X_arrays: dict[str, SOMASparseNDArray] = {}
        for _xname in all_X_names:
            if not isinstance(_xname, str) or not _xname:
                raise ValueError("X layer names must be specified as a string.")
            if _xname not in X_collection:
                raise ValueError("Unknown X layer name")
            # TODO: dense array slicing, which needs thought due to lack of point indexing
            if X_collection[_xname].soma_type != "SOMASparseNDArray" or not isinstance(
                X_collection[_xname], SOMASparseNDArray
            ):
                raise NotImplementedError("Dense array unsupported")
            all_X_arrays[_xname] = cast(SOMASparseNDArray, X_collection[_xname])

        if column_names is None:
            column_names = {"obs": None, "var": None}
        if "obs" not in column_names:
            column_names["obs"] = None
        if "var" not in column_names:
            column_names["var"] = None

        futures = (
            self.default_threadpool.submit(
                self._read_axis_dataframe,
                "obs",
                self.experiment.obs,
                column_names=column_names["obs"],
            ),
            self.default_threadpool.submit(
                self._read_axis_dataframe,
                "var",
                self.experiment.ms[self.ms].var,
                column_names=column_names["var"],
            ),
        )
        concurrent.futures.wait(futures)
        obs_table, var_table = (f.result() for f in futures)

        X_tables = {
            _xname: pa.concat_tables(self._fetchX(all_X_arrays[_xname], prefetch=True))
            for _xname in all_X_arrays
        }
        if use_position_indexing:
            X_tables = self._rewrite_X_for_positional_indexing(X_tables)

        X = X_tables.pop(X_name)
        query_result: ExperimentQueryReadArrowResult = {
            "obs": obs_table,
            "var": var_table,
            "X": X,
        }
        if len(X_layers) > 0:
            assert len(X_layers) == len(X_tables)
            query_result["X_layers"] = X_tables

        return query_result

    def read_as_anndata(
        self,
        X_name: str,
        *,
        column_names: Optional[AxisColumnNames] = None,
        X_layers: Optional[List[str]] = None,
    ) -> anndata.AnnData:
        """
        Execute the query and return result as an AnnData in-memory object.
        """
        query_result = self.read(
            X_name,
            column_names=column_names,
            X_layers=X_layers,
            use_position_indexing=True,
        )
        return make_anndata(query_result)

    def _rewrite_X_for_positional_indexing(
        self, X_tables: Dict[str, pa.Table]
    ) -> Dict[str, pa.Table]:
        """
        This is a private convenience function to convert axis dataframe to X matrix joins
        from `soma_joinid`-based joins to positionally indexed joins (like AnnData uses).

        Input is organized as:
            obs[i] annotates X[ obs[i].soma_joinid, : ]
        and
            var[j] annotates X[ :, var[j].soma_joinid ]

        Output is organized as:
            obs[i] annotates X[i, :]
        and
            var[j] annotates X[:, j]

        In addition, the `soma_joinid` column is dropped from the axis dataframes.
        """
        new_X_tables = {}
        indexer = self.get_indexer()
        for X_name, X_table in X_tables.items():
            new_X_tables[X_name] = pa.Table.from_arrays(
                (
                    indexer.obs_index(X_table["soma_dim_0"]),
                    indexer.var_index(X_table["soma_dim_1"]),
                    X_table[
                        "soma_data"
                    ].to_numpy(),  # as a side effect, consolidates chunks
                ),
                names=("_dim_0", "_dim_1", "soma_data"),
            )
        return new_X_tables

    def get_async(self) -> "AsyncExperimentQuery":
        return AsyncExperimentQuery(self)

    def get_indexer(self) -> "AxisIndexer":
        return self._indexer


class AsyncExperimentQuery:
    """
    An async proxy for ExperimentQuery, allowing use with coroutines
    [lifecycle: experimental].
    """

    query: ExperimentQuery

    def __init__(self, query: ExperimentQuery):
        self.query = query

    async def __aenter__(self) -> "AsyncExperimentQuery":
        return self

    async def __aexit__(self, *excinfo: Any) -> None:
        self.close()

    def close(self) -> None:
        self.query.close()

    @property
    def n_obs(self) -> int:
        return self.query.n_obs

    @property
    def n_vars(self) -> int:
        return self.query.n_vars

    async def obs(
        self, *, column_names: Optional[Sequence[str]] = None
    ) -> AsyncIterator[pa.Table]:
        return await to_thread(self.query.obs, column_names=column_names)

    async def var(
        self, *, column_names: Optional[Sequence[str]] = None
    ) -> AsyncIterator[pa.Table]:
        return await to_thread(self.query.var, column_names=column_names)

    async def obs_joinids(self) -> pa.Array:
        if self.query._joinids["obs"] is not None:
            return self.query._joinids["obs"]
        return await to_thread(self.query.obs_joinids)

    async def var_joinids(self) -> pa.Array:
        if self.query._joinids["var"] is not None:
            return self.query._joinids["var"]
        return await to_thread(self.query.var_joinids)

    async def X(self, layer: str, prefetch: bool = False) -> AsyncIterator[pa.Table]:
        chunk: pa.Table
        async for chunk in async_iter((i for i in self.query.X(layer, prefetch))):
            yield chunk

    async def read_as_anndata(
        self,
        X_name: str,
        *,
        column_names: Optional[AxisColumnNames] = None,
        X_layers: Optional[List[str]] = None,
    ) -> anndata.AnnData:
        return await to_thread(
            self.query.read_as_anndata,
            X_name,
            column_names=column_names,
            X_layers=X_layers,
        )


T = TypeVar("T")


async def async_iter(gen: Iterator[T]) -> AsyncIterator[T]:
    """
    Convert a generator into an async coroutine
    """
    fn = wrap_generator(gen)
    while True:
        value, done = await to_thread(fn)
        if done:
            return
        assert value is not None
        yield value


def wrap_generator(it: Iterator[T]) -> Callable[[], Tuple[Optional[T], bool]]:
    """
    Wrap a generator, making it a "normal" function that is amenable
    to running in a thread. Each time it is called, it returns a
    tuple:
        If there is another value: (next_value, False)
        If end of iteration:       (None, True)
    """

    def _next() -> Tuple[Optional[T], bool]:
        try:
            value = next(it)
            return value, False
        except StopIteration:
            return None, True

    return _next


_P = ParamSpec("_P")
_R = TypeVar("_R")


async def to_thread(
    __func: Callable[_P, _R], *args: _P.args, **kwargs: _P.kwargs
) -> _R:
    """
    Reimplementation of asyncio.to_thread, which was introduced in Py 3.9. Added
    here for support on earlier versions of Python.

    See https://docs.python.org/3/library/asyncio-task.html#asyncio.to_thread
    """
    loop = asyncio.events.get_running_loop()
    ctx = contextvars.copy_context()
    func_call = cast(
        Callable[..., _R], functools.partial(ctx.run, __func, *args, **kwargs)
    )
    return await loop.run_in_executor(None, func_call)


class AxisIndexer:
    """
    Given a query, providing index-bulding services for obs/var axis.
    """

    query: ExperimentQuery
    _obs_index: pd.Index
    _var_index: pd.Index

    def __init__(self, query: Union[ExperimentQuery, AsyncExperimentQuery]):
        if isinstance(query, AsyncExperimentQuery):
            query = query.query

        self.query = query
        self._obs_index = None
        self._var_index = None

    def obs_index(
        self, coords: Union[pa.Array, pa.ChunkedArray, npt.NDArray[np.int64]]
    ) -> npt.NDArray[np.intp]:
        if not isinstance(coords, np.ndarray):
            coords = coords.to_numpy()
        if self._obs_index is None:
            self._obs_index = pd.Index(data=self.query.obs_joinids().to_numpy())
        return cast(npt.NDArray[np.intp], self._obs_index.get_indexer(coords))

    def var_index(
        self, coords: Union[pa.Array, pa.ChunkedArray, npt.NDArray[np.int64]]
    ) -> npt.NDArray[np.intp]:
        if not isinstance(coords, np.ndarray):
            coords = coords.to_numpy()
        if self._var_index is None:
            self._var_index = pd.Index(data=self.query.var_joinids().to_numpy())
        return cast(npt.NDArray[np.intp], self._var_index.get_indexer(coords))
