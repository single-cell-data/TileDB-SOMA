import abc
import concurrent.futures
from typing import (
    TYPE_CHECKING,
    Any,
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
from typing_extensions import Literal, TypedDict

from ..dataframe import DataFrame as SOMADataFrame

if TYPE_CHECKING:
    from ..experiment import Experiment
from ..sparse_nd_array import SparseNDArray as SOMASparseNDArray
from .anndata import _make_anndata
from .axis import AxisQuery
from .eq_types import AxisColumnNames, ExperimentAxisQueryReadArrowResult

AxisJoinIds = TypedDict(
    "AxisJoinIds",
    {
        "obs": pa.Array,
        "var": pa.Array,
    },
)

MatrixAxisQuery = TypedDict(
    "MatrixAxisQuery",
    {
        "obs": AxisQuery,
        "var": AxisQuery,
    },
)


"""TEMPORARY

    with exp.query(...) as query:

        query.obs().tables() -> Iterator[pa.Table]
        query.obs().tables().concat() -> pa.Table
        query.

"""


class ExperimentAxisQuery(ContextManager["ExperimentAxisQuery"]):
    """
    ExperimentAxisQuery allows easy selection and extraction of data from a single soma.Measurement
    in a soma.Experiment, by obs/var (axis) coordinates and/or value filter [lifecycle: experimental].

    The primary use for this class is slicing Experiment ``X`` layers by obs or var value and/or coordinates.
    Slicing on SparseNDArray ``X`` matrices is support; DenseNDArray is not supported at this time.

    IMPORTANT: this class is not thread safe.

    IMPORTANT: this query class assumes it can store the full result of both axis dataframe
    queries in memory, and only provides incremental access to the underlying X NDArray. API
    features such as `n_obs` and `n_vars` codify this in the API.

    IMPORTANT: you must call `close()` on any instance of this class in order to release
    underlying resources. The ExperimentAxisQuery is a context manager, and it is recommended
    that you use the following pattern to make this easy and safe:
    ```
        with ExperimentAxisQuery(...) as query:
            ...
    ```
    """

    experiment: "Experiment"
    ms: str

    _query: MatrixAxisQuery
    _joinids: AxisJoinIds
    _indexer: "AxisIndexer"
    __threadpool: Optional[
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

        self.__threadpool = None

    @property
    def _threadpool(self) -> concurrent.futures.ThreadPoolExecutor:
        """Private threadpool cache"""
        if self.__threadpool is None:
            # TODO: the user should be able to set their own threadpool, a la asyncio's
            # loop.set_default_executor().  This is important for managing the level of
            # concurrency, etc.
            self.__threadpool = concurrent.futures.ThreadPoolExecutor()
        return self.__threadpool

    def close(self) -> None:
        """
        Cleanup and close all resources. This must be called or the thread pool
        will not be released.
        """
        if self.__threadpool is not None:
            self.__threadpool.shutdown()
            self.__threadpool = None

    def __enter__(self) -> "ExperimentAxisQuery":
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
        Private. Read the specified axis. Will load and save the resulting soma_joinids for the
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

    def obs_new(
        self, *, column_names: Optional[Sequence[str]] = None
    ) -> "DataFrameTableRead":
        query = self._query["obs"]
        return DataFrameTableRead(
            self,
            self.experiment.obs.read(
                ids=query.coords,
                value_filter=query.value_filter,
                column_names=column_names,
            ),
        )

    def var(self, *, column_names: Optional[Sequence[str]] = None) -> pa.Table:
        """Return var as an Arrow table."""
        return self._read_axis_dataframe(
            "var", self.experiment.ms[self.ms].var, column_names=column_names
        )

    def obs_joinids(self) -> pa.Array:
        """Return obs soma_joinids as an Arrow array"""
        return self._read_axis_joinids("obs", self.experiment.obs)

    def var_joinids(self) -> pa.Array:
        """Return var soma_joinids as an Arrow array"""
        return self._read_axis_joinids("var", self.experiment.ms[self.ms].var)

    @property
    def n_obs(self) -> int:
        """Return number of obs axis query results"""
        return len(self.obs_joinids())

    @property
    def n_vars(self) -> int:
        """Return number of var axis query results"""
        return len(self.var_joinids())

    def _load_joinids(self) -> None:
        """Private. Ensure joinids for both axis are in-memory."""
        futures = []
        if not self._joinids["obs"]:
            futures.append(self._threadpool.submit(self.obs_joinids))
        if not self._joinids["var"]:
            futures.append(self._threadpool.submit(self.var_joinids))
        if futures:
            concurrent.futures.wait(futures)
        assert self._joinids["obs"] is not None, "Internal error"
        assert self._joinids["var"] is not None, "Internal error"

    def _fetchX(self, X: SOMASparseNDArray) -> Iterator[pa.Table]:
        """Private helper for ``X``"""
        assert self._joinids["obs"] is not None
        assert self._joinids["var"] is not None

        obs_joinids = self._joinids["obs"]
        var_joinids = self._joinids["var"]

        if len(obs_joinids) == 0 or len(var_joinids) == 0:
            yield pa.Table.from_pylist([], schema=X.schema)
            return

        # yield for clarity
        yield from X.read_table((obs_joinids, var_joinids))

    def X(self, layer: str) -> Iterator[pa.Table]:
        """
        Return an X layer as an iterator of Arrow Tables.

        Parameters
        ----------
        layer : str
            The X layer name to return.

        Examples
        --------
        >>> with experiment.query(
        ...     "RNA",
        ...     obs_query=AxisQuery(value_filter='tissue == "lung"')
        ... ) as query:
        ...     X_result = pa.concat_tables(query.X("raw"))
        >>> X_result
        pyarrow.Table
        soma_dim_0: int64
        soma_dim_1: int64
        soma_data: float
        ----
        soma_dim_0: [[1538112,1538112,1538112,1538112,1538112,...,1538308,1538308,1538308,1538308,1538308], ...]
        soma_dim_1: [[5,19,26,34,37,...,10577,10603,10616,10617,10655], ...]
        soma_data: [[2,2,1,1,1,...,1,1,1,2,2], ...]
        """
        if not (layer and layer in self.experiment.ms[self.ms].X):
            raise ValueError("Must specify X layer")

        X = self.experiment.ms[self.ms].X[layer]
        if X.soma_type != "SOMASparseNDArray":
            raise NotImplementedError("Dense array unsupported")
        assert isinstance(X, SOMASparseNDArray)

        self._load_joinids()

        yield from self._fetchX(X)

    def X_new(self, layer: str) -> "SparseNDArrayTableRead":
        if not (layer and layer in self.experiment.ms[self.ms].X):
            raise ValueError("Must specify X layer")

        X = self.experiment.ms[self.ms].X[layer]
        if X.soma_type != "SOMASparseNDArray":
            raise NotImplementedError("Dense array unsupported")
        assert isinstance(X, SOMASparseNDArray)

        self._load_joinids()

        return SparseNDArrayTableRead(self, X.shape, self._fetchX(X))

    def _axisp_inner(
        self, axis: Literal["obs", "var"], layer: str
    ) -> Iterator[pa.Table]:
        assert axis in ["obs", "var"]
        key = f"{axis}p"

        ms = self.experiment.ms[self.ms]
        if key not in self.experiment.ms[self.ms]:
            raise ValueError(f"Measurement does not contain {key} data")

        axisp = ms.obsp if axis == "obs" else ms.varp
        if not (layer and layer in axisp):
            raise ValueError(f"Must specify '{key}' layer")

        if axisp[layer].soma_type != "SOMASparseNDArray":
            raise TypeError(f"Unexpected SOMA type stored in '{key}' layer")
        assert isinstance(axisp[layer], SOMASparseNDArray)

        self._load_joinids()
        assert self._joinids[axis] is not None

        joinids = self._joinids[axis]
        if len(joinids) == 0:
            yield pa.Table.from_pylist([], schema=axisp[layer].schema)
            return

        yield from axisp[layer].read_table((joinids, joinids))

    def obsp(self, layer: str) -> Iterator[pa.Table]:
        yield from self._axisp_inner("obs", layer)

    def varp(self, layer: str) -> Iterator[pa.Table]:
        yield from self._axisp_inner("var", layer)

    def read(
        self,
        X_name: str,
        *,
        use_position_indexing: bool = False,
        column_names: Optional[AxisColumnNames] = None,
        X_layers: Optional[List[str]] = None,
    ) -> ExperimentAxisQueryReadArrowResult:
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
            self._threadpool.submit(
                self._read_axis_dataframe,
                "obs",
                self.experiment.obs,
                column_names=column_names["obs"],
            ),
            self._threadpool.submit(
                self._read_axis_dataframe,
                "var",
                self.experiment.ms[self.ms].var,
                column_names=column_names["var"],
            ),
        )
        concurrent.futures.wait(futures)
        obs_table, var_table = (f.result() for f in futures)

        X_tables = {
            _xname: pa.concat_tables(self._fetchX(all_X_arrays[_xname]))
            for _xname in all_X_arrays
        }
        if use_position_indexing:
            X_tables = self._rewrite_X_for_positional_indexing(X_tables)

        X = X_tables.pop(X_name)
        query_result: ExperimentAxisQueryReadArrowResult = {
            "obs": obs_table,
            "var": var_table,
            "X": X,
        }
        if len(X_layers) > 0:
            assert len(X_layers) == len(X_tables)
            query_result["X_layers"] = X_tables

        return query_result

    def to_anndata(
        self,
        X_name: str,
        *,
        column_names: Optional[AxisColumnNames] = None,
        X_layers: Optional[List[str]] = None,
    ) -> anndata.AnnData:
        """
        Execute the query and return result as an AnnData in-memory object.

        Parameters
        ----------
        X_name : str
            The name of the X layer to read and return in the AnnData.X slot
        column_names : Optional[dict]
            Specify which column names in ``var`` and ``obs`` dataframes to read and return.
        X_layers : Optional[List[str]]
            Addtional X layers read read and return in the AnnData.layers slot

        Examples
        --------
        >>> with exp.query_by_axis("RNA", obs_query=AxisQuery(value_filter='tissue == "lung"')) as query:
        ...     ad = query.to_anndata("raw", column_names={"obs": ["cell_type", "tissue"]})
        >>> ad
        AnnData object with n_obs × n_vars = 127310 × 52373
        obs: 'cell_type', 'tissue'
        var: 'soma_joinid', 'feature_id', 'feature_name', 'feature_length'
        """
        query_result = self.read(
            X_name,
            column_names=column_names,
            X_layers=X_layers,
            use_position_indexing=True,
        )
        return _make_anndata(query_result)

    def _rewrite_X_for_positional_indexing(
        self, X_tables: Dict[str, pa.Table]
    ) -> Dict[str, pa.Table]:
        """
        Private convenience function to convert axis dataframe to X matrix joins
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

    def get_indexer(self) -> "AxisIndexer":
        return self._indexer


class AxisIndexer:
    """
    Given a query, providing index-bulding services for obs/var axis.
    """

    query: ExperimentAxisQuery
    _obs_index: pd.Index
    _var_index: pd.Index

    def __init__(self, query: ExperimentAxisQuery):
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


# ********* XXX: Begin temporary code **********
#

"""
TODO: This code is temporary mock until such time as the core classes implement
these interfaces (at which point we can remove all of this).
"""
_T = TypeVar("_T")


class ReadIter(Iterator[_T], metaclass=abc.ABCMeta):
    """SparseRead result iterator allowing users to flatten the iteration."""

    # __iter__ is already implemented as `return self` in Iterator
    # SOMA implementations must implement __next__.

    # XXX: Considering the name "flat" here too.
    @abc.abstractmethod
    def flat(self) -> _T:
        """Returns all the requested data in a single operation."""
        raise NotImplementedError()


class SparseRead(metaclass=abc.ABCMeta):
    """Intermediate type to allow users to format when reading a sparse array."""

    def coos(self) -> ReadIter[pa.SparseCOOTensor]:
        raise NotImplementedError()

    def cscs(self) -> ReadIter[pa.SparseCSCMatrix]:
        raise NotImplementedError()

    def csrs(self) -> ReadIter[pa.SparseCSRMatrix]:
        raise NotImplementedError()

    def record_batches(self) -> ReadIter[pa.RecordBatch]:
        raise NotImplementedError()

    def tables(self) -> ReadIter[pa.Table]:
        raise NotImplementedError()


class ReadTableIter(ReadIter[pa.Table]):
    def __init__(self, tables: Iterator[pa.Table]):
        self._tables = tables

    def __next__(self) -> pa.Table:
        return next(self._tables)

    def flat(self) -> pa.Table:
        return pa.concat_tables(self._tables)


class ReadCOOTensorIter(ReadIter[pa.SparseCOOTensor]):
    def __init__(self, shape: Tuple[int, ...], tables: Iterator[pa.Table]):
        self._tables = tables
        self._shape = shape

    def _make_coo_tensor(self, arrow_tbl: pa.Table) -> pa.SparseCOOTensor:
        coo_data = arrow_tbl.column("soma_data").to_numpy()
        coo_coords = np.array(
            [
                arrow_tbl.column(f"soma_dim_{n}").to_numpy()
                for n in range(len(self._shape))
            ]
        ).T
        return pa.SparseCOOTensor.from_numpy(coo_data, coo_coords, shape=self._shape)

    def __next__(self) -> pa.SparseCOOTensor:
        return self._make_coo_tensor(next(self._tables))

    def flat(self) -> pa.SparseCOOTensor:
        return self._make_coo_tensor(pa.concat_tables(self._tables))


class DataFrameTableRead(SparseRead):
    def __init__(self, query: ExperimentAxisQuery, tables: Iterator[pa.Table]):
        self._query = query
        self._tables = tables

    def tables(self) -> ReadIter[pa.Table]:
        return ReadTableIter(self._tables)


class SparseNDArrayTableRead(SparseRead):
    def __init__(
        self,
        query: ExperimentAxisQuery,
        shape: Tuple[int, ...],
        tables: Iterator[pa.Table],
    ):
        self._query = query
        self._tables = tables
        self._shape = shape

    def tables(self) -> ReadIter[pa.Table]:
        return ReadTableIter(self._tables)

    def coo(self) -> ReadIter[pa.SparseCOOTensor]:
        return ReadCOOTensorIter(self._shape, self._tables)


#
# ********* XXX End temporary copy **********
