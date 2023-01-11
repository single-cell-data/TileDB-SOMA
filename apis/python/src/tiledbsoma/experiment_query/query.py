import concurrent.futures
from typing import (
    TYPE_CHECKING,
    Any,
    ContextManager,
    Dict,
    Optional,
    Sequence,
    Tuple,
    Union,
    cast,
)

import anndata
import attrs
import numpy as np
import numpy.typing as npt
import pandas as pd
import pyarrow as pa
import scipy.sparse as sparse
from typing_extensions import Literal, TypedDict

from .. import somacore  # to be replaced by somacore package, when available
from ..dataframe import DataFrame as SOMADataFrame
from ..sparse_nd_array import SparseNDArrayRead

if TYPE_CHECKING:
    from ..experiment import Experiment

from ..sparse_nd_array import SparseNDArray as SOMASparseNDArray
from .axis import AxisQuery


@attrs.define
class _AxisJoinIds:
    """Private: cache per-axis join ids in the query"""

    obs: Optional[pa.Array] = attrs.field(default=None)
    var: Optional[pa.Array] = attrs.field(default=None)


@attrs.define
class _MatrixAxisQuery:
    """Private: store per-axis user query definition"""

    obs: AxisQuery
    var: AxisQuery


class AxisColumnNames(TypedDict, total=False):
    """Specify column names for the ExperimentAxisQuery API read operations"""

    obs: Optional[Sequence[str]]  # None is all columns
    var: Optional[Sequence[str]]  # None is all columns


class ExperimentAxisQuery(ContextManager["ExperimentAxisQuery"]):
    """
    ExperimentAxisQuery allows easy selection and extraction of data from a single soma.Measurement
    in a soma.Experiment, by obs/var (axis) coordinates and/or value filter [lifecycle: experimental].

    The primary use for this class is slicing Experiment ``X`` layers by obs or var value and/or coordinates.
    Slicing on SparseNDArray ``X`` matrices is supported; DenseNDArray is not supported at this time.

    IMPORTANT: this class is not thread-safe.

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
        self.measurement_name = measurement_name

        self._matrix_axis_query = _MatrixAxisQuery(
            obs=obs_query if obs_query is not None else AxisQuery(),
            var=var_query if var_query is not None else AxisQuery(),
        )
        self._joinids = _AxisJoinIds()
        self._indexer: "AxisIndexer" = AxisIndexer(self)
        self._threadpool_: Optional[concurrent.futures.ThreadPoolExecutor] = None

    def close(self) -> None:
        """
        Clean up and close all resources. This must be called or the thread pool
        will not be released.
        """
        if self._threadpool_:
            self._threadpool_.shutdown()
            self._threadpool_ = None

    def __enter__(self) -> "ExperimentAxisQuery":
        return self

    def __exit__(self, *excinfo: Any) -> None:
        self.close()

    def obs(
        self, *, column_names: Optional[Sequence[str]] = None
    ) -> somacore.ReadIter[pa.Table]:
        """Return obs as an Arrow table iterator."""
        obs_query = self._matrix_axis_query.obs
        return self.experiment.obs.read(
            ids=obs_query.coords,
            value_filter=obs_query.value_filter,
            column_names=column_names,
        )

    def var(
        self, *, column_names: Optional[Sequence[str]] = None
    ) -> somacore.ReadIter[pa.Table]:
        """Return obs as an Arrow table iterator."""
        var_query = self._matrix_axis_query.var
        return self.experiment.ms[self.measurement_name].var.read(
            ids=var_query.coords,
            value_filter=var_query.value_filter,
            column_names=column_names,
        )

    def obs_joinids(self) -> pa.Array:
        """Return obs soma_joinids as an Arrow array."""
        if self._joinids.obs is None:
            obs_query = self._matrix_axis_query.obs
            self._joinids.obs = (
                self.experiment.obs.read(
                    ids=obs_query.coords,
                    value_filter=obs_query.value_filter,
                    column_names=["soma_joinid"],
                )
                .concat()
                .column("soma_joinid")
                .combine_chunks()
            )
        return self._joinids.obs

    def var_joinids(self) -> pa.Array:
        """Return var soma_joinids as an Arrow array."""
        if self._joinids.var is None:
            var_query = self._matrix_axis_query.var
            self._joinids.var = (
                self.experiment.ms[self.measurement_name]
                .var.read(
                    ids=var_query.coords,
                    value_filter=var_query.value_filter,
                    column_names=["soma_joinid"],
                )
                .concat()
                .column("soma_joinid")
                .combine_chunks()
            )
        return self._joinids.var

    @property
    def n_obs(self) -> int:
        """Return number of obs axis query results"""
        return len(self.obs_joinids())

    @property
    def n_vars(self) -> int:
        """Return number of var axis query results"""
        return len(self.var_joinids())

    def X(self, layer_name: str) -> SparseNDArrayRead:
        """
        Return an X layer as SparseNDArrayRead.

        Parameters
        ----------
        layer_name : str
            The X layer name to return.

        Examples
        --------
        >>> with experiment.query(
        ...     "RNA",
        ...     obs_query=AxisQuery(value_filter='tissue == "lung"')
        ... ) as query:
        ...     X_result = query.X("raw").tables().concat()
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
        if not (
            layer_name and layer_name in self.experiment.ms[self.measurement_name].X
        ):
            raise ValueError("Must specify X layer name")

        x_layer = self.experiment.ms[self.measurement_name].X[layer_name]
        if not isinstance(x_layer, SOMASparseNDArray):
            raise TypeError("Dense array unsupported")

        self._load_joinids()
        return x_layer.read((self._joinids.obs, self._joinids.var))

    def obsp(self, layer: str) -> SparseNDArrayRead:
        """Return an ``obsp`` layer as a SparseNDArrayRead"""
        return self._axisp_inner("obs", layer)

    def varp(self, layer: str) -> SparseNDArrayRead:
        """Return an ``varp`` layer as a SparseNDArrayRead"""
        return self._axisp_inner("var", layer)

    def _read_axis_dataframe(
        self,
        axis: Literal["obs", "var"],
        axis_df: SOMADataFrame,
        axis_query: AxisQuery,
        *,
        column_names: Optional[Sequence[str]],
    ) -> pa.Table:
        """
        Private. Read the specified axis. Will load and save the resulting soma_joinids for the
        axis, if they are not already known.
        """
        need_joinids = getattr(self._joinids, axis, None) is None
        query_columns = column_names
        if (
            need_joinids
            and column_names is not None
            and "soma_joinid" not in column_names
        ):
            query_columns = ["soma_joinid"] + list(column_names)

        arrow_table = axis_df.read(
            ids=axis_query.coords,
            value_filter=axis_query.value_filter,
            column_names=query_columns,
        ).concat()

        if need_joinids:
            setattr(
                self._joinids, axis, arrow_table.column("soma_joinid").combine_chunks()
            )
        assert getattr(self._joinids, axis, None) is not None

        if column_names is not None:
            arrow_table = arrow_table.select(column_names)
        return arrow_table

    def _read_both_axes(
        self,
        column_names: AxisColumnNames,
    ) -> Tuple[pa.Table, pa.Table]:
        """Private. Read both axes in its entirety, ensure that soma_joinid is retained."""
        futures = (
            self._threadpool.submit(
                self._read_axis_dataframe,
                "obs",
                self.experiment.obs,
                self._matrix_axis_query.obs,
                column_names=column_names.get("obs", None),
            ),
            self._threadpool.submit(
                self._read_axis_dataframe,
                "var",
                self.experiment.ms[self.measurement_name].var,
                self._matrix_axis_query.var,
                column_names=column_names.get("var", None),
            ),
        )
        concurrent.futures.wait(futures)
        obs_table, var_table = (f.result() for f in futures)
        return obs_table, var_table

    def read(
        self,
        X_name: str,
        *,
        column_names: Optional[AxisColumnNames] = None,
        X_layers: Sequence[str] = (),
    ) -> "ExperimentAxisQueryReadArrowResult":
        """
        Read the _entire_ query result into (in-memory) Arrow Tables. This is a
        low-level routine intended to be used by loaders for other in-core
        formats, such as AnnData, which can be created from the resulting Tables.

        Parameters
        ----------
        X_name : str
            The name of the X layer to read and return in the AnnData.X slot
        column_names : Optional[dict]
            Specify which column names in ``var`` and ``obs`` dataframes to read and return.
        X_layers : Optional[Sequence[str]]
            Addtional X layers read read and return in the AnnData.layers slot

        """
        if column_names is None:
            column_names = AxisColumnNames(obs=None, var=None)
        X_collection = self.experiment.ms[self.measurement_name].X
        all_X_names = [X_name] + list(X_layers)
        all_X_arrays: Dict[str, SOMASparseNDArray] = {}
        for _xname in all_X_names:
            if not isinstance(_xname, str) or not _xname:
                raise ValueError("X layer names must be specified as a string.")
            if _xname not in X_collection:
                raise ValueError("Unknown X layer name")
            x_array = X_collection[_xname]
            if not isinstance(x_array, SOMASparseNDArray):
                raise NotImplementedError("Dense array unsupported")
            all_X_arrays[_xname] = x_array

        obs_table, var_table = self._read_both_axes(column_names)

        X_tables = {
            # TODO: could also be done concurrently
            _xname: all_X_arrays[_xname]
            .read((self.obs_joinids(), self.var_joinids()))
            .tables()
            .concat()
            for _xname in all_X_arrays
        }

        X = X_tables.pop(X_name)
        query_result = ExperimentAxisQueryReadArrowResult(
            obs=obs_table, var=var_table, X=X, X_layers=X_tables
        )
        return query_result

    def to_anndata(
        self,
        X_name: str,
        *,
        column_names: Optional[AxisColumnNames] = None,
        X_layers: Sequence[str] = (),
    ) -> anndata.AnnData:
        """
        Execute the query and return result as an AnnData in-memory object.

        Parameters
        ----------
        X_name : str
            The name of the X layer to read and return in the AnnData.X slot
        column_names : Optional[dict]
            Specify which column names in ``var`` and ``obs`` dataframes to read and return.
        X_layers : Optional[Sequence[str]]
            Addtional X layers read read and return in the AnnData.layers slot

        Examples
        --------
        >>> with exp.axis_query("RNA", obs_query=AxisQuery(value_filter='tissue == "lung"')) as query:
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
        )

        # AnnData uses positional indexing
        query_result.X = self._rewrite_matrix_for_positional_indexing(query_result.X)
        query_result.X_layers = {
            name: self._rewrite_matrix_for_positional_indexing(matrix)
            for name, matrix in query_result.X_layers.items()
        }
        return query_result.to_anndata()

    def get_indexer(self) -> "AxisIndexer":
        return self._indexer

    """
    Private helper methods
    """

    @property
    def _threadpool(self) -> concurrent.futures.ThreadPoolExecutor:
        """Private threadpool cache"""
        if self._threadpool_ is None:
            # TODO: the user should be able to set their own threadpool, a la asyncio's
            # loop.set_default_executor().  This is important for managing the level of
            # concurrency, etc.
            self._threadpool_ = concurrent.futures.ThreadPoolExecutor()
        return self._threadpool_

    def _load_joinids(self) -> None:
        """Private. Ensure joinids for both axis are in-memory."""
        futures = []
        if self._joinids.obs is None:
            futures.append(self._threadpool.submit(self.obs_joinids))
        if self._joinids.var is None:
            futures.append(self._threadpool.submit(self.var_joinids))
        if futures:
            concurrent.futures.wait(futures)
        assert self._joinids.obs is not None, "Internal error"
        assert self._joinids.var is not None, "Internal error"

    def _axisp_inner(
        self,
        axis: Literal["obs", "var"],
        layer: str,
    ) -> SparseNDArrayRead:
        assert axis in ["obs", "var"]
        key = f"{axis}p"

        ms = self.experiment.ms[self.measurement_name]
        if key not in self.experiment.ms[self.measurement_name]:
            raise ValueError(f"Measurement does not contain {key} data")

        axisp = ms.obsp if axis == "obs" else ms.varp
        if not (layer and layer in axisp):
            raise ValueError(f"Must specify '{key}' layer")
        if not isinstance(axisp[layer], SOMASparseNDArray):
            raise TypeError(f"Unexpected SOMA type stored in '{key}' layer")

        self._load_joinids()
        joinids = getattr(self._joinids, axis)
        return axisp[layer].read((joinids, joinids))

    def _rewrite_matrix_for_positional_indexing(self, x_table: pa.Table) -> pa.Table:
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
        indexer = self.get_indexer()
        return pa.Table.from_arrays(
            (
                indexer.obs_index(x_table["soma_dim_0"]),
                indexer.var_index(x_table["soma_dim_1"]),
                x_table[
                    "soma_data"
                ].to_numpy(),  # as a side effect, consolidates chunks
            ),
            names=("_dim_0", "_dim_1", "soma_data"),
        )


@attrs.define
class AxisIndexer:
    """
    Given a query, providing index-bulding services for obs/var axis.
    """

    query: ExperimentAxisQuery
    _obs_index: Optional[pd.Index] = attrs.field(init=False, default=None)
    _var_index: Optional[pd.Index] = attrs.field(init=False, default=None)

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


@attrs.define
class ExperimentAxisQueryReadArrowResult:
    """Return type for the ExperimentAxisQuery.read() method"""

    obs: pa.Table
    """Experiment.obs query slice, as an Arrow Table"""
    var: pa.Table
    """Experiment.ms[...].var query slice, as an Arrow Table"""
    X: pa.Table
    """Experiment.ms[...].X[...] query slice, as an Arrow Table"""
    X_layers: Dict[str, pa.Table] = attrs.field(factory=dict)
    """Any additional X layers requested, as Arrow Table(s)"""

    def to_anndata(self) -> anndata.AnnData:
        """Convert to AnnData"""
        obs = self.obs.to_pandas()
        obs.index = obs.index.map(str)

        var = self.var.to_pandas()
        var.index = var.index.map(str)

        shape = (len(obs), len(var))

        x = self.X
        if x is not None:
            x = _arrow_to_scipy_csr(x, shape)

        layers = {
            name: _arrow_to_scipy_csr(table, shape)
            for name, table in self.X_layers.items()
        }
        return anndata.AnnData(X=x, obs=obs, var=var, layers=(layers or None))


def _arrow_to_scipy_csr(
    arrow_table: pa.Table, shape: Tuple[int, int]
) -> sparse.csr_matrix:
    """
    Private utility which converts a table repesentation of X to a CSR matrix.

    IMPORTANT: by convention, assumes that the data is positionally indexed (hence
    the use of _dim_{n} rather than soma_dim{n}).

    See query.py::_rewrite_X_for_positional_indexing for more info.
    """
    assert "_dim_0" in arrow_table.column_names, "X must be positionally indexed"
    assert "_dim_1" in arrow_table.column_names, "X must be positionally indexed"

    return sparse.csr_matrix(
        (
            arrow_table["soma_data"].to_numpy(),
            (arrow_table["_dim_0"].to_numpy(), arrow_table["_dim_1"].to_numpy()),
        ),
        shape=shape,
    )
