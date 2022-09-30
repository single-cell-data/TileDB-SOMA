import math
import time
from typing import Any, List, Literal, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb

import tiledbsoma.eta as eta
import tiledbsoma.logging as logging
import tiledbsoma.util as util
import tiledbsoma.util_arrow as util_arrow
import tiledbsoma.util_tiledb as util_tiledb

from .soma_collection import SOMACollectionBase
from .tiledb_array import TileDBArray
from .types import Matrix, NTuple


class SOMADenseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    _shape: Tuple[int]

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[SOMACollectionBase[Any]] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(uri=uri, name=name, parent=parent, ctx=ctx)

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
    ) -> None:
        """
        Create a ``SOMADenseNdArray`` named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the uint64 range.
        """

        # checks on shape
        assert len(shape) > 0
        for e in shape:
            assert e > 0

        if not pa.types.is_primitive(type):
            raise TypeError(
                "Unsupported type - SOMADenseNdArray only supports primtive Arrow types"
            )

        level = self._tiledb_platform_config.string_dim_zstd_level

        dims = []
        for e in shape:
            dim = tiledb.Dim(
                # Use tiledb default names like ``__dim_0``
                domain=(0, e - 1),
                tile=min(e, 2048),  # TODO: PARAMETERIZE
                dtype=np.uint64,
                filters=[tiledb.ZstdFilter(level=level)],
            )
            dims.append(dim)
        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = [
            tiledb.Attr(
                name="data",
                dtype=util_arrow.tiledb_type_from_arrow_type(type),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
        ]

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=False,
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=100000,
            cell_order="row-major",
            tile_order="row-major",
            ctx=self._ctx,
        )

        tiledb.Array.create(self._uri, sch, ctx=self._ctx)

        self._common_create()  # object-type metadata etc

    def __repr__(self) -> str:
        """
        Default display of ``SOMADenseNdArray``.
        """
        return "\n".join(self._repr_aux())

    def _repr_aux(self) -> Sequence[str]:
        if not self.exists():
            return ["Unpopulated"]
        lines = [
            self.name
            + " "
            + self.__class__.__name__
            # Pending https://github.com/single-cell-data/TileDB-SOMA/issues/302
            # + " "
            # + str(self.shape)
        ]
        return lines

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        # TODO: cache read
        # return self._shape
        with self._tiledb_open() as A:
            # mypy says:
            # error: Returning Any from function declared to return "Tuple[int]"  [no-any-return]
            return A.schema.domain.shape  # type: ignore

    @property
    def ndims(self) -> int:
        """
        Return number of index columns
        """
        with self._tiledb_open() as A:
            # mypy says:
            # Returning Any from function declared to return "int"  [no-any-return]
            return A.schema.domain.ndim  # type: ignore

    @property
    def is_sparse(self) -> Literal[False]:
        """
        Returns ``False``.
        """
        return False

    def read(
        self,
        # TODO: partitions: Optional[SOMAReadPartitions] = None,,
        row_ids: Optional[Sequence[int]] = None,
        col_ids: Optional[Sequence[int]] = None,
        result_order: Optional[str] = None,
    ) -> Any:  # TODO: Iterator[DenseReadResult]
        """
        Read a user-specified subset of the object, and return as one or more Arrow.Tensor.

        :param ids: per-dimension slice, expressed as a scalar, a range, or a list of both.

        :param partitions: an optional [``SOMAReadPartitions``](#SOMAReadPartitions) hint to indicate how results should be organized.

        :param result_order: order of read results. Can be one of row-major or column-major.

        The ``read`` operation will return a language-specific iterator over one or more Arrow Tensor objects and information describing them, allowing the incremental processing of results larger than available memory. The actual iterator used is delegated to language-specific SOMA specs. The ``DenseReadResult`` should include:

        * The coordinates of the slice (e.g., origin, shape)
        * an Arrow.Tensor with the slice values
        """
        tiledb_result_order = (
            util_tiledb.tiledb_result_order_from_soma_result_order_indexed(result_order)
        )

        with self._tiledb_open("r") as A:
            query = A.query(
                return_arrow=True,
                return_incomplete=True,
                order=tiledb_result_order,
            )

            if row_ids is None:
                if col_ids is None:
                    iterator = query.df[:, :]
                else:
                    iterator = query.df[:, col_ids]
            else:
                if col_ids is None:
                    iterator = query.df[row_ids, :]
                else:
                    iterator = query.df[row_ids, col_ids]

            for df in iterator:
                batches = df.to_batches()
                for batch in batches:
                    yield batch

    def read_as_pandas(
        self,
        *,
        row_ids: Optional[Sequence[int]] = None,
        col_ids: Optional[Sequence[int]] = None,
        set_index: Optional[bool] = False,
    ) -> pd.DataFrame:
        """
        TODO: comment
        """
        with self._tiledb_open() as A:
            query = A.query(return_incomplete=True)

            if row_ids is None:
                if col_ids is None:
                    iterator = query.df[:, :]
                else:
                    iterator = query.df[:, col_ids]
            else:
                if col_ids is None:
                    iterator = query.df[row_ids, :]
                else:
                    iterator = query.df[row_ids, col_ids]

            for df in iterator:
                # Make this opt-in only.  For large arrays, this df.set_index is time-consuming
                # so we should not do it without direction.
                if set_index:
                    df.set_index(self._tiledb_dim_names(), inplace=True)
                yield df

    def read_all(
        self,
        *,
        # TODO: find the right syntax to get the typechecker to accept args like ``ids=slice(0,10)``
        # ids: Optional[Union[Sequence[int], Slice]] = None,
        row_ids: Optional[Sequence[int]] = None,
        col_ids: Optional[Sequence[int]] = None,
        result_order: Optional[str] = None,
        # TODO: batch_size
        # TODO: partition,
        # TODO: batch_format,
        # TODO: platform_config,
    ) -> pa.RecordBatch:
        """
        This is a convenience method around ``read``. It iterates the return value from ``read`` and returns a concatenation of all the record batches found. Its nominal use is to simply unit-test cases.
        """
        return util_arrow.concat_batches(
            self.read(
                row_ids=row_ids,
                col_ids=col_ids,
                result_order=result_order,
            )
        )

    def read_as_pandas_all(
        self,
        *,
        row_ids: Optional[Sequence[int]] = None,
        col_ids: Optional[Sequence[int]] = None,
        set_index: Optional[bool] = False,
    ) -> pa.RecordBatch:
        """
        This is a convenience method around ``read_as_pandas``. It iterates the return value from ``read_as_pandas`` and returns a concatenation of all the record batches found. Its nominal use is to simply unit-test cases.
        """
        dataframes = []
        generator = self.read_as_pandas(
            row_ids=row_ids,
            col_ids=col_ids,
            set_index=set_index,
        )
        for dataframe in generator:
            dataframes.append(dataframe)
        return pd.concat(dataframes)

    def write(
        self,
        # TODO: rework callsites with regard to the very latest spec rev
        # coords: Union[tuple, tuple[slice], NTuple, List[int]],
        coords: Any,
        values: pa.Tensor,
    ) -> None:
        """
        Write an Arrow.Tensor to the persistent object. As duplicate index values are not allowed, index values already present in the object are overwritten and new index values are added.

        :param coords: location at which to write the tensor

        :param values: an Arrow.Tensor containing values to be written. The type of elements in ``values`` must match the type of the SOMADenseNdArray.
        """

        with self._tiledb_open("w") as A:
            A[coords] = values.to_numpy()

    # ----------------------------------------------------------------
    def from_matrix(self, matrix: Matrix) -> None:
        """
        Imports a matrix -- nominally ``numpy.ndarray`` -- into a TileDB array which is used for ``obsp`` and ``varp`` matrices
        """

        s = util.get_start_stamp()
        logging.log_io(None, f"{self._indent}START  WRITING {self._nested_name}")

        if self.exists():
            logging.log_io(
                None, f"{self._indent}Re-using existing array {self._nested_name}"
            )
        else:
            self._create_empty_array(
                matrix_dtype=matrix.dtype,
                num_rows=matrix.shape[0],
                num_cols=matrix.shape[1],
            )

        if not self._tiledb_platform_config.write_X_chunked:
            self._ingest_data_whole(matrix)
        else:
            self._ingest_data_dense_rows_chunked(matrix)

        self._common_create()  # object-type metadata etc

        logging.log_io(
            f"Wrote {self._nested_name}",
            util.format_elapsed(s, f"{self._indent}FINISH WRITING {self._nested_name}"),
        )

    # ----------------------------------------------------------------
    def _create_empty_array(
        self, *, matrix_dtype: np.dtype, num_rows: int, num_cols: int
    ) -> None:
        """
        Create a TileDB 2D dense array with int dimensions and a single attribute.
        """

        dom = tiledb.Domain(
            tiledb.Dim(
                domain=(0, num_rows - 1),
                dtype=np.uint64,
                # TODO: filters=[tiledb.RleFilter()],
            ),
            tiledb.Dim(
                domain=(0, num_cols - 1),
                dtype=np.uint64,
                # TODO: filters=[tiledb.ZstdFilter(level=level)],
            ),
            ctx=self._ctx,
        )

        attrs = tiledb.Attr(
            dtype=matrix_dtype,
            filters=[tiledb.ZstdFilter()],
            ctx=self._ctx,
        )

        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=(attrs,),
            sparse=True,
            allows_duplicates=self._tiledb_platform_config.allows_duplicates,
            offsets_filters=[
                tiledb.DoubleDeltaFilter(),
                tiledb.BitWidthReductionFilter(),
                tiledb.ZstdFilter(),
            ],
            capacity=self._tiledb_platform_config.X_capacity,
            cell_order=self._tiledb_platform_config.X_cell_order,
            tile_order=self._tiledb_platform_config.X_tile_order,
            ctx=self._ctx,
        )

        tiledb.Array.create(self.uri, sch, ctx=self._ctx)

    def _ingest_data_whole(self, matrix: np.ndarray) -> None:
        raise NotImplementedError()

    def _ingest_data_dense_rows_chunked(
        self,
        matrix: np.ndarray,
    ) -> None:
        """
        Convert dense matrix to coo_matrix chunkwise and ingest into TileDB.

        :param uri: TileDB URI of the array to be written.
        :param matrix: dense matrix.
        """

        nrow, ncol = matrix.shape

        s = util.get_start_stamp()
        logging.log_io(
            None,
            f"{self._indent}START  ingest",
        )

        eta_tracker = eta.Tracker()
        with tiledb.open(self.uri, mode="w", ctx=self._ctx) as A:

            i = 0
            while i < nrow:
                t1 = time.time()
                # Find a number of dense rows which will result in a desired nnz for the chunk,
                # rounding up to the nearest integer. Example: goal_chunk_nnz is 120. ncol is 50;
                # 120/50 rounds up to 3; take 3 rows per chunk.
                chunk_size = int(
                    math.ceil(self._tiledb_platform_config.goal_chunk_nnz / ncol)
                )
                i2 = i + chunk_size

                # Convert the chunk to a COO matrix.
                chunk = matrix[i:i2]
                chunk_coo = sp.csr_matrix(chunk).tocoo()

                # Python ranges are (lo, hi) with lo inclusive and hi exclusive. But saying that
                # makes us look buggy if we say we're ingesting chunk 0:18 and then 18:32.
                # Instead, print doubly-inclusive lo..hi like 0..17 and 18..31.
                chunk_percent = min(100, 100 * (i2 - 1) / nrow)
                logging.log_io(
                    None,
                    "%sSTART  chunk rows %d..%d of %d (%.3f%%), nnz=%d"
                    % (
                        self._indent,
                        i,
                        i2 - 1,
                        nrow,
                        chunk_percent,
                        chunk_coo.nnz,
                    ),
                )

                # Write a TileDB fragment
                A[chunk_coo.row + i, chunk_coo.col] = chunk_coo.data

                t2 = time.time()
                chunk_seconds = t2 - t1
                eta_seconds = eta_tracker.ingest_and_predict(
                    chunk_percent, chunk_seconds
                )

                if chunk_percent < 100:
                    logging.log_io(
                        "... %s %7.3f%% done, ETA %s"
                        % (self._nested_name, chunk_percent, eta_seconds),
                        "%sFINISH chunk in %.3f seconds, %7.3f%% done, ETA %s"
                        % (self._indent, chunk_seconds, chunk_percent, eta_seconds),
                    )

                i = i2

        logging.log_io(
            None,
            util.format_elapsed(
                s,
                f"{self._indent}FINISH ingest",
            ),
        )
