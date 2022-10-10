import math
import time
from typing import Any, List, Literal, Optional, Tuple, Union, cast

import numpy as np
import pyarrow as pa
import scipy.sparse as sp
import tiledb

import tiledbsoma.eta as eta
import tiledbsoma.logging as logging
import tiledbsoma.util as util
import tiledbsoma.util_arrow as util_arrow
from tiledbsoma.util_tiledb import tiledb_result_order_from_soma_result_order

from .soma_collection import SOMACollectionBase
from .tiledb_array import TileDBArray
from .types import Matrix, NTuple, SOMADenseNdCoordinates, SOMAResultOrder


class SOMADenseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    _shape: Tuple[int]

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[SOMACollectionBase[Any]] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(uri=uri, parent=parent, ctx=ctx)

    @property
    def type(self) -> Literal["SOMADenseNdArray"]:
        return "SOMADenseNdArray"

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
    ) -> "SOMADenseNdArray":
        """
        Create a ``SOMADenseNdArray`` named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the uint64 range.
        """

        # check on shape
        if len(shape) == 0 or any(e <= 0 for e in shape):
            raise ValueError(
                "DenseNdArray shape must be non-zero length tuple of ints > 0"
            )

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

        return self

    @property
    def shape(self) -> NTuple:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        with self._tiledb_open() as A:
            return cast(NTuple, A.schema.domain.shape)

    def reshape(self, shape: NTuple) -> None:
        raise NotImplementedError("reshape operation not implemented.")

    @property
    def ndims(self) -> int:
        """
        Return number of index columns
        """
        with self._tiledb_open() as A:
            return cast(int, A.schema.domain.ndim)

    @property
    def is_sparse(self) -> Literal[False]:
        """
        Returns ``False``.
        """
        return False

    def read_tensor(
        self,
        coords: SOMADenseNdCoordinates,
        *,
        result_order: Optional[SOMAResultOrder] = "row-major",
    ) -> pa.Tensor:
        """
        Read a user-defined dense slice of the array and return as an Arrow ``Tensor``.
        """
        tiledb_result_order = tiledb_result_order_from_soma_result_order(
            result_order, accept=["column-major", "row-major"]
        )
        with self._tiledb_open("r") as A:
            target_shape = _dense_index_to_shape(coords, A.shape, result_order)
            query = A.query(return_arrow=True, order=tiledb_result_order)
            arrow_tbl = query.df[coords]
            return pa.Tensor.from_numpy(
                arrow_tbl.column("data").to_numpy().reshape(target_shape)
            )

    def read_numpy(
        self,
        coords: SOMADenseNdCoordinates,
        *,
        result_order: Optional[SOMAResultOrder] = None,
    ) -> np.ndarray:
        """
        Read a user-specified dense slice of the array and return as an Numpy ``ndarray``.
        """
        return cast(
            np.ndarray, self.read_tensor(coords, result_order=result_order).to_numpy()
        )

    def write_tensor(
        self,
        coords: SOMADenseNdCoordinates,
        values: pa.Tensor,
    ) -> None:
        """
        Write subarray, defined by ``coords`` and ``values``. Will overwrite existing
        values in the array.

        Parameters
        ----------
        coords - per-dimension tuple of scalar or slice
            Define the bounds of the subarray to be written.

        values - pyarrow.Tensor
            Define the values to be written to the subarray.  Must have same shape
            as defind by ``coords``, and the type must match the SOMADenseNdArray.
        """
        with self._tiledb_open("w") as A:
            A[coords] = values.to_numpy()

    def write_numpy(self, coords: SOMADenseNdCoordinates, values: np.ndarray) -> None:
        """ "
        Write a numpy ``ndarray`` to the user specified coordinates
        """
        self.write_tensor(coords, pa.Tensor.from_numpy(values))

    # ----------------------------------------------------------------
    #
    # TODO: this code seems obsolete given the current API. Can we port the `io` package
    # to use the above API, or move this into tiledbsoma.io as helper code?
    #
    #
    def from_matrix(self, matrix: Matrix) -> None:
        """
        Imports a matrix -- nominally ``numpy.ndarray`` -- into a TileDB array which is used for ``obsp`` and ``varp`` matrices
        """

        s = util.get_start_stamp()
        logging.log_io(None, f"{self._indent}START  WRITING")

        if self.exists():
            logging.log_io(None, f"{self._indent}Re-using existing array")
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
            f"Wrote {self.uri}",
            util.format_elapsed(s, f"{self._indent}FINISH WRITING"),
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
                        "... %7.3f%% done, ETA %s" % (chunk_percent, eta_seconds),
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


# module-private utility
def _dense_index_to_shape(
    coords: Tuple[Union[int, slice], ...],
    array_shape: Tuple[int, ...],
    result_order: Optional[SOMAResultOrder] = "row-major",
) -> Tuple[int, ...]:
    """
    Given a subarray index specified as a tuple of per-dimension slices or scalars
    (eg, ``([:], 1, [1:2])``), and the shape of the array, return the shape of
    the subarray.

    See read_tensor for usage.
    """
    shape: List[int] = []
    for n, idx in enumerate(coords):
        if type(idx) is int:
            shape.append(1)
        elif type(idx) is slice:
            start, stop, step = idx.indices(array_shape[n])
            if step != 1:
                raise ValueError("stepped slice ranges are not supported")
            shape.append(stop - start)
        else:
            raise ValueError("coordinates must be tuple of int or slice")

    if result_order == "row-major":
        return tuple(shape)

    return tuple(reversed(shape))
