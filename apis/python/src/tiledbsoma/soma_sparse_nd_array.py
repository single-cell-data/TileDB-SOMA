from typing import Any, Iterator, List, Literal, Optional, Union, cast

import numpy as np
import pandas as pd
import pyarrow as pa
import scipy.sparse as sp
import tiledb

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from . import util, util_arrow
from .soma_collection import SOMACollectionBase
from .tiledb_array import TileDBArray
from .tiledb_platform_config import TileDBPlatformConfig
from .types import NTuple, SOMASparseNdCoordinates


class SOMASparseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    def __init__(
        self,
        uri: str,
        *,
        parent: Optional[SOMACollectionBase[Any]] = None,
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """

        super().__init__(
            uri=uri,
            parent=parent,
            tiledb_platform_config=tiledb_platform_config,
            ctx=ctx,
        )

    @property
    def soma_type(self) -> Literal["SOMASparseNdArray"]:
        return "SOMASparseNdArray"

    def create(
        self,
        type: pa.DataType,
        shape: Union[NTuple, List[int]],
    ) -> "SOMASparseNdArray":
        """
        Create a ``SOMASparseNdArray`` named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in the positive int64 range.
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
        for n, e in enumerate(shape):
            dim = tiledb.Dim(
                name=f"soma_dim_{n}",
                domain=(0, e - 1),
                tile=min(e, 2048),  # TODO: PARAMETERIZE,
                dtype=np.int64,
                filters=[tiledb.ZstdFilter(level=level)],
            )
            dims.append(dim)
        dom = tiledb.Domain(dims, ctx=self._ctx)

        attrs = [
            tiledb.Attr(
                name="soma_data",
                dtype=util_arrow.tiledb_type_from_arrow_type(type),
                filters=[tiledb.ZstdFilter()],
                ctx=self._ctx,
            )
        ]

        # TODO: code-dedupe w/ regard to SOMADenseNdArray. The two creates are
        # almost identical & could share a common parent-class _create() method.
        sch = tiledb.ArraySchema(
            domain=dom,
            attrs=attrs,
            sparse=True,
            allows_duplicates=False,
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
        Return length of each dimension, always a list of length ``ndim``
        """
        with self._tiledb_open() as A:
            return cast(NTuple, A.schema.domain.shape)

    def reshape(self, shape: NTuple) -> None:
        raise NotImplementedError("reshape operation not implemented.")

    @property
    def ndim(self) -> int:
        """
        Return number of index columns
        """
        return len(self.shape)

    @property
    def is_sparse(self) -> Literal[True]:
        """
        Returns ``True``.
        """
        return True

    @property
    def nnz(self) -> int:
        """
        Return the number of stored values in the array, including explicitly stored zeros.
        """
        return cast(int, clib.SOMAReader(self.uri).nnz())

    def read_sparse_tensor(
        self,
        coords: SOMASparseNdCoordinates,
        *,
        format: Literal["coo", "csr", "csc"] = "coo",
    ) -> Iterator[Union[pa.SparseCOOTensor, pa.SparseCSCMatrix, pa.SparseCSRMatrix]]:
        """
        Read a use-defined slice of the SparseNdArray and return as an Arrow sparse tensor.

        Parameters
        ----------
        coords : Tuple[Union[int, slice, Tuple[int, ...], List[int], pa.IntegerArray], ...]
            Per-dimension tuple of scalar, slice, sequence of scalar or Arrow IntegerArray
            Arrow arrays currently uninimplemented.

        format - Literal["coo", "csr", "csc"]
            Requested return format:
            * ``coo`` - return an Arrow SparseCOOTensor (default)
            * ``csr`` - return an Arrow SparseCSRMatrix
            * ``csc`` - return an Arrow SparseCSCMatrix


        Returns
        -------
        The requested data in the sparse tensor format specified.
        """

        if format != "coo" and self.ndim != 2:
            raise ValueError(f"Format {format} only supported for 2D SparseNdArray")
        if format not in ("coo", "csr", "csc"):
            raise NotImplementedError("format not implemented")

        with self._tiledb_open("r") as A:
            shape = A.shape
        for arrow_tbl in self.read_table(coords):
            """
            In PyArrow 9.0.0, there is a bug preventing the creation of "empty"
            (zero element) SparseCOOTensor objects.

            See https://issues.apache.org/jira/browse/ARROW-17933

            Just stop the iteration when we run out of results. The caller must be
            prepared to have a StopIteration, rather than an empty tensor, as the result of
            an empty query.
            """
            if arrow_tbl.num_rows == 0:
                return

            if format == "coo":

                coo_data = arrow_tbl.column("soma_data").to_numpy()
                coo_coords = np.array(
                    [
                        arrow_tbl.column(f"soma_dim_{n}").to_numpy()
                        for n in range(self.ndim)
                    ]
                ).T
                yield pa.SparseCOOTensor.from_numpy(coo_data, coo_coords, shape=shape)

            elif format in ("csr", "csc"):
                # Temporary: as these must be 2D, convert to scipy COO and use
                # scipy to perform conversions.  C++ reader will be nicer!
                data = arrow_tbl.column("soma_data").to_numpy()
                row = arrow_tbl.column("soma_dim_0").to_numpy()
                col = arrow_tbl.column("soma_dim_1").to_numpy()
                scipy_coo = sp.coo_array((data, (row, col)), shape=shape)
                if format == "csr":
                    yield pa.SparseCSRMatrix.from_scipy(scipy_coo.tocsr())
                if format == "csc":
                    yield pa.SparseCSCMatrix.from_scipy(scipy_coo.tocsc())

    def read_table(self, coords: SOMASparseNdCoordinates) -> Iterator[pa.Table]:
        """
        Read a user-defined slice of the sparse array and return in COO format
        as an Arrow Table
        """
        with self._tiledb_open("r") as A:
            sr = clib.SOMAReader(
                self._uri,
                name=self.__class__.__name__,
                schema=A.schema,
            )

            # TODO: make a util function to be shared with SOMADenseNdArray
            # coords are a tuple of (int or slice-of-int)
            if len(coords) == 1 and coords[0] == slice(None):
                # Special case which tiledb-py supports, so we should too
                coords = coords * A.schema.domain.ndim
            elif len(coords) != A.schema.domain.ndim:
                raise ValueError(
                    f"coordinate length {len(coords)} != array ndim {A.schema.domain.ndim}"
                )
            for i in range(A.schema.domain.ndim):
                coord = coords[i]

                if coord is None:
                    continue

                dim_name = A.schema.domain.dim(i).name
                if isinstance(coord, int):
                    sr.set_dim_points(dim_name, [coord])
                elif isinstance(coord, list):
                    sr.set_dim_points(dim_name, coord)
                elif isinstance(coord, pa.ChunkedArray):
                    sr.set_dim_points(dim_name, coord)
                elif isinstance(coord, pa.Array):
                    sr.set_dim_points(dim_name, pa.chunked_array(coord))
                elif isinstance(coord, slice):
                    lo_hi = util.slice_to_range(coord)
                    if lo_hi is not None:
                        sr.set_dim_ranges(dim_name, [lo_hi])
                else:
                    raise ValueError(
                        f"could not handle coordinate with value {coord} of type {type(coord)}"
                    )

            sr.submit()

            # This requires careful handling in the no-data case.
            #
            # When there is at least one table which has non-zero length, we could use the following:
            #   while arrow_table := sr.read_next():
            #     yield arrow_table
            # This respects the caller's expectation that there will always be at least one table in
            # the iterator. (For example, pd.concat(self.read_as_pandas(coords)) will raise an
            # exception if the iterator has zero elements.)
            #
            # But when there is only zero-length data available, the above does *not* respect
            # caller expectations: because a zero-length pyarrow.Table is falsy (not truthy)
            # the while-loop becomes zero-pass.
            #
            # A tempting alternative is to instead write:
            #   for arrow_table in sr.read_next():
            #       yield arrow_table
            # This is correctly one-pass. However, it yields an iterator of pyarrow.ChunkedArray,
            # not an iterator of pyarrow.Table.
            #
            # For this reason, we use the following i > 0 check to guarantee a minimum of one
            # pass through the yielded iterator even when the resulting table is zero-length.
            i = 0
            while True:
                arrow_table = sr.read_next()
                if not arrow_table and i > 0:
                    break
                i += 1
                yield arrow_table

    def read_as_pandas(self, coords: SOMASparseNdCoordinates) -> Iterator[pd.DataFrame]:
        """
        Read a user-defined slice of the sparse array and return as a Pandas DataFrame
        containing COO data.
        """
        for arrow_tbl in self.read_table(coords):
            yield arrow_tbl.to_pandas()

    def read_as_pandas_all(
        self, coords: Optional[SOMASparseNdCoordinates] = None
    ) -> pd.DataFrame:
        """
        Return the sparse array as a single Pandas DataFrame containing COO data.
        """
        if coords is None:
            coords = (slice(None),) * self.ndim
        return pd.concat(self.read_as_pandas(coords))

    def write_sparse_tensor(
        self,
        tensor: Union[
            pa.SparseCOOTensor,
            pa.SparseCSRMatrix,
            pa.SparseCSCMatrix,
        ],
    ) -> None:
        """
        Write an Arrow sparse tensor to the SOMASparseNdArray. The coordinates in the Arrow
        SparseTensor will be interpreted as the coordinates to write to.

        Currently supports the _experimental_ Arrow SparseCOOTensor, SparseCSRMatrix and
        SparseCSCMatrix. There is currently no support for Arrow SparseCSFTensor or dense
        Tensor.
        """
        if isinstance(tensor, pa.SparseCOOTensor):
            data, coords = tensor.to_numpy()
            with self._tiledb_open("w") as A:
                A[tuple(c for c in coords.T)] = data
            return

        if isinstance(tensor, (pa.SparseCSCMatrix, pa.SparseCSRMatrix)):
            if self.ndim != 2:
                raise ValueError(
                    f"Unable to write 2D Arrow sparse matrix to {self.ndim}D SOMASparseNdArray"
                )
            # TODO: the `to_scipy` function is not zero copy. Need to explore zero-copy options.
            sp = tensor.to_scipy().tocoo()
            with self._tiledb_open("w") as A:
                A[sp.row, sp.col] = sp.data
            return

        raise TypeError("Unsuppoted tensor type")

    def write_table(self, arrow_table: pa.Table) -> None:
        """
        Write a COO table, with columns named ``soma_dim_0``, ..., ``soma_dim_N`` and ``soma_data``
        to the dense nD array.
        """
        data = arrow_table.column("soma_data").to_numpy()
        coord_tbl = arrow_table.drop(["soma_data"])
        coords = tuple(
            coord_tbl.column(f"soma_dim_{n}").to_numpy()
            for n in range(coord_tbl.num_columns)
        )
        with self._tiledb_open("w") as A:
            A[coords] = data
