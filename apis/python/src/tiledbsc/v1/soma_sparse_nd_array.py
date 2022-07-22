from typing import List, Optional, Tuple, Union

import numpy as np
import pyarrow as pa
import tiledb

from .tiledb_array import TileDBArray
from .tiledb_group import TileDBGroup
from .util import tiledb_type_from_arrow_type


class SOMASparseNdArray(TileDBArray):
    """
    Represents ``X`` and others.
    """

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        parent: Optional[TileDBGroup] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """

        super().__init__(uri=uri, name=name, parent=parent)

    def create(
        self,
        type: pa.DataType,
        shape: Union[Tuple, List[int]],
    ) -> None:
        """
        Create a SOMASparseNdArray named with the URI.

        :param type: an Arrow type defining the type of each element in the array. If the type is
        unsupported, an error will be raised.

        :param shape: the length of each domain as a list, e.g., [100, 10]. All lengths must be in
        the uint64 range.
        """

        # checks on shape
        assert len(shape) > 0
        for e in shape:
            assert e > 0

        level = self._tiledb_platform_config.string_dim_zstd_level

        dims = []
        for e in shape:
            dim = tiledb.Dim(
                # Use tiledb default names like `__dim_0`
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
                dtype=tiledb_type_from_arrow_type(type),
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

    #        # ----------------------------------------------------------------
    #        # TODO: type,
    #        shape: Tuple,
    #
    #        # Check that ndims, and each dimension, are positive
    #        assert len(shape) > 0
    #        for e in shape:
    #            assert e > 0

    # TODO: static/class method?
    #    def delete(uri: str) -> None
    #        """
    #        Delete the SOMASparseNdArray specified with the URI.
    #        """

    # TODO: static/class method?
    #    def exists(uri: str) -> bool
    #        """
    #        Return true if object exists and is a SOMASparseNdArray.
    #        """

    #    def get_metadata():
    #        """
    #        Access the metadata as a mutable [`SOMAMetadataMapping`](#SOMAMetadataMapping)
    #        """

    # get_type() is inherited from TileDBObject

    def get_shape(self) -> Tuple:
        """
        Return length of each dimension, always a list of length ``ndims``
        """
        return self._shape

    def get_ndims(self) -> int:
        """
        Return number of index columns
        """
        return len(self._shape)

    #    def get_schema(self) -> Arrow.Schema:
    #        """
    #        Return data schema, in the form of an Arrow Schema
    #        """

    def get_is_sparse(self) -> bool:
        """
        Returns ``True``.
        """
        return True

    #    def get_nnz(self) -> wint:
    #        """
    #        Return the number of non-zero values in the array
    #        """
    #        return 999

    # ----------------------------------------------------------------
    #    def read():
    #        """
    #        Read a slice of data from the SOMASparseNdArray
    #        """

    # ### Operation: read()
    #
    # Read a user-specified subset of the object, and return as one or more Arrow.SparseTensor.
    #
    # Summary:
    #
    # ```
    # read(
    #     [slice, ...],
    #     partitions,
    #     result_order
    # ) -> delayed iterator over Arrow.SparseTensor
    # ```
    #
    # - slice - per-dimension slice, expressed as a scalar, a range, or a list of both.
    # - partitions - an optional [`SOMAReadPartitions`](#SOMAReadPartitions) hint to indicate how
    #   results should be organized.
    # - result_order - order of read results. Can be one of row-major, column-major and unordered.
    #
    # The `read` operation will return a language-specific iterator over one or more Arrow SparseTensor
    # objects, allowing the incremental processing of results larger than available memory. The actual
    # iterator used is delegated to language-specific SOMA specs.

    def write(
        self,
        # TODO: define a compound type for this in types.py
        tensor: Union[pa.SparseCOOTensor, pa.SparseCSFTensor],
    ) -> None:
        """
        Write an Arrow.Tensor to the persistent object. As duplicate index values are not allowed, index
        values already present in the object are overwritten and new index values are added.

        :param values: an Arrow.SparseTensor containing values to be written. The type of elements in `values`
        must match the type of the `SOMASparseNdArray`.
        """

        # Example caller-side:
        #
        # tensor = pa.SparseCOOTensor.from_numpy(
        #     data=np.asarray([4,5,6]),
        #     coords=[[1,2], [3,4], [5,6]],
        #     shape=(nr, nc),
        # )
        #
        # Here on the callee side: need to figure out how to extrat the coords separate from the
        # data since tiledb write needs `A[coords] = data`.

        # >>> tensor.data
        # AttributeError: 'pyarrow.lib.SparseCOOTensor' object has no attribute 'data'

        # >>> tensor.coords
        # AttributeError: 'pyarrow.lib.SparseCOOTensor' object has no attribute 'coords'

        # :headdesk:

        nt = tensor.to_numpy()
        assert len(nt) == 2
        coords = nt[1]
        data = nt[0]

        # The coords come in as a list of [i,j] pairs. We need a list of i's and a list of j's.
        # TODO: for now, only support 2D matrices. We need to complexify this code to handle n-d arrays.
        assert len(coords) > 0
        assert len(coords[0]) == 2
        icoords = coords[:, 0]
        jcoords = coords[:, 1]

        with self._open("w") as A:
            A[icoords, jcoords] = data
