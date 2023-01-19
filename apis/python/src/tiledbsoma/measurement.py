from typing import Any, Dict, Optional, Tuple, Union, cast

from typing_extensions import Final

from .collection import CollectionBase
from .dataframe import DataFrame
from .dense_nd_array import DenseNDArray
from .options import SOMATileDBContext
from .sparse_nd_array import SparseNDArray
from .tiledb_object import TileDBObject


class Measurement(CollectionBase[TileDBObject]):
    """
    A ``Measurement`` is a sub-element of a ``Experiment``, and is otherwise a specialized ``Collection`` with pre-defined fields:

    ``var``: ``DataFrame``

    Primary annotations on the variable axis, for variables in this measurement (i.e., annotates columns of ``X``). The contents of the ``soma_joinid`` column define the variable index domain, AKA var_id. All variables for this measurement must be defined in this dataframe.

    ``X``: ``Collection`` of ``SparseNDArray``

    A collection of sparse matrices, each containing measured feature values. Each matrix is indexed by ``[obsid, varid]``.

    ``obsm``: ``Collection`` of ``DenseNDArray``

    A collection of dense matrices containing annotations of each ``obs`` row. Has the same shape as ``obs``, and is indexed with ``obsid``.

    ``obsp``: ``Collection`` of ``SparseNDArray``

    A collection of sparse matrices containing pairwise annotations of each ``obs`` row. Indexed with ``[obsid_1, obsid_2]``.

    ``varm``: ``Collection`` of ``DenseNDArray``

    A collection of dense matrices containing annotations of each ``var`` row. Has the same shape as ``var``, and is indexed with ``varid``.

    ``varp``: ``Collection`` of ``SparseNDArray``

    A collection of sparse matrices containing pairwise annotations of each ``var`` row. Indexed with ``[varid_1, varid_2]``
    """

    _subclass_constrained_soma_types: Dict[str, Tuple[str, ...]] = {
        "var": ("SOMADataFrame", "SOMADataFrame"),
        "X": ("SOMACollection",),
        "obsm": ("SOMACollection",),
        "obsp": ("SOMACollection",),
        "varm": ("SOMACollection",),
        "varp": ("SOMACollection",),
    }

    def __init__(
        self,
        uri: str,
        *,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional[CollectionBase[Any]] = None,
        # Top-level objects should specify this:
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(uri=uri, parent=parent, context=context)

    # Inherited from somacore
    soma_type: Final = "SOMAMeasurement"

    def create(self) -> "Measurement":
        """
        Creates the data structure on disk/S3/cloud.
        """
        self._create(self.soma_type)
        return self

    @property
    def var(self) -> DataFrame:
        return cast(DataFrame, self["var"])

    @property
    def X(self) -> CollectionBase[Union[DenseNDArray, SparseNDArray]]:
        return cast(CollectionBase[Union[DenseNDArray, SparseNDArray]], self["X"])

    @property
    def obsm(self) -> CollectionBase[DenseNDArray]:
        return cast(CollectionBase[DenseNDArray], self["obsm"])

    @property
    def obsp(self) -> CollectionBase[SparseNDArray]:
        return cast(CollectionBase[SparseNDArray], self["obsp"])

    @property
    def varm(self) -> CollectionBase[DenseNDArray]:
        return cast(CollectionBase[DenseNDArray], self["varm"])

    @property
    def varp(self) -> CollectionBase[SparseNDArray]:
        return cast(CollectionBase[SparseNDArray], self["varp"])
