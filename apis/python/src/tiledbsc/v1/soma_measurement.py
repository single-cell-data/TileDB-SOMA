from typing import Any, Optional

import tiledb

from .soma_collection import SOMACollection
from .soma_dataframe import SOMADataFrame
from .tiledb_platform_config import TileDBPlatformConfig


class SOMAMeasurement(SOMACollection):
    """
    The `SOMAMeasurement` is a sub-element of a SOMAExperiment, and is otherwise a specialized
    SOMACollection with pre-defined fields:

    var: SOMADataFrame

    Primary annotations on the _variable_ axis, for variables in this measurement (i.e., annotates
    columns of X). The contents of the soma_rowid pseudo-column define the _variable_ index domain,
    aka varid. All variables for this measurement _must_ be defined in this dataframe.

    X: SOMACollection # of SOMASparseNdArray

    A collection of sparse matrices, each containing measured feature values. Each matrix is indexed
    by [obsid, varid]

    obsm: SOMACollection # of SOMADenseNdArray

    A collection of dense matrices containing annotations of each _obs_ row. Has the same shape as
    obs, and is indexed with obsid.

    obsp: SOMACollection # of SOMASparseNdArray

    A collection of sparse matrices containing pairwise annotations of each _obs_ row. Indexed with
    [obsid_1, obsid_2].

    varm: SOMACollection # of SOMADenseNdArray

    A collection of dense matrices containing annotations of each _var_ row. Has the same shape as
    var, and is indexed with varid

    varp: SOMACollection # of SOMASparseNdArray

    A collection of sparse matrices containing pairwise annotations of each _var_ row. Indexed with
    [varid_1, varid_2]
    """

    _cached_var: Optional[SOMADataFrame]
    _cached_X: Optional[SOMACollection]  # of SOMASparseNdArray
    _cached_obsm: Optional[SOMACollection]  # of SOMADenseNdArray
    _cached_obsp: Optional[SOMACollection]  # of SOMASparseNdArray
    _cached_varm: Optional[SOMACollection]  # of SOMADenseNdArray
    _cached_varp: Optional[SOMACollection]  # of SOMASparseNdArray

    # TODO: check more constraints at runtime:
    # `obs`, `var`
    # Field type is a `SOMADataFrame`

    # `obsp`, `varp`, `X`
    # Field type is a `SOMACollection`, and each element in the collection has a value of type `SOMASparseNdArray`

    # `obsm`, `varm`
    # Field type is a `SOMACollection`, and each element in the collection has a value of type `SOMADenseNdArray`

    # `obsm`, `obsp`, `varm`, `varp`
    # Fields may be empty collections.

    # `X` collection values
    # All matrices must have the shape `(#obs, #var)`. The domain of the first dimension is the values
    # of `obs.soma_rowid`, and the index domain of the second dimension is the values of `var.soma_rowid` in
    # the containing `SOMAMeasurement`.

    # `obsm` collection values
    # All matrices must have the shape `(#obs, M)`, where `M` is user-defined. The domain of the first
    # dimension is the values of `obs.soma_rowid`.

    # `obsp` collection values
    # All matrices must have the shape `(#obs, #obs)`. The domain of both dimensions is the values of
    # `obs.soma_rowid`.

    # `varm` collection values
    # All matrices must have the shape `(#var, M)`, where `M` is user-defined. The domain of the first
    # dimension is the values of `var.soma_rowid`.

    # `varp` collection values
    # All matrices must have the shape `(#var, #var)`. The domain of both dimensions is the values of
    # `var.soma_rowid`.

    def __init__(
        self,
        uri: str,
        *,
        name: Optional[str] = None,
        # Non-top-level objects can have a parent to propagate context, depth, etc.
        parent: Optional[SOMACollection] = None,
        # Top-level objects should specify these:
        tiledb_platform_config: Optional[TileDBPlatformConfig] = None,
        ctx: Optional[tiledb.Ctx] = None,
    ):
        """
        Also see the :class:`TileDBObject` constructor.
        """
        super().__init__(
            uri=uri,
            name=name,
            parent=parent,
            tiledb_platform_config=tiledb_platform_config,
            ctx=ctx,
        )

        self._cached_var = None
        self._cached_X = None
        self._cached_obsm = None
        self._cached_obsp = None
        self._cached_varm = None
        self._cached_varp = None

    def create(self) -> None:
        """
        Creates the data structure on disk/S3/cloud.
        """
        super().create()

    def __getattr__(self, name: str) -> Any:
        """
        TODO: COMMENT
        """
        if name == "var":
            if self._cached_var is None:
                child_uri = self._get_child_uri("var")
                self._cached_var = SOMADataFrame(uri=child_uri, name="var", parent=self)
            return self._cached_var

        elif name == "X":
            if self._cached_X is None:
                child_uri = self._get_child_uri("X")
                self._cached_X = SOMACollection(uri=child_uri, name="X", parent=self)
            return self._cached_X

        elif name == "obsm":
            if self._cached_obsm is None:
                child_uri = self._get_child_uri("obsm")
                self._cached_obsm = SOMACollection(
                    uri=child_uri, name="obsm", parent=self
                )
            return self._cached_obsm

        elif name == "obsp":
            if self._cached_obsp is None:
                child_uri = self._get_child_uri("obsp")
                self._cached_obsp = SOMACollection(
                    uri=child_uri, name="obsp", parent=self
                )
            return self._cached_obsp

        elif name == "varm":
            if self._cached_varm is None:
                child_uri = self._get_child_uri("varm")
                self._cached_varm = SOMACollection(
                    uri=child_uri, name="varm", parent=self
                )
            return self._cached_varm

        elif name == "varp":
            if self._cached_varp is None:
                child_uri = self._get_child_uri("varp")
                self._cached_varp = SOMACollection(
                    uri=child_uri, name="varp", parent=self
                )
            return self._cached_varp

        else:
            # Unlike __getattribute__ this is _only_ called when the member isn't otherwise
            # resolvable. So raising here is the right thing to do.
            raise AttributeError(f"unrecognized attribute: {name}")
