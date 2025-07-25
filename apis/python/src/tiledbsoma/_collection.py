# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Implementation of a SOMA Collection."""

from __future__ import annotations

import itertools
from typing import (
    Any,
    Callable,
    ClassVar,
    TypeVar,
    cast,
    overload,
)

import somacore
import somacore.collection
from somacore import options
from typing_extensions import Self

from . import _funcs, _tdb_handles
from . import pytiledbsoma as clib
from ._common_nd_array import NDArray
from ._dataframe import DataFrame
from ._dense_nd_array import DenseNDArray
from ._exception import (
    SOMAError,
    map_exception_for_create,
)
from ._funcs import typeguard_ignore
from ._soma_group import SOMAGroup
from ._soma_object import AnySOMAObject, SOMAObject
from ._sparse_nd_array import SparseNDArray
from ._types import OpenTimestamp
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context

# A collection can hold any sub-type of SOMAObject
CollectionElementType = TypeVar("CollectionElementType", bound=AnySOMAObject)
_TDBO = TypeVar("_TDBO", bound=SOMAObject)  # type: ignore[type-arg]
_Coll = TypeVar("_Coll", bound="CollectionBase[AnySOMAObject]")
_NDArr = TypeVar("_NDArr", bound=NDArray)


class CollectionBase(
    SOMAGroup[CollectionElementType],
    somacore.collection.BaseCollection[CollectionElementType],
):
    """Contains a key-value mapping where the keys are string names and the values
    are any SOMA-defined foundational or composed type, including :class:`Collection`,
    :class:`DataFrame`, :class:`DenseNDArray`, :class:`SparseNDArray` or :class:`Experiment`.
    """

    __slots__ = ()

    # TODO: Implement additional creation of members on collection subclasses.
    @classmethod
    def create(
        cls,
        uri: str,
        *,
        platform_config: options.PlatformConfig | None = None,  # noqa: ARG003
        context: SOMATileDBContext | None = None,
        tiledb_timestamp: OpenTimestamp | None = None,
    ) -> Self:
        """Creates and opens a new SOMA collection in storage.

        This creates a new SOMA collection of the current type in storage and
        returns it opened for writing.

        Args:
            uri:
                The location to create this SOMA collection at.
            platform_config:
                Platform-specific options used to create this collection.
                This may be provided as settings in a dictionary, with options
                located in the ``{'tiledb': {'create': ...}}`` key,
                or as a :class:`~tiledbsoma.TileDBCreateOptions` object.
                (Currently unused for collections.)
            context:
                If provided, the :class:`SOMATileDBContext` to use when creating and
                opening this collection.
            tiledb_timestamp:
                If specified, overrides the default timestamp
                used to open this object. If unset, uses the timestamp provided by
                the context.

        Raises:
            tiledbsoma.AlreadyExistsError:
                If the underlying object already exists at the given URI.
            tiledbsoma.NotCreateableError:
                If the URI is malformed for a particular storage backend.
            TileDBError:
                If unable to create the underlying object.

        Lifecycle:
            Maturing.
        """
        context = _validate_soma_tiledb_context(context)
        try:
            wrapper = cast("_tdb_handles.SOMAGroupWrapper[Any]", cls._wrapper_type)
            timestamp_ms = context._open_timestamp_ms(tiledb_timestamp)
            clib.SOMAGroup.create(
                uri=uri,
                soma_type=wrapper._WRAPPED_TYPE.__name__,
                ctx=context.native_context,
                timestamp=(0, timestamp_ms),
            )
            handle = wrapper.open(uri, "w", context, tiledb_timestamp)
            return cls(
                handle,
                _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
            )
        except SOMAError as e:
            raise map_exception_for_create(e, uri) from None

    # Subclass protocol to constrain which SOMA objects types  may be set on a
    # particular collection key. Used by Experiment and Measurement.
    _subclass_constrained_soma_types: ClassVar[dict[str, tuple[str, ...]]] = {}
    """A map limiting what types may be set to certain keys.

    Map keys are the key of the collection to constrain; values are the SOMA
    type names of the types that may be set to the key.  See :class:`Experiment` and
    :class:`Measurement` for details.
    """

    # Overloads to allow type inference to work when doing:
    #
    #     some_coll.add_new_collection("key")  # -> Collection
    # and
    #     some_coll.add_new_collection("key", Experiment)  # -> Experiment
    #
    # These are only used in type inference to provide better type-checking and
    # autocompletion etc. in static analysis, not at runtime.

    @overload  # type: ignore[override]  # intentionally stricter
    def add_new_collection(
        self,
        key: str,
        kind: None = None,
        *,
        uri: str | None = ...,
        platform_config: options.PlatformConfig | None = ...,
        **kwargs: Any,  # noqa: ANN401
    ) -> "Collection[AnySOMAObject]": ...

    @overload
    def add_new_collection(
        self,
        key: str,
        kind: type[_Coll],
        *,
        uri: str | None = ...,
        platform_config: options.PlatformConfig | None = ...,
        **kwargs: Any,  # noqa: ANN401
    ) -> _Coll: ...

    def add_new_collection(
        self,
        key: str,
        kind: type[CollectionBase] | None = None,  # type: ignore[type-arg]
        *,
        uri: str | None = None,
        platform_config: options.PlatformConfig | None = None,
        **kwargs: Any,
    ) -> AnyTileDBCollection:
        """Adds a new sub-collection to this collection.

        Args:
            key:
                The key to add.
            kind:
                Optionally, the specific type of sub-collection to create.
                For instance, passing ``tiledbsoma.Experiment`` here will create a
                ``SOMAExperiment`` as the sub-entry. By default, a basic
                :class:`Collection` will be created.
            uri:
                If provided, the sub-collection will be created at this URI.
                This can be absolute, in which case the sub-collection will be
                linked to by absolute URI in the stored collection, or relative,
                in which case the sub-collection will be linked to by relative URI.
                The default is to use a relative URI generated based on the key.
                Absolute example: ``uri="s3://mybucket/myexperiment/ms/RNA/newchild"``.
                Relative example: ``uri="newchild"``.
            platform_config:
                Platform configuration options to use when
                creating this sub-collection. This is passed directly to
                ``[CurrentCollectionType].create()``.
            kwargs:
                Other keyword arguments to pass to the ``create`` method of the
                new sub-collection type.

        Examples:
            >>> with tiledbsoma.Collection.create("/tmp/parent") as parent_collection:
            ...     # Create a Collection, with the key ``child_collection``
            ...     parent_collection.add_new_collection("child_collection")
            ...     # And an Experiment, with the key ``child_experiment``
            ...     parent_collection.add_new_collection("child_experiment", tiledbsoma.Experiment)
            ...
            >>> with tiledbsoma.open("/tmp/parent") as parent_collection:
            ...     print(parent_collection['child_collection'].uri)
            ...     print(parent_collection['child_experiment'].uri)
            ...
            file:///tmp/parent/child_collection
            file:///tmp/parent/exp

        Lifecycle:
            Maturing.
        """
        child_cls = kind or Collection
        return self._add_new_element(
            key,
            child_cls,
            lambda create_uri: child_cls.create(
                create_uri,
                platform_config=platform_config,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                **kwargs,
            ),
            uri,
        )

    @_funcs.forwards_kwargs_to(DataFrame.create, exclude=("context", "tiledb_timestamp"))
    def add_new_dataframe(self, key: str, *, uri: str | None = None, **kwargs: Any) -> DataFrame:  # noqa: ANN401
        """Adds a new DataFrame to this collection.

        For details about the behavior of ``key`` and ``uri``, see
        :meth:`add_new_collection`. The remaining parameters are passed to
        :meth:`DataFrame.create` unchanged.

        Examples:
            >>> import tiledbsoma
            >>> import pandas as pd
            >>> import pyarrow as pa
            >>> df = pd.DataFrame(data={"soma_joinid": [0, 1], "col1": [1, 2], "col2": [3, 4]})
            ... with tiledbsoma.Collection.create("/tmp/collection") as soma_collection:
            ...     soma_df = soma_collection.add_new_dataframe(
            ...         "a_dataframe", schema=pa.Schema.from_pandas(df), domain=[[0,9]],
            ...     )
            ...     soma_df.write(pa.Table.from_pandas(df, preserve_index=False))
            ...
            >>> with tiledbsoma.open("/tmp/collection") as soma_collection:
            ...     data = soma_collection['a_dataframe'].read().concat().to_pandas()
            ...
            >>> data
               soma_joinid  col1  col2
            0            0     1     3
            1            1     2     4

        Lifecycle:
            Maturing.
        """
        return self._add_new_element(
            key,
            DataFrame,
            lambda create_uri: DataFrame.create(
                create_uri,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                **kwargs,
            ),
            uri,
        )

    @_funcs.forwards_kwargs_to(NDArray.create, exclude=("context", "tiledb_timestamp"))
    def _add_new_ndarray(self, cls: type[_NDArr], key: str, *, uri: str | None = None, **kwargs: Any) -> _NDArr:  # noqa: ANN401
        """Internal implementation of common NDArray-adding operations."""
        return self._add_new_element(
            key,
            cls,
            lambda create_uri: cls.create(
                create_uri,
                context=self.context,
                tiledb_timestamp=self.tiledb_timestamp_ms,
                **kwargs,
            ),
            uri,
        )

    @_funcs.forwards_kwargs_to(_add_new_ndarray, exclude=("kind",))
    def add_new_dense_ndarray(self, key: str, **kwargs: Any) -> DenseNDArray:  # noqa: ANN401
        """Adds a new DenseNDArray to this Collection.

        For details about the behavior of ``key`` and ``uri``, see
        :meth:`add_new_collection`. The remaining parameters are passed to
        the :meth:`DenseNDArray.create` method unchanged.

        Examples:
            >>> # create a collection and add a (4, 4) dense matrix to it
            >>> with tiledbsoma.Collection.create("./test_collection") as my_collection:
            ...     # collection created. You can now add SOMA objects, e.g., a DenseNDArray
            ...     my_dense_ndarray = my_collection.add_new_dense_ndarray(
            ...         "my_dense_ndarray", type=pa.int32(), shape=(4, 4)
            ...     )
            ...     data = pa.Tensor.from_numpy(np.eye(4, 4, dtype=np.int32))
            ...     my_dense_ndarray.write((slice(None), slice(None)), data)
            ...
            ... # example of opening collection to read an object back
            ... with tiledbsoma.open("./test_collection") as my_collection:
            ...     data = my_collection["my_dense_ndarray"].read()
            ...
            >>> data
            <pyarrow.Tensor>
            type: int32
            shape: (4, 4)
            strides: (16, 4)
            >>> data.to_numpy()
            array([[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]], dtype=int32)

        Lifecycle:
            Maturing.
        """
        return self._add_new_ndarray(DenseNDArray, key, **kwargs)

    @_funcs.forwards_kwargs_to(_add_new_ndarray, exclude=("kind",))
    def add_new_sparse_ndarray(self, key: str, **kwargs: Any) -> SparseNDArray:  # noqa: ANN401
        """Adds a new SparseNDArray to this Collection.

        For details about the behavior of ``key`` and ``uri``, see
        :meth:`add_new_collection`. The remaining parameters are passed to
        the :meth:`SparseNDArray.create` method unchanged.

        Examples:
            >>> with tiledbsoma.Collection.create("./test_collection") as my_collection:
            ...     a_sparse_ndarray = my_collection.add_new_sparse_ndarray(
            ...         "a_sparse_ndarray", type=pa.float32(), shape=(100, 100)
            ...     )
            ...     data = pa.SparseCOOTensor.from_scipy(
            ...         scipy.sparse.random(100, 100, dtype=np.float32)
            ...     )
            ...     a_sparse_ndarray.write(data)
            ...
            >>> with tiledbsoma.open("./test_collection") as my_collection:
            ...     data = my_collection["a_sparse_ndarray"].read().coos().concat()
            ...
            >>> data
            <pyarrow.SparseCOOTensor>
            type: float
            shape: (100, 100)
            >>> data.to_scipy()
            <100x100 sparse matrix of type '<class 'numpy.float32'>'
                    with 100 stored elements in COOrdinate format>

        Lifecycle:
            Maturing.
        """
        return self._add_new_ndarray(SparseNDArray, key, **kwargs)

    def _add_new_element(
        self,
        key: str,
        kind: type[_TDBO],
        factory: Callable[[str], _TDBO],
        user_uri: str | None,
    ) -> _TDBO:
        """Handles the common parts of adding new elements.

        Args:
            key:
                The key to be added.
            kind:
                The type of the element to be added.
            factory:
                A callable that, given the full URI to be added,
                will create the backing storage at that URI and return
                the reified SOMA object.
            user_uri:
                If set, the URI to use for the child
                instead of the default.
        """
        self._check_allows_child(key, kind)
        return super()._add_new_element(key, kind=kind, factory=factory, user_uri=user_uri)

    def members(self) -> dict[str, tuple[str, str]]:
        """Get a mapping of {member_name: (uri, soma_object_type)}."""
        handle = cast("_tdb_handles.SOMAGroupWrapper[Any]", self._handle)
        return handle.members()

    def __repr__(self) -> str:
        """Default display for :class:`Collection`."""
        lines = itertools.chain((self._my_repr(),), self._contents_lines(""))
        return "<" + "\n".join(lines) + ">"

    # ================================================================
    # PRIVATE METHODS FROM HERE ON DOWN
    # ================================================================

    def _my_repr(self) -> str:
        start = super()._my_repr()
        if self.closed:
            return start
        n = len(self)
        if n == 0:
            count = "empty"
        elif n == 1:
            count = "1 item"
        else:
            count = f"{n} items"
        return f"{start} ({count})"

    def _set_element(
        self,
        key: str,
        *,
        uri: str,
        relative: bool,
        soma_object: CollectionElementType,
    ) -> None:
        """Internal implementation of element setting.

        Args:
            key:
                The key to set.
            uri:
                The resolved URI to pass to :meth:`clib.SOMAGroup.add`.
            relative:
                The ``relative`` parameter to pass to ``add``.
            value:
                The reified SOMA object to store locally.
        """
        self._check_allows_child(key, type(soma_object))
        super()._set_element(key, uri=uri, relative=relative, soma_object=soma_object)

    @classmethod
    def _check_allows_child(cls, key: str, child_cls: type) -> None:
        real_child = _real_class(child_cls)
        if not issubclass(real_child, SOMAObject):
            raise TypeError(f"only TileDB objects can be added as children of {cls}, not {child_cls}")
        constraint = cls._subclass_constrained_soma_types.get(key)
        if constraint is not None and real_child.soma_type not in constraint:
            raise TypeError(f"cannot add {child_cls} at {cls}[{key!r}]; only {constraint}")


AnyTileDBCollection = CollectionBase[Any]


class Collection(CollectionBase[CollectionElementType], somacore.Collection[CollectionElementType]):
    """:class:`Collection` is a persistent container of named SOMA objects, stored as
    a mapping of string keys and SOMA object values. Values may be any
    persistent ``tiledbsoma`` object, including :class:`DataFrame`,
    :class:`SparseNDArray`, :class:`DenseNDArray`, :class:`Experiment`, :class:`Measurement`,
    or another :class:`Collection`. A :class:`Collection` refers to elements by a
    per-element URI. A :class:`Collection` may store its reference to an
    element by absolute URI or relative URI.

    Lifecycle:
        Maturing.

    Examples:
        >>> import tiledbsoma
        >>> import pyarrow as pa
        >>> import numpy as np
        >>> # create a collection and add a (4, 4) dense matrix to it
        >>> with tiledbsoma.Collection.create("./test_collection") as my_collection:
        ...     # collection created. You can now add SOMA objects, e.g., a DenseNDArray.
        ...     # New objects are returned open for write.
        ...     my_dense_ndarray = my_collection.add_new_dense_ndarray(
        ...         "my_dense_ndarray", type=pa.int32(), shape=(4, 4)
        ...     )
        ...     data = pa.Tensor.from_numpy(np.eye(4, 4, dtype=np.int32))
        ...     my_dense_ndarray.write((slice(None), slice(None)), data)
        ...
        ... # example of opening collection to read an object back
        ... with tiledbsoma.open("./test_collection") as my_collection:
        ...     data = my_collection["my_dense_ndarray"].read()
        ...
        >>> data
        <pyarrow.Tensor>
        type: int32
        shape: (4, 4)
        strides: (16, 4)
        >>> data.to_numpy()
        array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, 1]], dtype=int32)
    """

    __slots__ = ()

    _wrapper_type = _tdb_handles.CollectionWrapper


@typeguard_ignore
def _real_class(cls: type[Any]) -> type:
    """Extracts the real class from a generic alias.

    Generic aliases like ``Collection[whatever]`` cannot be used in instance or
    subclass checks because they are not actual types present at runtime.
    This extracts the real type from a generic alias::

        _real_class(Collection[whatever])  # -> Collection
        _real_class(List[whatever])  # -> list
    """
    try:
        # If this is a generic alias (e.g. List[x] or list[x]), this will fail.
        issubclass(object, cls)  # Ordering intentional here.
        # Do some extra checking because later Pythons get weird.
        if issubclass(cls, object) and isinstance(cls, type):
            return cls
    except TypeError:
        pass
    err = TypeError(f"{cls} cannot be turned into a real type")
    try:
        # All types of generic alias have this.
        origin = cls.__origin__
        # Other special forms, like Union, also have an __origin__ that is not
        # an actual type.  Verify that the origin is a real, instantiable type.
        issubclass(object, origin)  # Ordering intentional here.
        if issubclass(origin, object) and isinstance(origin, type):
            return origin
    except (AttributeError, TypeError) as exc:
        raise err from exc
    raise err
