# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
from __future__ import annotations

import abc
from collections.abc import MutableMapping, Sequence
from typing import Any, Final, TypeVar, overload

import pyarrow as pa
from typing_extensions import Self

from . import base, data, options

_Elem = TypeVar("_Elem", bound=base.SOMAObject)
"""Element Type for a SOMA collection."""
_CT = TypeVar("_CT", bound="BaseCollection")  # type: ignore[type-arg]
"""Any implementation of a Collection."""


class BaseCollection(base.SOMAObject, MutableMapping[str, _Elem], metaclass=abc.ABCMeta):
    """A generic string-keyed collection of :class:`base.SOMAObject`s.

    The generic type specifies what type the Collection may contain. At its
    most generic, a Collection may contain any SOMA object, but a subclass
    may specify that it is a Collection of a specific type of SOMA object.

    Lifecycle: maturing
    """

    __slots__ = ()

    @classmethod
    @abc.abstractmethod
    def create(
        cls,
        uri: str,
        *,
        platform_config: options.PlatformConfig | None = None,
        context: Any | None = None,  # noqa: ANN401
    ) -> Self:
        """Creates a new collection of this type at the given URI.

        Args:
            uri: The URI where the collection will be created.

        Returns:
            The newly created collection, opened for writing.
        Lifecycle: maturing
        """
        raise NotImplementedError

    @overload
    @abc.abstractmethod
    def add_new_collection(
        self,
        key: str,
        kind: None = None,
        *,
        uri: str | None = ...,
        platform_config: options.PlatformConfig | None = ...,
    ) -> Collection: ...  # type: ignore[type-arg]

    @overload
    @abc.abstractmethod
    def add_new_collection(
        self,
        key: str,
        kind: type[_CT],
        *,
        uri: str | None = ...,
        platform_config: options.PlatformConfig | None = ...,
    ) -> _CT: ...

    @abc.abstractmethod
    def add_new_collection(
        self,
        key: str,
        kind: type[BaseCollection] | None = None,  # type: ignore[type-arg]
        *,
        uri: str | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> BaseCollection:  # type: ignore[type-arg]
        """Creates a new sub-collection of this collection.
        To add an existing collection as a sub-element of this collection,
        use :meth:`set` or indexed access (``coll[name] = value``) instead.

        The type provided is used to create the skeleton of this collection
        as in :meth:`create`. By default, this will create a basic collection::

            # Create a child Measurement object at the key "density"
            # with default settings.
            density = the_collection.add_new_collection("density", somacore.Measurement)

            # This will create a basic Collection as a child at the location
            # storage://authority/path/to/sub-child.
            sub_child = density.add_new_collection(
                "sub_child", uri="storage://authority/path/to/sub-child")

        The URI provided may be absolute or relative. If a child URI is not
        provided, the collection should construct a default child URI based
        on the key of the new entry, making a relative URI (when possible)::

            # coll URI is "file:///path/to/coll"
            coll.add_new_collection("new child!")
            # The URI of the child collection might be:
            # "file:///path/to/coll/new_child"

            # flat_ns_coll URI is "flat://authority/key" in a flat namespace
            # where relative paths are unsupported.
            flat_ns_coll.add_new_collection("flat child")
            # The URI of the child collection might be:
            # "flat://authority/flat_child"

        The way the URI is constructed is left unspecified so that an
        implementation can create a URI based on its own needs. Users should
        directly get the URI of the new child using ``new_child.uri`` if needed;
        they should never assume what it will be.

        Args:
            key: The key that this child should live at
                (i.e., it will be accessed via ``the_collection[key]``).
            kind: The type of child that should be added.
            uri: If provided, overrides the default URI that would be used
                to create this object. This may be absolute or relative.
                If not provided,
            platform_config: Platform-specific configuration options used
                when creating the child.

        Returns:
            The newly created collection, opened for writing.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @abc.abstractmethod
    def add_new_dataframe(
        self,
        key: str,
        *,
        uri: str | None = None,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (options.SOMA_JOINID,),
        domain: Sequence[tuple[Any, Any] | None] | None = None,
        platform_config: options.PlatformConfig | None = None,
    ) -> data.DataFrame:
        """Creates a new DataFrame as a child of this collection.

        Parameters are as in :meth:`data.DataFrame.create`.
        See :meth:`add_new_collection` for details about child URIs.

        Returns:
            The newly created DataFrame, opened for writing.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @abc.abstractmethod
    def add_new_dense_ndarray(
        self,
        key: str,
        *,
        uri: str | None = None,
        type: pa.DataType,
        shape: Sequence[int],
        platform_config: options.PlatformConfig | None = None,
    ) -> data.DenseNDArray:
        """Creates a new dense NDArray as a child of this collection.

        Parameters are as in :meth:`data.DenseNDArray.create`.
        See :meth:`add_new_collection` for details about child URIs.

        Returns:
            The newly created dense NDArray, opened for writing.

        Lifecycle: maturing
        """
        raise NotImplementedError

    @abc.abstractmethod
    def add_new_sparse_ndarray(
        self,
        key: str,
        *,
        uri: str | None = None,
        type: pa.DataType,
        shape: Sequence[int],
        platform_config: options.PlatformConfig | None = None,
    ) -> data.SparseNDArray:
        """Creates a new sparse NDArray as a child of this collection.

        Parameters are as in :meth:`data.SparseNDArray.create`.
        See :meth:`add_new_collection` for details about child creation.

        Returns:
            The newly created sparse NDArray, opened for writing.

        Lifecycle: maturing
        """
        raise NotImplementedError

    def __setitem__(self, key: str, value: _Elem) -> None:
        """Sets an entry into this collection. See :meth:`set` for details."""
        self.set(key, value)

    @abc.abstractmethod
    def set(self, key: str, value: _Elem, *, use_relative_uri: bool | None = None) -> Self:
        """Sets an entry of this collection.

        Important note: Because parent objects may need to share
        implementation-internal state with children, when you set an item in a
        collection, it is not guaranteed that the SOMAObject instance available
        by accessing the collection is the same as the one that was set::

            some_collection["thing"] = my_soma_object
            added_soma_object = some_collection["thing"]
            my_soma_object is added_soma_object  # could be False

        The two objects *will* refer to the same stored data.

        Args:
            key: The string key to set.
            value: The SOMA object to insert into the collection.
            use_relative_uri: Determines whether to store the collection
                entry with a relative URI (provided the storage engine
                supports it).
                If ``None`` (the default), will automatically determine whether
                to use an absolute or relative URI based on their relative
                location.
                If ``True``, will always use a relative URI. If the new child
                does not share a relative URI base, or use of relative URIs
                is not possible at all, the collection should raise an error.
                If ``False``, will always use an absolute URI.

        Returns: ``self``, to enable method chaining.

        Lifecycle: maturing
        """
        raise NotImplementedError


class Collection(BaseCollection[_Elem]):
    """SOMA Collection imposing no semantics on the contained values.

    Lifecycle: maturing
    """

    soma_type: Final = "SOMACollection"  # type: ignore[misc]
    __slots__ = ()
