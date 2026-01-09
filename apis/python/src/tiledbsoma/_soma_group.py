# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import warnings
from collections.abc import Iterable, Iterator
from threading import Lock
from typing import Any, Callable, Generic, TypeVar, cast

import attrs
from somacore import options
from typing_extensions import Self

from . import _tdb_handles

# This package's pybind11 code
from . import pytiledbsoma as clib
from ._exception import SOMAError, UnsupportedOperationError, is_does_not_exist_error
from ._soma_context import SOMAContext
from ._soma_object import SOMAObject
from ._types import OpenTimestamp, SOMABaseTileDBType
from ._util import is_relative_uri, make_relative_path, sanitize_key, uri_joinpath

CollectionElementType = TypeVar("CollectionElementType", bound=SOMAObject)
_TDBO = TypeVar("_TDBO", bound=SOMAObject)


@attrs.define()
class _CachedElement:
    """Item we have loaded in the cache of a collection."""

    uri: str
    tiledb_type: SOMABaseTileDBType | None = None
    soma: SOMAObject | None = None
    """The reified object, if it has been opened."""

    @classmethod
    def from_handle_entry(cls, obj: tuple[str, str]) -> _CachedElement:
        uri, type = obj[0], obj[1]
        if type == "SOMAArray":
            return _CachedElement(uri, SOMABaseTileDBType.SOMAArray)
        if type == "SOMAGroup":
            return _CachedElement(uri, SOMABaseTileDBType.SOMAGroup)
        raise SOMAError(f"internal error: unknown object type {uri}")


class SOMAGroup(SOMAObject, Generic[CollectionElementType]):
    """Base class for all SOMAGroups: CollectionBase and MultiscaleImage.

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_contents", "_mutated_keys", "_reify_lock")

    def __init__(
        self,
        handle: _tdb_handles.RawHandle,
        *,
        uri: str,
        context: SOMAContext,
        **kwargs: Any,  # noqa: ANN401
    ) -> None:
        super().__init__(handle, uri=uri, context=context, **kwargs)
        self._contents = {
            name: _CachedElement.from_handle_entry(entry) for name, entry in self._handle.members().items()
        }
        """The contents of the persisted TileDB Group.

        This is loaded at startup when we have a read handle.
        """
        self._mutated_keys: set[str] = set()
        self._reify_lock = Lock()

    def __len__(self) -> int:
        """Return the number of members in the collection."""
        return len(self._contents)

    def __getitem__(self, key: str) -> CollectionElementType:
        """Gets the value associated with the key."""
        err_str = f"{self.__class__.__name__} has no item {key!r}"

        try:
            entry = self._contents[key]
        except KeyError:
            raise KeyError(err_str) from None

        with self._reify_lock:
            if entry.soma is None:
                from . import _factory  # Delayed binding to resolve circular import.

                entry.soma = _factory._open_soma_object(
                    uri=entry.uri,
                    mode=self.mode,
                    context=self.context,
                    tiledb_timestamp=self.tiledb_timestamp_ms,
                    clib_type=None if entry.tiledb_type is None else entry.tiledb_type.name,
                )

                # Since we just opened this object, we own it and should close it.
                self._close_stack.enter_context(entry.soma)

        return cast("CollectionElementType", entry.soma)

    def __setitem__(self, key: str, value: CollectionElementType) -> None:
        """Default collection __setattr__."""
        self.set(key, value, use_relative_uri=None)

    def __delitem__(self, key: str) -> None:
        """Removes a member from the collection, when invoked as ``del collection["namegoeshere"]``.

        Raises:
            SOMAError:
                Upon deletion of a mutated key.
        """
        self._del_element(key)

    def __iter__(self) -> Iterator[str]:
        return iter(self._contents)

    def _contents_lines(self, last_indent: str) -> Iterable[str]:
        indent = last_indent + "    "
        if self.closed:
            return
        for key, entry in self._contents.items():
            obj = entry.soma
            if obj is None:
                # We haven't reified this SOMA object yet. Don't try to open it.
                yield f"{indent}{key!r}: {entry.uri!r} (unopened)"
            else:
                yield f"{indent}{key!r}: {obj._my_repr()}"
                if isinstance(obj, SOMAGroup):
                    yield from obj._contents_lines(indent)

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
        if key in self._mutated_keys.union(self._contents):
            # TileDB groups currently do not support replacing elements.
            # If we use a hack to flush writes, corruption is possible.
            raise SOMAError(f"replacing key {key!r} is unsupported")
        clib_collection = self._handle
        relative_type = clib.URIType.relative if relative else clib.URIType.absolute
        if self.context.is_tiledbv2_uri(self.uri):
            clib_collection.add(
                uri=uri,
                uri_type=relative_type,
                name=key,
                soma_type=soma_object.soma_type,
            )
        self._contents[key] = _CachedElement(uri=soma_object.uri, tiledb_type=None, soma=soma_object)
        self._mutated_keys.add(key)

    def _del_element(self, key: str) -> None:
        if key in self._mutated_keys:
            raise SOMAError(f"Cannot delete previously-mutated key '{key!r}'.")
        try:
            if self.closed:
                raise SOMAError(f"Cannot delete '{key!r}'. {self} is closed")
            if self.mode == "d":
                self._handle.remove(key)
            elif self.mode == "w":
                warnings.warn(
                    f"Deleting in write mode is deprecated. {self} should be reopened with mode='d'.",
                    DeprecationWarning,
                    stacklevel=3,
                )
                self._handle.remove(key)
            else:
                raise SOMAError(
                    f"Deleting is not allowed in mode '{self.mode}'. {self} should be reopened with mode='d'."
                )
        except Exception as tdbe:
            if is_does_not_exist_error(tdbe):
                raise KeyError(tdbe) from tdbe
            raise
        self._contents.pop(key, None)
        self._mutated_keys.add(key)

    def _add_new_element(
        self,
        key: str,
        kind: type[_TDBO],  # noqa: ARG002
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
        if key in self:
            raise KeyError(f"{key!r} already exists in {type(self)}")
        child_uri = self._new_child_uri(key=key, user_uri=user_uri)

        if self.context.is_tiledbv3_uri(self.uri) and (not child_uri.relative or key != child_uri.add_uri):
            # Carrara data model requires relative URI, and that member name == member uri
            raise UnsupportedOperationError(
                "TileDB Carrara data model requires Collection member name and uri to be equal."
            )

        child = factory(child_uri.full_uri)
        # The resulting element may not be the right type for this collection,
        # but we can't really handle that within the type system.
        self._set_element(
            key,
            uri=child_uri.add_uri,
            relative=child_uri.relative,
            soma_object=child,  # type: ignore[arg-type]
        )
        self._close_stack.enter_context(child)
        return child

    def _new_child_uri(self, *, key: str, user_uri: str | None) -> _ChildURI:
        maybe_relative_uri = user_uri or sanitize_key(key, data_protocol=self.context.data_protocol(self.uri))
        if not is_relative_uri(maybe_relative_uri):
            # It's an absolute URI.
            return _ChildURI(
                add_uri=maybe_relative_uri,
                full_uri=maybe_relative_uri,
                relative=False,
            )
        if not self.uri.startswith("tiledb://") or self.context.is_tiledbv3_uri(self.uri):
            # We don't need to post-process anything - URI schema handles relative paths.
            return _ChildURI(
                add_uri=maybe_relative_uri,
                full_uri=uri_joinpath(self.uri, maybe_relative_uri),
                relative=True,
            )
        # TileDB Cloud requires absolute URIs; we need to calculate the absolute URI to pass to Group.add
        # based on our creation URI.
        # TODO: Handle the case where we reopen a TileDB Cloud Group, but by name rather than creation path.
        absolute_uri = uri_joinpath(self.uri, maybe_relative_uri)
        return _ChildURI(add_uri=absolute_uri, full_uri=absolute_uri, relative=False)

    def reopen(self, mode: options.OpenMode, tiledb_timestamp: OpenTimestamp | None = None) -> Self:
        """Return a new copy of the SOMAObject with the given mode at the current
        Unix timestamp.

        Args:
            mode:
                The mode to open the object in.
                - ``r``: Open for reading only (cannot write or delete).
                - ``w``: Open for writing only (cannot read or delete).
                - ``d``: Open for deleting only (cannot read or write).
            tiledb_timestamp:
                The TileDB timestamp to open this object at, either an int representing milliseconds since the Unix
                epoch or a datetime.datetime object. When not provided (the default), the current time is used. A
                value of zero results in current time.

        Raises:
            ValueError:
                If the user-provided ``mode`` is invalid.
            SOMAError:
                If the object has unwritten metadata.

        Lifecycle:
            Experimental.
        """
        super().reopen(mode, tiledb_timestamp)
        self._contents = {
            name: _CachedElement.from_handle_entry(entry) for name, entry in self._handle.members().items()
        }
        self._mutated_keys = set()
        return self

    def set(
        self,
        key: str,
        value: CollectionElementType,
        *,
        use_relative_uri: bool | None = None,
    ) -> Self:
        """Adds an element to the collection.

        Args:
            key:
                The key of the element to be added.
            value:
                The value to be added to this collection.
            use_relative_uri:
                By default (None), the collection will determine whether the
                element should be stored by relative URI.
                If True, the collection will store the child by absolute URI.
                If False, the collection will store the child by relative URI.

        Raises:
            SOMAError:
                If an existing key is set (replacement is unsupported).

        Lifecycle:
            Maturing.
        """
        uri_to_add = value.uri
        if self.context.is_tiledbv3_uri(uri_to_add):
            raise UnsupportedOperationError(
                "TileDB Carrara data model does not support the set or __setitem__ operation on this object."
            )

        # The SOMA API supports use_relative_uri in [True, False, None].
        # The TileDB-Py API supports use_relative_uri in [True, False].
        # Map from the former to the latter -- and also honor our somacore contract for None --
        # using the following rule.
        if use_relative_uri is None and value.uri.startswith("tiledb://"):
            # TileDB-Cloud does not use relative URIs, ever.
            use_relative_uri = False

        if use_relative_uri is not False:
            try:
                uri_to_add = make_relative_path(value.uri, relative_to=self.uri)
                use_relative_uri = True
            except ValueError:
                if use_relative_uri:
                    # We couldn't construct a relative URI, but we were asked
                    # to use one, so raise the error.
                    raise
                use_relative_uri = False

        self._set_element(key, uri=uri_to_add, relative=use_relative_uri, soma_object=value)
        return self


@attrs.define(frozen=True, kw_only=True)
class _ChildURI:
    add_uri: str
    """The URI of the child for passing to :meth:``clib.SOMAGroup.add``."""
    full_uri: str
    """The full URI of the child, used to create a new element."""
    relative: bool
    """The ``relative`` value to pass to :meth:``clib.SOMAGroup.add``."""
