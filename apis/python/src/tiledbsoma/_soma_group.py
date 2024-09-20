# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc.
#
# Licensed under the MIT License.

import re
from typing import (
    Any,
    Callable,
    Generic,
    Iterable,
    Iterator,
    Optional,
    Set,
    Type,
    TypeVar,
    cast,
)

import attrs
from typing_extensions import Self

from . import _tdb_handles

# This package's pybind11 code
from . import pytiledbsoma as clib  # noqa: E402
from ._exception import (
    SOMAError,
    is_does_not_exist_error,
)
from ._soma_object import AnySOMAObject, SOMAObject
from ._util import (
    is_relative_uri,
    make_relative_path,
    uri_joinpath,
)

CollectionElementType = TypeVar("CollectionElementType", bound=AnySOMAObject)
_TDBO = TypeVar("_TDBO", bound=SOMAObject)  # type: ignore[type-arg]


@attrs.define()
class _CachedElement:
    """Item we have loaded in the cache of a collection."""

    entry: _tdb_handles.GroupEntry
    soma: Optional[AnySOMAObject] = None
    """The reified object, if it has been opened."""


class SOMAGroup(
    SOMAObject[_tdb_handles.SOMAGroupWrapper[Any]], Generic[CollectionElementType]
):
    """Base class for all SOMAGroups: CollectionBase and MultiscaleImage.

    Lifecycle:
        Experimental.
    """

    __slots__ = ("_contents", "_mutated_keys")

    def __init__(
        self,
        handle: _tdb_handles.SOMAGroupWrapper[Any],
        **kwargs: Any,
    ):
        super().__init__(handle, **kwargs)
        self._contents = {
            key: _CachedElement(entry) for key, entry in handle.initial_contents.items()
        }
        """The contents of the persisted TileDB Group.

        This is loaded at startup when we have a read handle.
        """
        self._mutated_keys: Set[str] = set()

    def __len__(self) -> int:
        """Return the number of members in the collection"""
        return len(self._contents)

    def __getitem__(self, key: str) -> CollectionElementType:
        """Gets the value associated with the key."""
        err_str = f"{self.__class__.__name__} has no item {key!r}"

        try:
            entry = self._contents[key]
        except KeyError:
            raise KeyError(err_str) from None
        if entry.soma is None:
            from . import _factory  # Delayed binding to resolve circular import.

            uri = entry.entry.uri
            mode = self.mode
            context = self.context
            timestamp = self.tiledb_timestamp_ms
            clib_type = entry.entry.wrapper_type.clib_type

            wrapper = _tdb_handles.open(uri, mode, context, timestamp, clib_type)
            entry.soma = _factory.reify_handle(wrapper)

            # Since we just opened this object, we own it and should close it.
            self._close_stack.enter_context(entry.soma)
        return cast(CollectionElementType, entry.soma)

    def __setitem__(self, key: str, value: CollectionElementType) -> None:
        """Default collection __setattr__"""
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
                yield f"{indent}{key!r}: {entry.entry.uri!r} (unopened)"
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
        clib_collection = self._handle._handle
        relative_type = clib.URIType.relative if relative else clib.URIType.absolute
        clib_collection.add(
            uri=uri,
            uri_type=relative_type,
            name=key,
            soma_type=clib_collection.type,
        )
        self._contents[key] = _CachedElement(
            entry=_tdb_handles.GroupEntry(soma_object.uri, soma_object._wrapper_type),
            soma=soma_object,
        )
        self._mutated_keys.add(key)

    def _del_element(self, key: str) -> None:
        if key in self._mutated_keys:
            raise SOMAError(f"cannot delete previously-mutated key {key!r}")
        try:
            self._handle.writer.remove(key)
        except RuntimeError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise KeyError(tdbe) from tdbe
            raise
        self._contents.pop(key, None)
        self._mutated_keys.add(key)

    def _add_new_element(
        self,
        key: str,
        kind: Type[_TDBO],
        factory: Callable[[str], _TDBO],
        user_uri: Optional[str],
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

    def _new_child_uri(self, *, key: str, user_uri: Optional[str]) -> "_ChildURI":
        maybe_relative_uri = user_uri or _sanitize_for_path(key)
        if not is_relative_uri(maybe_relative_uri):
            # It's an absolute URI.
            return _ChildURI(
                add_uri=maybe_relative_uri,
                full_uri=maybe_relative_uri,
                relative=False,
            )
        if not self.uri.startswith("tiledb://"):
            # We don't need to post-process anything.
            return _ChildURI(
                add_uri=maybe_relative_uri,
                full_uri=uri_joinpath(self.uri, maybe_relative_uri),
                relative=True,
            )
        # Our own URI is a ``tiledb://`` URI. Since TileDB Cloud requires absolute
        # URIs, we need to calculate the absolute URI to pass to Group.add
        # based on our creation URI.
        # TODO: Handle the case where we reopen a TileDB Cloud Group, but by
        # name rather than creation path.
        absolute_uri = uri_joinpath(self.uri, maybe_relative_uri)
        return _ChildURI(add_uri=absolute_uri, full_uri=absolute_uri, relative=False)

    def set(
        self,
        key: str,
        value: CollectionElementType,
        *,
        use_relative_uri: Optional[bool] = None,
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

        self._set_element(
            key, uri=uri_to_add, relative=use_relative_uri, soma_object=value
        )
        return self


_NON_WORDS = re.compile(r"[\W_]+")


def _sanitize_for_path(key: str) -> str:
    """Prepares the given key for use as a path component."""
    sanitized = "_".join(_NON_WORDS.split(key))
    return sanitized


@attrs.define(frozen=True, kw_only=True)
class _ChildURI:
    add_uri: str
    """The URI of the child for passing to :meth:``clib.SOMAGroup.add``."""
    full_uri: str
    """The full URI of the child, used to create a new element."""
    relative: bool
    """The ``relative`` value to pass to :meth:``clib.SOMAGroup.add``."""
