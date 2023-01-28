from __future__ import annotations

import time
from contextlib import ExitStack
from dataclasses import dataclass
from typing import (
    Any,
    Callable,
    Dict,
    Generic,
    Iterator,
    List,
    Optional,
    Tuple,
    Type,
    TypeVar,
    cast,
)

import attrs
import somacore
import tiledb
from somacore import options
from typing_extensions import NoReturn

from .exception import SOMAError
from .options import SOMATileDBContext
from .tiledb_object import AnyTileDBObject, TileDBObject
from .types import StorageType, TDBHandle
from .util import make_relative_path
from .util_tiledb import (
    ReadWriteHandle,
    is_does_not_exist_error,
    is_duplicate_group_key_error,
)

# A collection can hold any sub-type of TileDBObject
CollectionElementType = TypeVar("CollectionElementType", bound=AnyTileDBObject)
_Self = TypeVar("_Self", bound="CollectionBase[AnyTileDBObject]")


@attrs.define()
class _CachedElement:
    """Item we have loaded in the cache of a collection."""

    uri: str
    cls: Type[TDBHandle]
    soma: Optional[AnyTileDBObject] = None
    """The reified object, if it has been opened."""


class CollectionBase(
    TileDBObject[tiledb.Group],
    somacore.Collection[CollectionElementType],
    Generic[CollectionElementType],
):
    """
    Contains a key-value mapping where the keys are string names and the values
    are any SOMA-defined foundational or composed type, including ``Collection``,
    ``DataFrame``, ``DenseNDArray``, ``SparseNDArray`` or ``Experiment``.
    """

    _tiledb_type = tiledb.Group

    # TODO: Implement additional creation of members on collection subclasses.
    @classmethod
    def create(
        cls: Type[_Self],
        uri: str,
        *,
        platform_config: Optional[options.PlatformConfig] = None,
        context: Optional[SOMATileDBContext] = None,
    ) -> _Self:
        context = context or SOMATileDBContext()
        tiledb.group_create(uri=uri, ctx=context.tiledb_ctx)
        handle = ReadWriteHandle.open_group(uri, "w", context)
        cls._set_create_metadata(handle)
        return cls(handle, _this_is_internal_only="tiledbsoma-internal-code")

    @classmethod
    def open(
        cls: Type[_Self],
        uri: str,
        mode: options.OpenMode = "r",
        *,
        context: Optional[SOMATileDBContext] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> _Self:
        """Opens this specific type of collection."""
        del platform_config  # unused
        context = context or SOMATileDBContext()
        handle = ReadWriteHandle.open_group(uri, mode, context)
        # TODO: Verify that we have the right type using metadata.
        return cls(handle, _this_is_internal_only="tiledbsoma-internal-code")

    # Subclass protocol to constrain which SOMA objects types  may be set on a
    # particular collection key. Used by Experiment and Measurement.
    _subclass_constrained_soma_types: Dict[str, Tuple[str, ...]] = {}

    def __init__(
        self,
        handle: ReadWriteHandle[tiledb.Group],
        *,
        _this_is_internal_only: str = "",
    ):
        super().__init__(handle, _this_is_internal_only=_this_is_internal_only)
        self._cached_group_contents: Optional[Dict[str, _CachedElement]] = None
        """The contents of the persisted TileDB Group, cached for performance.

        If None, the cache has not been loaded.
        """

    def _not_implemented(self, *args: Any, **kwargs: Any) -> NoReturn:
        raise NotImplementedError()

    add_new_collection = _not_implemented
    add_new_dataframe = _not_implemented
    add_new_dense_ndarray = _not_implemented
    add_new_sparse_ndarray = _not_implemented

    def __len__(self) -> int:
        """
        Return the number of members in the collection
        """
        return len(self._group_contents)

    def __getitem__(self, key: str) -> CollectionElementType:
        """
        Gets the value associated with the key.
        """

        err_str = f"{self.__class__.__name__} has no item {key!r}"

        try:
            entry = self._group_contents[key]
        except KeyError:
            raise KeyError(err_str) from None
        if entry.soma is None:
            from . import factory  # Delayed binding to resolve circular import.

            storage_type: StorageType
            if issubclass(entry.cls, tiledb.Array):
                storage_type = "array"
            elif issubclass(entry.cls, tiledb.Group):
                storage_type = "group"
            else:
                raise SOMAError(
                    f"internal error: unrecognized group entry type {entry.cls}"
                )
            entry.soma = factory._open_internal(
                uri=entry.uri,
                mode=self.mode,
                context=self.context,
                tiledb_type=storage_type,
                soma_type=None,
            )
        return cast(CollectionElementType, entry.soma)

    def set(
        self,
        key: str,
        value: CollectionElementType,
        *,
        use_relative_uri: Optional[bool] = None,
    ) -> None:
        """
        Adds an element to the collection.  This interface allows explicit control over
        `relative` URI, and uses the member's default name.

        [lifecycle: experimental]
        """
        self._set_element(key, value, relative=use_relative_uri)

    def __setitem__(self, key: str, value: CollectionElementType) -> None:
        """
        Default collection __setattr__
        """
        self._set_element(key, value)

    def __delitem__(self, key: str) -> None:
        """
        Removes a member from the collection, when invoked as ``del collection["namegoeshere"]``.
        """
        self._del_element(key)

    def __iter__(self) -> Iterator[str]:
        return iter(self._group_contents)

    def __repr__(self) -> str:
        """
        Default display for ``Collection``.
        """
        return "\n".join(self._get_collection_repr())

    # ================================================================
    # PRIVATE METHODS FROM HERE ON DOWN
    # ================================================================

    @classmethod
    def _get_element_repr(
        cls, args: Tuple[CollectionBase[CollectionElementType], str]
    ) -> List[str]:
        collection, key = args
        value = collection.__getitem__(key)
        if isinstance(value, CollectionBase):
            return value._get_collection_repr()
        else:
            return [value.__repr__()]

    def _get_collection_repr(self) -> List[str]:
        me = super().__repr__()
        keys = list(self._group_contents.keys())
        me += ":" if len(keys) > 0 else ""
        lines = [me]

        for elmt_key in keys:
            elmt_repr_lines = CollectionBase._get_element_repr((self, elmt_key))
            lines.append(f'  "{elmt_key}": {elmt_repr_lines[0]}')
            for line in elmt_repr_lines[1:]:
                lines.append(f"    {line}")

        return lines

    @property
    def _group_contents(self) -> Dict[str, _CachedElement]:
        """
        Load all objects in the persistent tiledb group. Discard any anonymous objects,
        as all Collection elements must have a key.

        Update the cache with the group contents:
        * delete cached items not in tdb group
        * update/overwrite cache items where there is a URI or type mismatch
        * add any new elements to cache
        """
        if self._cached_group_contents is None:
            self._cached_group_contents = {
                o.name: _CachedElement(uri=o.uri, cls=o.type)
                for o in self._handle.reader
                if o.name is not None
            }
        return self._cached_group_contents

    def _determine_default_relative(self, uri: str) -> Optional[bool]:
        """Defaulting for the relative parameter."""
        if self.context.member_uris_are_relative is not None:
            return self.context.member_uris_are_relative
        if uri.startswith("tiledb://"):
            # TileDB-Cloud does not use relative URIs, ever.
            return False
        return None

    def _set_element(
        self, key: str, value: CollectionElementType, relative: Optional[bool] = None
    ) -> None:
        if relative is None:
            relative = self._determine_default_relative(value.uri)

        if key in self._subclass_constrained_soma_types:
            # Implement the sub-class protocol constraining the value type of certain item keys
            accepted_types = self._subclass_constrained_soma_types[key]
            if (
                not isinstance(value, TileDBObject)
                or value.soma_type not in accepted_types
            ):
                NL = "\n"
                raise TypeError(
                    f"{self.__class__.__name__} field '{key}' only accepts values of type(s) {NL.join(accepted_types)}."
                )

        # Set has update semantics. Add if missing, delete/add if not. The TileDB Group
        # API only has add/delete. Assume add will succeed, and deal with delete/retry
        # if we get an error on add.

        maybe_relative_uri = (
            make_relative_path(value.uri, relative_to=self.uri)
            if relative
            else value.uri
        )
        for retry in [True, False]:
            try:
                self._handle.writer.add(
                    uri=maybe_relative_uri, relative=relative, name=key
                )
                break
            except tiledb.TileDBError as e:
                if not is_duplicate_group_key_error(e):
                    raise
            if retry:
                self._del_element(key)

                # There can be timestamp overlap in a very-rapid-fire unit-test environment.  When
                # that happens, we effectively fall back to filesystem file order, which will be the
                # lexical ordering of the group-metadata filenames. Since the timestamp components
                # are the same, that will be the lexical order of the UUIDs.  So if the new metadata
                # file is sorted before the old one, the group will look like the old state.
                #
                # The standard solution is a negligible but non-zero delay.
                time.sleep(0.001)
        # HACK: There is no way to change a group entry without deleting it and
        # re-adding it, but you can't do both of those in the same transaction.
        # You get "member already set for removal" in an error.
        #
        # This also means that if, in one transaction, you do
        #     grp["x"] = y
        #     del grp["x"]
        # you would also get an error without this hack.
        self._handle._flush_hack()

        self._group_contents[key] = _CachedElement(
            uri=value.uri, cls=value._tiledb_type, soma=value
        )

    def _del_element(self, key: str) -> None:
        try:
            self._handle.writer.remove(key)
            # HACK: see note above
            self._handle._flush_hack()
            self._group_contents.pop(key, None)
        except tiledb.TileDBError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise KeyError(f"{key!r} does not exist in {self}") from tdbe
            raise


AnyTileDBCollection = CollectionBase[AnyTileDBObject]


class Collection(CollectionBase[CollectionElementType], Generic[CollectionElementType]):
    """
    A persistent collection of SOMA objects, mapping string keys to any SOMA object.

    [lifecycle: experimental]
    """
