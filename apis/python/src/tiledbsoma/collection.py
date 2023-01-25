from __future__ import annotations

import time
from dataclasses import dataclass
from typing import (
    Any,
    Dict,
    Iterator,
    List,
    Optional,
    Tuple,
    Type,
    TypeVar,
    cast,
)

import somacore
import tiledb
from somacore import options
from typing_extensions import NoReturn

from tiledbsoma.types import StorageType

from .exception import DoesNotExistError, SOMAError
from .options import SOMATileDBContext
from .tiledb_object import AnyTileDBObject, TileDBObject
from .util import make_relative_path
from .util_tiledb import ReadWriteHandle, is_duplicate_group_key_error

# A collection can hold any sub-type of TileDBObject
CollectionElementType = TypeVar("CollectionElementType", bound=AnyTileDBObject)
_Self = TypeVar("_Self", bound="CollectionBase[AnyTileDBObject]")


@dataclass
class _CachedElement:
    """
    Private class representing cached state. Includes the Collection's
    backing tiledb.Group state, which reduces I/O overhead, and
    SOMA objects pointing to the collection elements.
    """

    TdbInfo = TypeVar("TdbInfo", bound="_CachedElement._TdbInfo")

    @dataclass
    class _TdbInfo:
        type: Type[tiledb.Object]
        uri: str
        name: str

        @classmethod
        def from_tdb_object(
            cls: Type[_CachedElement.TdbInfo], o: tiledb.Object
        ) -> _CachedElement._TdbInfo:
            if o.name is None:
                raise SOMAError("expected name to be provided")
            return cls(type=o.type, uri=o.uri, name=o.name)

        def storage_type(self) -> StorageType:
            typename = self.type.__name__.lower()
            assert typename in ("array", "group")
            return cast(StorageType, typename)

    tdb: _TdbInfo
    soma: Optional[AnyTileDBObject] = None


class CollectionBase(
    TileDBObject[tiledb.Group], somacore.Collection[CollectionElementType]
):
    """
    Contains a key-value mapping where the keys are string names and the values
    are any SOMA-defined foundational or composed type, including ``Collection``,
    ``DataFrame``, ``DenseNDArray``, ``SparseNDArray`` or ``Experiment``.
    """

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
        cls._set_create_metadata(handle.writer)
        handle.flush(update_read=True)
        return cls(
            uri, "w", handle, context, _this_is_internal_only="tiledbsoma-internal-code"
        )

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
        return cls(
            uri,
            mode,
            handle,
            context,
            _this_is_internal_only="tiledbsoma-internal-code",
        )

    # Subclass protocol to constrain which SOMA objects types  may be set on a
    # particular collection key. Used by Experiment and Measurement.
    _subclass_constrained_soma_types: Dict[str, Tuple[str, ...]] = {}

    def __init__(
        self,
        uri: str,
        mode: options.OpenMode,
        handle: ReadWriteHandle[tiledb.Group],
        context: SOMATileDBContext,
        *,
        _this_is_internal_only: str = "",
    ):
        super().__init__(
            uri, mode, handle, context, _this_is_internal_only=_this_is_internal_only
        )
        self._cached_values: Optional[Dict[str, _CachedElement]] = None
        """The contents of the persisted TileDB Group, cached for performance.

        The cache is re-validated on a read miss or a collection update (add or
        delete). If None, the cache has not been loaded.
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
        self._load_tdb_group_cache()
        return 0 if self._cached_values is None else len(self._cached_values)

    def __getitem__(self, key: str) -> CollectionElementType:
        """
        Gets the value associated with the key.
        """
        from . import factory  # Delayed binding to resolve circular import.

        err_str = f"{self.__class__.__name__} has no item {key!r}"

        # Load cached TileDB Group if not yet loaded or if key not in cache.
        if self._cached_values is None or key not in self._cached_values:
            self._load_tdb_group_cache()
            if self._cached_values is None:
                raise SOMAError("internal SOMA error loading cached values")

        # if element is in the TileDB Group, so make a SOMA in-memory object to represent it.
        if key in self._cached_values:
            soma = self._cached_values[key].soma
            if soma is None:
                tdb = self._cached_values[key].tdb
                try:
                    soma = factory._open_internal(
                        tdb.uri,
                        self.mode,
                        tiledb_type=tdb.storage_type(),
                        soma_type=None,
                        context=self.context,
                    )
                except (DoesNotExistError, TypeError) as exc:
                    raise KeyError(err_str) from exc

                self._cached_values[key].soma = soma

            return cast(CollectionElementType, self._cached_values[key].soma)

        raise KeyError(err_str)

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
        if self._cached_values is None:
            self._load_tdb_group_cache()
        if self._cached_values is None:
            raise SOMAError("internal error: _cached_values is None")
        return iter(self._cached_values)

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
        if not self.exists():
            return [me]

        self._load_tdb_group_cache()
        if self._cached_values is None:
            raise SOMAError("internal error: _cached_values is None")
        keys = list(self._cached_values.keys())
        me += ":" if len(keys) > 0 else ""
        lines = [me]

        for elmt_key in keys:
            elmt_repr_lines = CollectionBase._get_element_repr((self, elmt_key))
            lines.append(f'  "{elmt_key}": {elmt_repr_lines[0]}')
            for line in elmt_repr_lines[1:]:
                lines.append(f"    {line}")

        return lines

    def _load_tdb_group_cache(self) -> None:
        """
        Load all objects in the persistent tiledb group. Discard any anonymous objects,
        as all Collection elements must have a key.

        Update the cache with the group contents:
        * delete cached items not in tdb group
        * update/overwrite cache items where there is a URI or type mismatch
        * add any new elements to cache
        """
        tdb_group = {o.name: o for o in self._handle.reader if o.name is not None}

        if self._cached_values is None:
            self._cached_values = {}

        # first, remove anything in cache but not in the group
        for key in self._cached_values:
            if key not in tdb_group:
                del self._cached_values[key]

        # then, update cache with the group contents
        for key in tdb_group:
            tdb = tdb_group[key]
            if key in self._cached_values:
                cached = self._cached_values[key]
                if tdb.uri == cached.tdb.uri and tdb.type == cached.tdb.type:
                    continue

            self._cached_values[key] = _CachedElement(
                soma=None, tdb=_CachedElement._TdbInfo.from_tdb_object(tdb)
            )

        if self._cached_values is None:
            raise SOMAError("internal error: _cached_values is None")

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

        uri = (
            make_relative_path(value.uri, relative_to=self.uri)
            if relative
            else value.uri
        )
        for retry in [True, False]:
            try:
                self._handle.writer.add(uri=uri, relative=relative, name=key)
                break
            except tiledb.TileDBError as e:
                if not is_duplicate_group_key_error(e):
                    raise
            if retry:
                self._del_element(key, skip_cache_reload=True)
                # There can be timestamp overlap in a very-rapid-fire unit-test environment.  When
                # that happens, we effectively fall back to filesystem file order, which will be the
                # lexical ordering of the group-metadata filenames. Since the timestamp components
                # are the same, that will be the lexical order of the UUIDs.  So if the new metadata
                # file is sorted before the old one, the group will look like the old state.
                #
                # The standard solution is a negligible but non-zero delay.
                time.sleep(0.001)
                self._handle.flush()
        self._handle.flush(update_read=True)

        # TODO: Timestamp management -- do we need to remain open from above?
        self._load_tdb_group_cache()
        if self._cached_values is None:
            # This can only happen if _load_tdb_group_cache failed
            raise SOMAError("Internal error -- unable to save TileDB Group")
        # Replace the wrapper object synthesized by self.load_tdb_group_cache() with the value
        # passed in to us.
        self._cached_values[key].soma = value

    def _del_element(self, key: str, skip_cache_reload: bool = False) -> None:
        if self._cached_values is not None:
            del self._cached_values[key]
        self._handle.writer.remove(key)
        self._handle.flush(update_read=True)

        if not skip_cache_reload:
            self._load_tdb_group_cache()


AnyTileDBCollection = CollectionBase[AnyTileDBObject]


class Collection(AnyTileDBCollection):
    """
    A persistent collection of SOMA objects, mapping string keys to any SOMA object.

    [lifecycle: experimental]
    """
