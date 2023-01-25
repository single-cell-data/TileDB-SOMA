from __future__ import annotations

import time
from contextlib import ExitStack
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
    Union,
    cast,
)

import somacore
import tiledb
from typing_extensions import NoReturn

from .exception import DoesNotExistError, SOMAError
from .options import SOMATileDBContext
from .tiledb_object import TileDBObject
from .util import make_relative_path
from .util_tiledb import is_does_not_exist_error, is_duplicate_group_key_error

# A collection can hold any sub-type of TileDBObject
CollectionElementType = TypeVar("CollectionElementType", bound=TileDBObject)


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
        type: Type[Any]
        uri: str
        name: str

        @classmethod
        def from_tdb_object(
            cls: Type[_CachedElement.TdbInfo], o: tiledb.Object
        ) -> _CachedElement._TdbInfo:
            if o.name is None:
                raise SOMAError("expected name to be provided")
            return cls(type=o.type, uri=o.uri, name=o.name)

    tdb: _TdbInfo
    soma: Optional[TileDBObject] = None


class CollectionBase(TileDBObject, somacore.Collection[CollectionElementType]):
    """
    Contains a key-value mapping where the keys are string names and the values
    are any SOMA-defined foundational or composed type, including ``Collection``,
    ``DataFrame``, ``DenseNDArray``, ``SparseNDArray`` or ``Experiment``.
    """

    # Subclass protocol to constrain which SOMA objects types  may be set on a
    # particular collection key. Used by Experiment and Measurement.
    _subclass_constrained_soma_types: Dict[str, Tuple[str, ...]] = {}

    # The collection is persisted as a TileDB Group. The group contents are
    # cached for read performance on higher-latency storage systems such as
    # S3 or TileDB-Cloud. The cache is re-validated on a read miss or on a
    # collection update (add or delete).
    #
    # A value of None implies that the cache has not been loaded, either
    # due to the Collection not existing (i.e., TileDB Group does not exist),
    # or because we have not yet tried to read it.
    #
    _cached_values: Union[Dict[str, _CachedElement], None]

    def __init__(
        self,
        uri: str,
        *,
        context: Optional[SOMATileDBContext] = None,
    ):
        """
        Also see the ``TileDBObject`` constructor.
        """
        super().__init__(uri=uri, context=context)
        self._cached_values = None

    def _not_implemented(self, *args: Any, **kwargs: Any) -> NoReturn:
        raise NotImplementedError()

    add_new_collection = _not_implemented
    add_new_dataframe = _not_implemented
    add_new_dense_ndarray = _not_implemented
    add_new_sparse_ndarray = _not_implemented

    def create_legacy(self) -> "CollectionBase[CollectionElementType]":
        """
        Creates the data structure on disk/S3/cloud.
        """
        return self._create(self.soma_type)

    def _create(self, soma_type: str) -> "CollectionBase[CollectionElementType]":
        """
        Helper for `create`. Ensures that the type name of a child class, not
        its parent class, is written to object-type metadata in storage.
        """
        tiledb.group_create(uri=self._uri, ctx=self.context.tiledb_ctx)
        self._common_create(soma_type)  # object-type metadata etc
        self._cached_values = {}
        return self

    # tiledb.Group handle for [re]use while self is open
    _open_tiledb_group: Optional[tiledb.Group] = None

    def _sub_open(self) -> None:
        assert self._open_mode in ("r", "w") and self._open_tiledb_group is None
        self._open_tiledb_group = self._close_stack.enter_context(
            tiledb.Group(self._uri, mode=self._open_mode, ctx=self.context.tiledb_ctx)
        )

    @property
    def _tiledb_obj(self) -> tiledb.Group:
        "get the open tiledb.Group handle"
        assert self._open_tiledb_group is not None  # => not self.closed
        return self._open_tiledb_group

    def close(self) -> None:
        self._open_tiledb_group = None
        super().close()  # closes self._open_tiledb_group via self._close_stack

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
        err_str = f"{self.__class__.__name__} has no attribute '{key}'"

        # Load cached TileDB Group if not yet loaded or if key not in cache.
        if self._cached_values is None or key not in self._cached_values:
            self._load_tdb_group_cache()
            if self._cached_values is None:
                # This collection was not yet created
                # TODO: SOMA needs better exception types
                raise DoesNotExistError("Collection has not been created")

        # if element is in the TileDB Group, so make a SOMA in-memory object to represent it.
        if key in self._cached_values:
            soma = self._cached_values[key].soma
            if soma is None:
                from .factory import _construct_member

                tdb = self._cached_values[key].tdb
                soma = _construct_member(
                    tdb.uri,
                    context=self.context,
                    object_type=tdb.type.__name__.lower(),
                )
                if soma is None:
                    # if we were unable to create an object, it wasn't actually a SOMA object
                    raise KeyError(err_str)

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

    def _load_tdb_group(self) -> Optional[Dict[str, tiledb.Object]]:
        try:
            with ExitStack() as cleanup:
                cleanup.enter_context(self._ensure_open())
                if self._open_mode == "r":
                    group = self._tiledb_obj
                else:
                    # If self._tiledb_obj is already open in w mode, it won't let us read from
                    # it so we need to open a separate r-mode handle (with no timestamps).
                    # TODO: keep both w- and r-mode handles open if warranted (TBD)
                    group = cleanup.enter_context(
                        tiledb.Group(self._uri, mode="r", ctx=self.context.tiledb_ctx)
                    )
                ans = {o.name: o for o in group if o.name is not None}
            return ans

        except tiledb.TileDBError as e:
            if is_does_not_exist_error(e):
                raise DoesNotExistError("Collection not created") from e
            raise

    def _load_tdb_group_cache(self) -> None:
        """
        Load all objects in the persistent tiledb group. Discard any anonymous objects,
        as all Collection elements must have a key.

        Update the cache with the group contents:
        * delete cached items not in tdb group
        * update/overwrite cache items where there is a URI or type mismatch
        * add any new elements to cache
        """
        tdb_group = self._load_tdb_group()
        if tdb_group is None:
            return

        if self._cached_values is None:
            self._cached_values = {}

        # first, remove anything in cache but not in the group
        for key in list(self._cached_values):
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
                with self._ensure_open("w"):
                    self._tiledb_obj.add(uri=uri, relative=relative, name=key)
                break
            except tiledb.TileDBError as e:
                if is_does_not_exist_error(e):
                    raise DoesNotExistError("Collection not created") from e
                if not is_duplicate_group_key_error(e):
                    raise e
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

        with self._ensure_open("w"):
            # We need to be open "w" here to make sure we're reading our own writes.
            # TODO: keep self open from above, if retry wasn't necessary
            self._load_tdb_group_cache()
        if self._cached_values is not None:
            # Replace the wrapper object synthesized by self.load_tdb_group_cache() with the value
            # passed in to us.
            self._cached_values[key].soma = value
        else:
            # This can only happen if _load_tdb_group_cache failed
            raise SOMAError("Internal error -- unable to save TileDB Group")

    def _del_element(self, key: str, skip_cache_reload: bool = False) -> None:
        if self._cached_values is not None:
            del self._cached_values[key]
        try:
            with self._ensure_open("w"):
                self._tiledb_obj.remove(key)
        except tiledb.TileDBError as e:
            if is_does_not_exist_error(e):
                raise DoesNotExistError("Collection has not been created") from e
            raise

        if not skip_cache_reload:
            self._load_tdb_group_cache()


class Collection(CollectionBase[TileDBObject]):
    """
    A persistent collection of SOMA objects, mapping string keys to any SOMA object.
    """

    pass
    # Inherited from somacore
    # soma_type: Final = "SOMACollection"
