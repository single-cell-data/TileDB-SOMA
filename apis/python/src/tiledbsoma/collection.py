from __future__ import annotations

import functools
import re
import time
from typing import (
    Any,
    ClassVar,
    Dict,
    Generic,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    cast,
    overload,
)

import attrs
import pyarrow as pa
import somacore
import tiledb
from somacore import options

from .common_nd_array import NDArray
from .constants import SOMA_JOINID
from .dataframe import DataFrame
from .dense_nd_array import DenseNDArray
from .exception import SOMAError
from .handles import ReadWriteHandle, StorageType, TDBHandle
from .options import SOMATileDBContext
from .sparse_nd_array import SparseNDArray
from .tiledb_object import AnyTileDBObject, TileDBObject
from .util import is_relative_uri, make_relative_path, uri_joinpath
from .util_tiledb import is_does_not_exist_error, is_duplicate_group_key_error

# A collection can hold any sub-type of TileDBObject
CollectionElementType = TypeVar("CollectionElementType", bound=AnyTileDBObject)
_Coll = TypeVar("_Coll", bound="CollectionBase[AnyTileDBObject]")
_Self = TypeVar("_Self", bound="CollectionBase[AnyTileDBObject]")
_NDArr = TypeVar("_NDArr", bound=NDArray)


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

    __slots__ = ()
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
        return cls(
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
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
            handle,
            _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
        )

    # Subclass protocol to constrain which SOMA objects types  may be set on a
    # particular collection key. Used by Experiment and Measurement.
    _subclass_constrained_soma_types: ClassVar[Dict[str, Tuple[str, ...]]] = {}

    def __init__(
        self,
        handle: ReadWriteHandle[tiledb.Group],
        *,
        _dont_call_this_use_create_or_open_instead: str = "",
    ):
        super().__init__(
            handle,
            _dont_call_this_use_create_or_open_instead=_dont_call_this_use_create_or_open_instead,
        )
        self._cached_group_contents: Optional[Dict[str, _CachedElement]] = None
        """The contents of the persisted TileDB Group, cached for performance.

        If None, the cache has not been loaded.
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
        cls: None = None,
        *,
        uri: Optional[str] = ...,
        platform_config: Optional[options.PlatformConfig] = ...,
    ) -> "Collection[AnyTileDBObject]":
        ...

    @overload
    def add_new_collection(
        self,
        key: str,
        cls: Type[_Coll],
        *,
        uri: Optional[str] = ...,
        platform_config: Optional[options.PlatformConfig] = ...,
    ) -> _Coll:
        ...

    def add_new_collection(
        self,
        key: str,
        cls: Optional[Type[AnyTileDBCollection]] = None,
        *,
        uri: Optional[str] = None,
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> "AnyTileDBCollection":
        if key in self:
            raise KeyError(f"{key!r} already exists in {type(self)}")
        child_cls = cls or Collection  # set the default
        self._check_allows_child(key, child_cls)
        absolute_uri, was_relative = self._new_absolute_uri(key=key, uri=uri)
        # mypy gets confused by the inferred union type of child_cls,
        # but this is valid.
        new_group: CollectionBase[Any] = child_cls.create(  # type: ignore[misc]
            absolute_uri, platform_config=platform_config, context=self.context
        )
        self._set_element(
            key, cast(CollectionElementType, new_group), use_relative_uri=was_relative
        )
        # Since we own this object, set it to be closed when we are closed.
        self._close_stack.enter_context(new_group)
        return new_group

    def add_new_dataframe(
        self,
        key: str,
        *,
        uri: Optional[str] = None,
        schema: pa.Schema,
        index_column_names: Sequence[str] = (SOMA_JOINID,),
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> DataFrame:
        if key in self:
            raise KeyError(f"{key!r} already exists in {type(self)}")
        self._check_allows_child(key, DataFrame)
        absolute_uri, was_relative = self._new_absolute_uri(key=key, uri=uri)
        new_df = DataFrame.create(
            absolute_uri,
            index_column_names=index_column_names,
            schema=schema,
            platform_config=platform_config,
            context=self.context,
        )
        # A dataframe might not be the declared type of this collection,
        # but we can't really handle that within the type system.
        self._set_element(
            key, new_df, use_relative_uri=was_relative  # type: ignore[arg-type]
        )
        self._close_stack.enter_context(new_df)
        return new_df

    def _add_new_ndarray(
        self,
        cls: Type[_NDArr],
        key: str,
        *,
        uri: Optional[str] = None,
        type: pa.DataType,
        shape: Sequence[int],
        platform_config: Optional[options.PlatformConfig] = None,
    ) -> _NDArr:
        if key in self:
            raise KeyError(f"{key!r} already exists in {type(self)}")
        self._check_allows_child(key, cls)
        absolute_uri, was_relative = self._new_absolute_uri(key=key, uri=uri)
        new_arr = cls.create(
            absolute_uri,
            type=type,
            shape=shape,
            platform_config=platform_config,
            context=self.context,
        )
        # An NDArray might not be the declared type of this collection,
        # but we can't really handle that within the type system.
        self._set_element(
            key, new_arr, use_relative_uri=was_relative  # type: ignore[arg-type]
        )
        self._close_stack.enter_context(new_arr)
        return new_arr

    # These user-facing functions forward all parameters directly to
    # self._add_new_ndarray, with the specific class type substituted in the
    # first parameter, but without having to duplicate the entire arg list.
    # (mypy doesn't yet understand that these are of the correct type.)
    add_new_dense_ndarray = functools.partialmethod(  # type: ignore[assignment]
        _add_new_ndarray, DenseNDArray
    )
    """Creates a new dense NDArray as a child of this collection.

    Parameters are as in :meth:`DenseNDArray.create`.
    See :meth:`add_new_collection` for details about child creation.
    """
    add_new_sparse_ndarray = functools.partialmethod(  # type: ignore[assignment]
        _add_new_ndarray, SparseNDArray
    )
    """Creates a new sparse NDArray as a child of this collection.

    Parameters are as in :meth:`SparseNDArray.create`.
    See :meth:`add_new_collection` for details about child creation.
    """

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
            # Since we just opened this object, we own it and should close it.
            self._close_stack.enter_context(entry.soma)
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
        self._set_element(key, value, use_relative_uri=use_relative_uri)

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

    def _set_element(
        self,
        key: str,
        value: CollectionElementType,
        use_relative_uri: Optional[bool] = False,
    ) -> None:

        # The SOMA API supports use_relative_uri in [True, False, None].
        # The TileDB-Py API supports use_relative_uri in [True, False].
        # Map from the former to the latter -- and also honor our somacore contract for None --
        # using the following rule.
        if use_relative_uri is None:
            if value.uri.startswith("tiledb://"):
                # TileDB-Cloud does not use relative URIs, ever.
                use_relative_uri = False
            else:
                use_relative_uri = True

        self._check_allows_child(key, type(value))

        # Set has update semantics. Add if missing, delete/add if not. The TileDB Group
        # API only has add/delete. Assume add will succeed, and deal with delete/retry
        # if we get an error on add.

        maybe_relative_uri = (
            make_relative_path(value.uri, relative_to=self.uri)
            if use_relative_uri
            else value.uri
        )
        for retry in [True, False]:
            try:
                self._handle.writer.add(
                    uri=maybe_relative_uri, relative=use_relative_uri, name=key
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

    def _new_absolute_uri(self, *, key: str, uri: Optional[str]) -> Tuple[str, bool]:
        """Builds the full URI for a given child element.

        Returns the full URI, and a boolean indicating whether the URI was
        originally relative.
        """
        maybe_relative_uri = uri or _sanitize_for_path(key)
        if is_relative_uri(maybe_relative_uri):
            return uri_joinpath(self.uri, maybe_relative_uri), True
        return maybe_relative_uri, False

    @classmethod
    def _check_allows_child(cls, key: str, child_cls: type) -> None:
        real_child = _real_class(child_cls)
        if not issubclass(real_child, TileDBObject):
            raise TypeError(
                f"only TileDB objects can be added as children of {cls}, not {child_cls}"
            )
        constraint = cls._subclass_constrained_soma_types.get(key)
        if constraint is not None and real_child.soma_type not in constraint:
            raise TypeError(
                f"cannot add {child_cls} at {cls}[{key!r}]; only {constraint}"
            )


AnyTileDBCollection = CollectionBase[AnyTileDBObject]


class Collection(CollectionBase[CollectionElementType], Generic[CollectionElementType]):
    """
    A persistent collection of SOMA objects, mapping string keys to any SOMA object.

    [lifecycle: experimental]
    """


def _real_class(cls: Type[Any]) -> type:
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
        origin = getattr(cls, "__origin__")
        # Other special forms, like Union, also have an __origin__ that is not
        # an actual type.  Verify that the origin is a real, instantiable type.
        issubclass(object, origin)  # Ordering intentional here.
        if issubclass(origin, object) and isinstance(origin, type):
            return origin
    except (AttributeError, TypeError) as exc:
        raise err from exc
    raise err


_NON_WORDS = re.compile(r"[\W_]+")


def _sanitize_for_path(key: str) -> str:
    """Prepares the given key for use as a path component."""
    sanitized = "_".join(_NON_WORDS.split(key))
    return sanitized
