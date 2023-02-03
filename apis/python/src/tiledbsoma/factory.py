"""
This module exists to avoid what would otherwise be cyclic-module-import issues within
Collection.
"""

from typing import Dict, Optional, Type, TypeVar, cast, overload

import tiledb
from somacore import options

from . import (
    collection,
    dataframe,
    dense_nd_array,
    experiment,
    handles,
    measurement,
    sparse_nd_array,
    tiledb_array,
    tiledb_object,
)
from .constants import (
    SOMA_ENCODING_VERSION,
    SOMA_ENCODING_VERSION_METADATA_KEY,
    SOMA_OBJECT_TYPE_METADATA_KEY,
)
from .exception import DoesNotExistError, SOMAError
from .options import SOMATileDBContext
from .util import typeguard_ignore

_Obj = TypeVar("_Obj", bound="tiledb_object.AnyTileDBObject")
_Arr = TypeVar("_Arr", bound="tiledb_array.TileDBArray")
_Coll = TypeVar("_Coll", bound="collection.AnyTileDBCollection")


@overload
def open(
    uri: str,
    mode: options.OpenMode = ...,
    *,
    soma_type: None = None,  # TODO: Accept strings here.
    context: Optional[SOMATileDBContext] = None,
) -> "tiledb_object.AnyTileDBObject":
    ...


@overload
def open(
    uri: str,
    mode: options.OpenMode,
    *,
    soma_type: Type[_Obj] = ...,  # TODO: Accept strings here.
    context: Optional[SOMATileDBContext] = None,
) -> _Obj:
    ...


@typeguard_ignore
def open(
    uri: str,
    mode: options.OpenMode = "r",
    *,
    soma_type: Optional[
        Type["tiledb_object.AnyTileDBObject"]
    ] = None,  # TODO: Accept strings here.
    context: Optional[SOMATileDBContext] = None,
) -> "tiledb_object.AnyTileDBObject":
    context = context or SOMATileDBContext()
    return _open_internal(
        uri, mode, tiledb_type=None, soma_type=soma_type, context=context
    )


@overload
def _open_internal(
    uri: str,
    mode: options.OpenMode,
    *,
    tiledb_type: Optional[handles.StorageType],
    soma_type: None,
    context: SOMATileDBContext,
) -> "tiledb_object.AnyTileDBObject":
    ...


@overload
def _open_internal(
    uri: str,
    mode: options.OpenMode,
    *,
    tiledb_type: Optional[handles.StorageType],
    soma_type: Type[_Obj],
    context: SOMATileDBContext,
) -> _Obj:
    ...


@typeguard_ignore
def _open_internal(
    uri: str,
    mode: options.OpenMode,
    *,
    tiledb_type: Optional[handles.StorageType],
    soma_type: Optional[Type["tiledb_object.AnyTileDBObject"]],
    context: SOMATileDBContext,
) -> "tiledb_object.AnyTileDBObject":
    if soma_type:
        soma_tdb_type = _storage_of(soma_type)
        if tiledb_type and (tiledb_type != soma_tdb_type):
            raise TypeError(
                f"SOMA type {soma_type.soma_type} is stored as {soma_tdb_type},"
                f" but {tiledb_type} was requested"
            )
        tiledb_type = soma_tdb_type
    elif not tiledb_type:
        obj_type = tiledb.object_type(uri, ctx=context.tiledb_ctx)
        if not obj_type:
            raise DoesNotExistError(f"{uri!r} does not exist")
        tiledb_type = cast(handles.StorageType, obj_type)

    if tiledb_type == "array":
        if not soma_type or issubclass(soma_type, tiledb_array.TileDBArray):
            return _open_array(uri, mode, soma_type, context)  # type: ignore[type-var]

        raise TypeError(f"stored array object cannot be opened as a {soma_type}")

    if tiledb_type == "group":
        if not soma_type or issubclass(soma_type, collection.CollectionBase):
            return _open_group(uri, mode, soma_type, context)  # type: ignore[type-var]
        raise TypeError(f"stored group object cannot be opened as a {soma_type}")

    raise SOMAError(f"unknown stored TileDB object type {tiledb_type!r}")


@typeguard_ignore
def _open_array(
    uri: str,
    mode: options.OpenMode,
    soma_type: Optional[Type[_Arr]],
    context: SOMATileDBContext,
) -> _Arr:
    handle = handles.ArrayWrapper.open(uri, mode, context)
    try:
        data_type = _read_soma_type(handle)
        cls = _to_array_class(data_type)
        want_type = soma_type or cls
        if cls != want_type:
            raise TypeError(f"stored data is {cls}; expected {want_type}")
        return cast(
            _Arr,
            cls(
                handle,
                _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
            ),
        )
    except Exception:
        handle.close()
        raise


@typeguard_ignore
def _open_group(
    uri: str,
    mode: options.OpenMode,
    soma_type: Optional[Type[_Coll]],
    context: SOMATileDBContext,
) -> _Coll:
    handle = handles.GroupWrapper.open(uri, mode, context)
    try:
        data_type = _read_soma_type(handle)
        cls = _to_group_class(data_type)
        want_type = soma_type or cls
        if cls != want_type:
            raise TypeError(f"stored data is {cls}; expected {want_type}")
        return cast(
            _Coll,
            cls(
                handle,
                _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code",
            ),
        )
    except Exception:
        handle.close()
        raise


@typeguard_ignore
def _storage_of(cls: Type["tiledb_object.AnyTileDBObject"]) -> handles.StorageType:
    if issubclass(cls, tiledb_array.TileDBArray):
        return "array"
    if issubclass(cls, collection.CollectionBase):
        return "group"
    raise TypeError(f"{cls} is not a concrete stored object")


def _read_soma_type(hdl: handles.AnyWrapper) -> str:
    obj_type = hdl.metadata.get(SOMA_OBJECT_TYPE_METADATA_KEY)
    encoding_version = hdl.metadata.get(SOMA_ENCODING_VERSION_METADATA_KEY)

    if obj_type is None:
        raise SOMAError(
            f"stored TileDB object does not have {SOMA_OBJECT_TYPE_METADATA_KEY!r}"
        )
    if not isinstance(obj_type, str):
        raise SOMAError(
            f"stored TileDB object {SOMA_OBJECT_TYPE_METADATA_KEY!r}"
            f" is {type(obj_type)}"
        )
    if encoding_version is None:
        raise SOMAError(
            f"stored TileDB object does not have {SOMA_ENCODING_VERSION_METADATA_KEY}"
        )
    if encoding_version != SOMA_ENCODING_VERSION:
        raise ValueError(f"Unsupported SOMA object encoding version {encoding_version}")

    return obj_type


def _to_array_class(name: str) -> Type["tiledb_array.TileDBArray"]:
    options: Dict[str, Type[tiledb_array.TileDBArray]] = {
        cls.soma_type.lower(): cls
        for cls in (
            dataframe.DataFrame,
            dense_nd_array.DenseNDArray,
            sparse_nd_array.SparseNDArray,
        )
    }
    try:
        return options[name.lower()]
    except KeyError:
        raise SOMAError(f"{name!r} is not a valid SOMA data type") from None


@typeguard_ignore
def _to_group_class(name: str) -> Type["collection.AnyTileDBCollection"]:
    options: Dict[str, Type[collection.AnyTileDBCollection]] = {
        cls.soma_type.lower(): cls  # type: ignore[attr-defined, misc]
        for cls in (
            collection.Collection,
            experiment.Experiment,
            measurement.Measurement,
        )
    }
    try:
        return options[name.lower()]
    except KeyError:
        raise SOMAError(f"{name!r} is not a valid SOMA collection type") from None
