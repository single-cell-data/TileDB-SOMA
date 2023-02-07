"""
This module exists to avoid what would otherwise be cyclic-module-import issues within
Collection.
"""

from typing import Callable, Dict, Optional, Type, TypeVar, Union, cast, overload

import somacore
from somacore import options

from . import (
    collection,
    dataframe,
    dense_nd_array,
    experiment,
    measurement,
    sparse_nd_array,
    tdb_handles,
    tiledb_object,
)
from .constants import (
    SOMA_ENCODING_VERSION,
    SOMA_ENCODING_VERSION_METADATA_KEY,
    SOMA_OBJECT_TYPE_METADATA_KEY,
)
from .exception import SOMAError
from .funcs import typeguard_ignore
from .options import SOMATileDBContext

_Obj = TypeVar("_Obj", bound="tiledb_object.AnyTileDBObject")
_Wrapper = TypeVar("_Wrapper", bound=tdb_handles.AnyWrapper)


@overload
def open(
    uri: str,
    mode: options.OpenMode = ...,
    *,
    soma_type: Optional[str] = None,
    context: Optional[SOMATileDBContext] = None,
) -> "tiledb_object.AnyTileDBObject":
    ...


@overload
def open(
    uri: str,
    mode: options.OpenMode,
    *,
    soma_type: Type[_Obj],
    context: Optional[SOMATileDBContext] = None,
) -> _Obj:
    ...


@typeguard_ignore
def open(
    uri: str,
    mode: options.OpenMode = "r",
    *,
    soma_type: Union[Type["tiledb_object.AnyTileDBObject"], str, None] = None,
    context: Optional[SOMATileDBContext] = None,
) -> "tiledb_object.AnyTileDBObject":
    """Opens a TileDB SOMA object.

    :param uri: The URI to open.
    :param mode: The mode to open in: ``r`` to read (default), ``w`` to write.
    :param soma_type: If set, the SOMA class you are expecting to get back.
        This can be provided as a SOMA type name
        If the stored SOMA object is not of the correct type, an error will be
        raised.
    :param context: If set, the ``SOMATileDBContext`` data to use.
    """
    context = context or SOMATileDBContext()
    obj = _open_internal(tdb_handles.open, uri, mode, context)
    try:
        if soma_type:
            if isinstance(soma_type, str):
                soma_type_name = soma_type
            elif issubclass(soma_type, somacore.SOMAObject):
                soma_type_name = soma_type.soma_type
            else:
                raise TypeError(f"cannot convert {soma_type!r} to expected SOMA type")
            if obj.soma_type.lower() != soma_type_name.lower():
                raise TypeError(
                    f"type of URI {uri!r} was {obj.soma_type}; expected {soma_type_name}"
                )
        return obj
    except Exception:
        obj.close()
        raise


def _open_internal(
    opener: Callable[[str, options.OpenMode, SOMATileDBContext], _Wrapper],
    uri: str,
    mode: options.OpenMode,
    context: SOMATileDBContext,
) -> "tiledb_object.TileDBObject[_Wrapper]":
    """Lower-level open function for internal use only."""
    handle = opener(uri, mode, context)
    try:
        return _reify_handle(handle)
    except Exception:
        handle.close()
        raise


@typeguard_ignore
def _reify_handle(hdl: _Wrapper) -> "tiledb_object.TileDBObject[_Wrapper]":
    """Picks out the appropriate SOMA class for a handle and wraps it."""
    typename = _read_soma_type(hdl)
    cls = _type_name_to_cls(typename)
    if cls._wrapper_type != type(hdl):
        raise SOMAError(
            f"cannot open {hdl.uri!r}: a {type(hdl._handle)}"
            f" cannot be converted to a {typename}"
        )
    return cast(
        tiledb_object.TileDBObject[_Wrapper],
        cls(hdl, _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code"),
    )


def _read_soma_type(hdl: tdb_handles.AnyWrapper) -> str:
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


@typeguard_ignore
def _type_name_to_cls(type_name: str) -> Type["tiledb_object.AnyTileDBObject"]:
    type_map: Dict[str, Type["tiledb_object.AnyTileDBObject"]] = {
        t.soma_type.lower(): t  # type: ignore[attr-defined, misc]  # spurious
        for t in (
            collection.Collection,
            dataframe.DataFrame,
            dense_nd_array.DenseNDArray,
            experiment.Experiment,
            measurement.Measurement,
            sparse_nd_array.SparseNDArray,
        )
    }
    try:
        return type_map[type_name.lower()]
    except KeyError as ke:
        options = sorted(type_map)
        raise SOMAError(
            f"{type_name!r} is not a recognized SOMA type; expected one of {options}"
        ) from ke
