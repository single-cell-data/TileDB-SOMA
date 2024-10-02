# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""This module exists to avoid what would otherwise be cyclic-module-import issues within
Collection.
"""

from typing import (
    Callable,
    Dict,
    Optional,
    Type,
    TypeVar,
    Union,
    cast,
    no_type_check,
    overload,
)

import somacore
from somacore import options

from . import (
    _collection,
    _dataframe,
    _dense_nd_array,
    _experiment,
    _measurement,
    _multiscale_image,
    _point_cloud_dataframe,
    _scene,
    _soma_object,
    _sparse_nd_array,
    _tdb_handles,
)
from ._constants import (
    SOMA_ENCODING_VERSION_METADATA_KEY,
    SOMA_OBJECT_TYPE_METADATA_KEY,
)
from ._exception import SOMAError
from ._funcs import typeguard_ignore
from ._soma_object import AnySOMAObject, SOMAObject
from ._types import OpenTimestamp
from .options import SOMATileDBContext
from .options._soma_tiledb_context import _validate_soma_tiledb_context

_Obj = TypeVar("_Obj", bound="_soma_object.AnySOMAObject")
_Wrapper = TypeVar("_Wrapper", bound=_tdb_handles.AnyWrapper)


@overload
def open(
    uri: str,
    mode: options.OpenMode = ...,
    *,
    soma_type: Optional[str] = None,
    context: Optional[SOMATileDBContext] = None,
    tiledb_timestamp: Optional[OpenTimestamp] = None,
) -> AnySOMAObject: ...


@overload
def open(
    uri: str,
    mode: options.OpenMode,
    *,
    soma_type: Type[_Obj],
    context: Optional[SOMATileDBContext] = None,
    tiledb_timestamp: Optional[OpenTimestamp] = None,
) -> _Obj: ...


@typeguard_ignore
def open(
    uri: str,
    mode: options.OpenMode = "r",
    *,
    soma_type: Union[Type[SOMAObject], str, None] = None,  # type: ignore[type-arg]
    context: Optional[SOMATileDBContext] = None,
    tiledb_timestamp: Optional[OpenTimestamp] = None,
) -> AnySOMAObject:
    """Opens a TileDB SOMA object.

    Args:
        uri:
            The URI to open.
        mode:
            The mode to open in: ``r`` to read (default), ``w`` to write.
        soma_type:
            If set, the SOMA class you are expecting to get back.
            This can be provided as a SOMA type name.
            If the stored SOMA object is not of the correct type, an error will be
            raised.
        context:
            If set, the :class:`SOMATileDBContext` data to use.
        tiledb_timestamp:
            If specified, overrides the default timestamp
            used to open this object. If unset, uses the timestamp provided by
            the context.

    Returns:
        The TileDB SOMA object.

    Raises:
        DoesNotExistError:
            If the object named by URI can not be accessed.
        SOMAError:
            If the underlying TileDB object specified by ``uri`` is
            not recognized as a SOMA object.
        TypeError:
            If the opened SOMA object type does not match the user-
            specified ``soma_type`` parameter.
        TypeError:
            If the user-provided ``soma_type`` parameter is not a
            recognizable type name or value.
        ValueError:
            If the user-provided ``mode`` is invalid.

    Lifecycle:
        Maturing.
    """
    context = _validate_soma_tiledb_context(context)
    obj: SOMAObject[_Wrapper] = _open_internal(  # type: ignore[valid-type]
        _tdb_handles.open, uri, mode, context, tiledb_timestamp
    )
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
    opener: Callable[
        [str, options.OpenMode, SOMATileDBContext, Optional[OpenTimestamp]], _Wrapper
    ],
    uri: str,
    mode: options.OpenMode,
    context: SOMATileDBContext,
    timestamp: Optional[OpenTimestamp],
) -> SOMAObject[_Wrapper]:
    """Lower-level open function for internal use only."""
    handle = opener(uri, mode, context, timestamp)
    try:
        return reify_handle(handle)
    except Exception:
        handle.close()
        raise


@typeguard_ignore
def reify_handle(hdl: _Wrapper) -> SOMAObject[_Wrapper]:
    """Picks out the appropriate SOMA class for a handle and wraps it."""
    typename = _read_soma_type(hdl)
    cls = _type_name_to_cls(typename)  # type: ignore[no-untyped-call]
    if not isinstance(hdl, cls._wrapper_type):
        raise SOMAError(
            f"cannot open {hdl.uri!r}: a {type(hdl._handle)}"
            f" cannot be converted to a {typename}"
        )
    return cast(
        _soma_object.SOMAObject[_Wrapper],
        cls(hdl, _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code"),
    )


def _read_soma_type(hdl: _tdb_handles.AnyWrapper) -> str:
    obj_type = hdl.metadata.get(SOMA_OBJECT_TYPE_METADATA_KEY)
    encoding_version = hdl.metadata.get(SOMA_ENCODING_VERSION_METADATA_KEY)

    if obj_type is None:
        raise SOMAError(
            f"stored TileDB object does not have {SOMA_OBJECT_TYPE_METADATA_KEY!r}"
        )

    if isinstance(obj_type, bytes):
        obj_type = str(obj_type, "utf-8")

    if not isinstance(obj_type, str):
        raise SOMAError(
            f"stored TileDB object {SOMA_OBJECT_TYPE_METADATA_KEY!r}"
            f" is {type(obj_type)}"
        )
    if encoding_version is None:
        raise SOMAError(
            f"stored TileDB object does not have {SOMA_ENCODING_VERSION_METADATA_KEY}"
        )

    if isinstance(encoding_version, bytes):
        encoding_version = str(encoding_version, "utf-8")

    if encoding_version not in {"1", "1.1.0"}:
        raise ValueError(f"Unsupported SOMA object encoding version {encoding_version}")

    return obj_type


@no_type_check
def _type_name_to_cls(type_name: str) -> Type[AnySOMAObject]:
    type_map: Dict[str, Type[AnySOMAObject]] = {
        t.soma_type.lower(): t
        for t in (
            _collection.Collection,
            _dataframe.DataFrame,
            _dense_nd_array.DenseNDArray,
            _experiment.Experiment,
            _measurement.Measurement,
            _multiscale_image.MultiscaleImage,
            _sparse_nd_array.SparseNDArray,
            _scene.Scene,
            _point_cloud_dataframe.PointCloudDataFrame,
        )
    }
    try:
        return type_map[type_name.lower()]
    except KeyError as ke:
        options = sorted(type_map)
        raise SOMAError(
            f"{type_name!r} is not a recognized SOMA type; expected one of {options}"
        ) from ke
