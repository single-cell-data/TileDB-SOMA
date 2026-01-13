# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""This module exists to avoid what would otherwise be cyclic-module-import issues within
Collection.
"""

from __future__ import annotations

from typing import (
    TypeVar,
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
    _geometry_dataframe,
    _measurement,
    _multiscale_image,
    _point_cloud_dataframe,
    _scene,
    _sparse_nd_array,
    _tdb_handles,
)
from . import pytiledbsoma as clib
from ._exception import DoesNotExistError, SOMAError, is_does_not_exist_error
from ._funcs import typeguard_ignore
from ._soma_context import SOMAContext
from ._soma_object import SOMAObject
from ._types import OpenTimestamp
from ._util import tiledb_timestamp_to_ms
from .options import SOMATileDBContext, _update_context_and_timestamp

_Obj = TypeVar("_Obj", bound="SOMAObject")


@overload
def open(
    uri: str,
    mode: options.OpenMode = ...,
    *,
    soma_type: str | None = None,
    context: SOMAContext | SOMATileDBContext | None = None,
    tiledb_timestamp: OpenTimestamp | None = None,
) -> SOMAObject: ...


@overload
def open(
    uri: str,
    mode: options.OpenMode,
    *,
    soma_type: type[_Obj],
    context: SOMAContext | SOMATileDBContext | None = None,
    tiledb_timestamp: OpenTimestamp | None = None,
) -> _Obj: ...


@typeguard_ignore
def open(
    uri: str,
    mode: options.OpenMode = "r",
    *,
    soma_type: type[SOMAObject] | str | None = None,
    context: SOMAContext | SOMATileDBContext | None = None,
    tiledb_timestamp: OpenTimestamp | None = None,
) -> SOMAObject:
    """Opens a TileDB SOMA object.

    Args:
        uri:
            The URI to open.
        mode:
            The mode to open the object in.
            - ``r``: Open to read.
            - ``w``: Open to write.
            - ``d``: Open to delete.
        soma_type:
            If set, the SOMA class you are expecting to get back.
            This can be provided as a SOMA type name.
            If the stored SOMA object is not of the correct type, an error will be
            raised.
        context:
            If provided, the :class:`SOMAContext` to use when creating and opening this collection. If not,
            provide the default context will be used and possibly initialized.
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
    if soma_type is None:
        context, tiledb_timestamp = _update_context_and_timestamp(context, tiledb_timestamp)
        return _open_soma_object(uri, mode, context, tiledb_timestamp)

    if isinstance(soma_type, str):
        soma_type_name = soma_type
    elif issubclass(soma_type, somacore.SOMAObject):
        soma_type_name = soma_type.soma_type
    else:
        raise TypeError(f"Cannot convert soma_type {soma_type!r} to expected SOMA type.")

    obj: SOMAObject = _type_name_to_cls(soma_type_name).open(
        uri=uri, mode=mode, context=context, tiledb_timestamp=tiledb_timestamp
    )
    if soma_type and obj.soma_type.lower() != soma_type_name.lower():
        obj.close()
        raise TypeError(f"Type of URI {uri!r} was {obj.soma_type}; expected {soma_type_name}.")
    return obj


def _open_soma_object(
    uri: str,
    mode: options.OpenMode,
    context: SOMAContext,
    tiledb_timestamp: OpenTimestamp | None,
    clib_type: str | None = None,
) -> SOMAObject:
    """Picks out the appropriate SOMA class for a handle and wraps it."""
    if clib_type is None or clib_type.lower() in ["somaarray", "somagroup"]:
        timestamp_ms = tiledb_timestamp_to_ms(tiledb_timestamp)
        open_mode = _tdb_handles._open_mode_to_clib_mode(mode)
        try:
            handle = clib.SOMAObject.open(
                uri=uri,
                mode=open_mode,
                context=context.native_context,
                timestamp=(0, timestamp_ms),
                clib_type=clib_type,
            )
        except Exception as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(tdbe) from tdbe
            raise
        try:
            cls: type[SOMAObject] = _type_name_to_cls(handle.type.lower())
            return cls(
                handle, uri=uri, context=context, _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code"
            )
        except KeyError:
            raise SOMAError(f"{uri!r} has unknown storage type {clib_type!r}") from None

    try:
        cls = _type_name_to_cls(clib_type.lower())
    except KeyError:
        raise SOMAError(f"{uri!r} has unknown storage type {clib_type!r}") from None
    return cls.open(uri=uri, mode=mode, context=context, tiledb_timestamp=timestamp_ms)


@no_type_check
def _type_name_to_cls(type_name: str) -> type[SOMAObject]:
    type_map: dict[str, type[SOMAObject]] = {
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
            _geometry_dataframe.GeometryDataFrame,
        )
    }
    try:
        return type_map[type_name.lower()]
    except KeyError as ke:
        options = sorted(type_map)
        raise SOMAError(f"{type_name!r} is not a recognized SOMA type; expected one of {options}") from ke
