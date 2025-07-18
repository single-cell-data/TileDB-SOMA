# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""This module exists to avoid what would otherwise be cyclic-module-import issues within
Collection.
"""

from __future__ import annotations

from typing import (
    TypeVar,
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
    _geometry_dataframe,
    _measurement,
    _multiscale_image,
    _point_cloud_dataframe,
    _scene,
    _soma_object,
    _sparse_nd_array,
    _tdb_handles,
)
from ._constants import (
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
    soma_type: str | None = None,
    context: SOMATileDBContext | None = None,
    tiledb_timestamp: OpenTimestamp | None = None,
) -> AnySOMAObject: ...


@overload
def open(
    uri: str,
    mode: options.OpenMode,
    *,
    soma_type: type[_Obj],
    context: SOMATileDBContext | None = None,
    tiledb_timestamp: OpenTimestamp | None = None,
) -> _Obj: ...


@typeguard_ignore
def open(
    uri: str,
    mode: options.OpenMode = "r",
    *,
    soma_type: type[SOMAObject] | str | None = None,  # type: ignore[type-arg]
    context: SOMATileDBContext | None = None,
    tiledb_timestamp: OpenTimestamp | None = None,
) -> AnySOMAObject:
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
    if soma_type is None:
        context = _validate_soma_tiledb_context(context)
        handle = _tdb_handles.open_handle_wrapper(uri=uri, mode=mode, context=context, timestamp=tiledb_timestamp)
        return reify_handle(handle)

    if isinstance(soma_type, str):
        soma_type_name = soma_type
    elif issubclass(soma_type, somacore.SOMAObject):
        soma_type_name = soma_type.soma_type
    else:
        raise TypeError(f"Cannot convert soma_type {soma_type!r} to expected SOMA type.")

    obj: AnySOMAObject = _type_name_to_cls(soma_type_name).open(
        uri=uri, mode=mode, context=context, tiledb_timestamp=tiledb_timestamp
    )
    if soma_type and obj.soma_type.lower() != soma_type_name.lower():
        obj.close()
        raise TypeError(f"Type of URI {uri!r} was {obj.soma_type}; expected {soma_type_name}.")
    return obj


@typeguard_ignore
def reify_handle(hdl: _Wrapper) -> SOMAObject[_Wrapper]:
    """Picks out the appropriate SOMA class for a handle and wraps it."""
    typename = hdl.metadata.get(SOMA_OBJECT_TYPE_METADATA_KEY)
    cls = _type_name_to_cls(typename)  # type: ignore[no-untyped-call]
    if not isinstance(hdl, cls._wrapper_type):
        raise SOMAError(f"cannot open {hdl.uri!r}: a {type(hdl._handle)} cannot be converted to a {typename}")
    return cast(
        "_soma_object.SOMAObject[_Wrapper]",
        cls(hdl, _dont_call_this_use_create_or_open_instead="tiledbsoma-internal-code"),
    )


@no_type_check
def _type_name_to_cls(type_name: str) -> type[AnySOMAObject]:
    type_map: dict[str, type[AnySOMAObject]] = {
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
