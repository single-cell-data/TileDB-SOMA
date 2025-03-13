# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

"""Abstractions to more easily manage read and write access to TileDB data.

``open``, ``ArrayWrapper.open``, ``GroupWrapper.open`` are the important parts.
"""

from __future__ import annotations

import abc
import enum
from typing import (
    Any,
    Dict,
    Generic,
    Iterator,
    List,
    Mapping,
    MutableMapping,
    Sequence,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
)

import attrs
import numpy as np
import pyarrow as pa
from somacore import options
from typing_extensions import Literal, Self

from . import pytiledbsoma as clib
from ._exception import DoesNotExistError, SOMAError, is_does_not_exist_error
from ._types import METADATA_TYPES, Metadatum, OpenTimestamp, StatusAndReason
from .options._soma_tiledb_context import SOMATileDBContext

AxisDomain = Union[None, Tuple[Any, Any], List[Any]]
Domain = Sequence[AxisDomain]

RawHandle = Union[
    clib.SOMAArray,
    clib.SOMADataFrame,
    clib.SOMAPointCloudDataFrame,
    clib.SOMAGeometryDataFrame,
    clib.SOMASparseNDArray,
    clib.SOMADenseNDArray,
    clib.SOMAGroup,
    clib.SOMACollection,
    clib.SOMAMeasurement,
    clib.SOMAExperiment,
    clib.SOMAScene,
    clib.SOMAMultiscaleImage,
]
_RawHdl_co = TypeVar("_RawHdl_co", bound=RawHandle, covariant=True)
"""A raw TileDB object. Covariant because Handles are immutable enough."""

_SOMAObjectType = TypeVar("_SOMAObjectType", bound=clib.SOMAObject)


def open(
    uri: str,
    mode: options.OpenMode,
    context: SOMATileDBContext,
    timestamp: OpenTimestamp | None,
    clib_type: str | None = None,
) -> "Wrapper[RawHandle]":
    """Determine whether the URI is an array or group, and open it."""
    open_mode = clib.OpenMode.read if mode == "r" else clib.OpenMode.write

    timestamp_ms = context._open_timestamp_ms(timestamp)

    _type_to_open = {
        "somaarray": clib.SOMAArray.open,
        "somagroup": clib.SOMAGroup.open,
        "somadataframe": clib.SOMADataFrame.open,
        "somapointclouddataframe": clib.SOMAPointCloudDataFrame.open,
        "somageometrydataframe": clib.SOMAGeometryDataFrame.open,
        "somadensendarray": clib.SOMADenseNDArray.open,
        "somasparsendarray": clib.SOMASparseNDArray.open,
        "somacollection": clib.SOMACollection.open,
        "somaexperiment": clib.SOMAExperiment.open,
        "somameasurement": clib.SOMAMeasurement.open,
        "somascene": clib.SOMAScene.open,
        "somamultiscaleimage": clib.SOMAMultiscaleImage.open,
        None: clib.SOMAObject.open,
    }

    if clib_type is not None:
        clib_type = clib_type.lower()

    try:
        soma_object = _type_to_open[clib_type](
            uri=uri,
            mode=open_mode,
            context=context.native_context,
            timestamp=(0, timestamp_ms),
        )
    except (RuntimeError, SOMAError) as tdbe:
        if is_does_not_exist_error(tdbe):
            raise DoesNotExistError(tdbe) from tdbe
        raise SOMAError(tdbe) from tdbe

    _type_to_class = {
        "somadataframe": DataFrameWrapper,
        "somapointclouddataframe": PointCloudDataFrameWrapper,
        "somageometrydataframe": GeometryDataFrameWrapper,
        "somadensendarray": DenseNDArrayWrapper,
        "somasparsendarray": SparseNDArrayWrapper,
        "somacollection": CollectionWrapper,
        "somaexperiment": ExperimentWrapper,
        "somameasurement": MeasurementWrapper,
        "somascene": SceneWrapper,
        "somamultiscaleimage": MultiscaleImageWrapper,
    }

    try:
        return _type_to_class[soma_object.type.lower()]._from_soma_object(
            soma_object, context
        )
    except KeyError:
        if soma_object.type.lower() == "somageometrydataframe":
            raise NotImplementedError(
                f"Support for {soma_object.type!r} is not yet implemented."
            )
        raise SOMAError(f"{uri!r} has unknown storage type {soma_object.type!r}")


@attrs.define(eq=False, hash=False, slots=False)
class Wrapper(Generic[_RawHdl_co], metaclass=abc.ABCMeta):
    """Wrapper for TileDB handles to manage lifecycle and metadata.

    Callers may read and use (non-underscored) members but should never set
    attributes on instances.
    """

    uri: str
    mode: options.OpenMode
    context: SOMATileDBContext
    timestamp_ms: int
    _handle: _RawHdl_co
    closed: bool = attrs.field(default=False, init=False)
    clib_type: str | None = None

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: OpenTimestamp | None,
    ) -> Self:
        if mode not in ("r", "w"):
            raise ValueError(f"Invalid open mode {mode!r}")
        timestamp_ms = context._open_timestamp_ms(timestamp)
        try:
            tdb = cls._opener(uri, mode, context, timestamp_ms)
            handle = cls(uri, mode, context, timestamp_ms, tdb)
            if mode == "w":
                with cls._opener(uri, "r", context, timestamp_ms) as auxiliary_reader:
                    handle._do_initial_reads(auxiliary_reader)
            else:
                handle._do_initial_reads(tdb)

        except RuntimeError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(tdbe) from tdbe
            raise
        return handle

    @classmethod
    def _from_soma_object(
        cls, soma_object: clib.SOMAObject, context: SOMATileDBContext
    ) -> Self:
        uri = soma_object.uri
        mode = soma_object.mode
        timestamp = context._open_timestamp_ms(soma_object.timestamp)
        try:
            handle = cls(uri, mode, context, timestamp, soma_object)
            if handle.mode == "w":
                with cls._opener(uri, mode, context, timestamp) as auxiliary_reader:
                    handle._do_initial_reads(auxiliary_reader)
            else:
                handle._do_initial_reads(soma_object)

        except RuntimeError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(tdbe) from tdbe
            raise
        return handle

    @classmethod
    @abc.abstractmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> _RawHdl_co:
        """Opens and returns a TileDB object specific to this type."""
        raise NotImplementedError()

    def reopen(
        self, mode: options.OpenMode, timestamp: OpenTimestamp | None
    ) -> clib.SOMAObject:
        if mode not in ("r", "w"):
            raise ValueError(
                f"Invalid mode '{mode}' passed. " "Valid modes are 'r' and 'w'."
            )
        ts = self.context._open_timestamp_ms(timestamp)
        self.metadata._write()
        return self._handle.reopen(
            clib.OpenMode.read if mode == "r" else clib.OpenMode.write, (0, ts)
        )

    # Covariant types should normally not be in parameters, but this is for
    # internal use only so it's OK.
    def _do_initial_reads(self, reader: _RawHdl_co) -> None:  # type: ignore[misc]
        """Final setup step before returning the Handle.

        This is passed a raw TileDB object opened in read mode, since writers
        will need to retrieve data from the backing store on setup.
        """
        # non–attrs-managed field
        self.metadata = MetadataWrapper(self, dict(reader.meta))

    @property
    def reader(self) -> _RawHdl_co:
        """Accessor to assert that you are working in read mode."""
        if self.closed:
            raise SOMAError(f"{self} is closed")
        if self.mode == "r":
            return self._handle
        raise SOMAError(f"cannot read from {self}; it is open for writing")

    @property
    def writer(self) -> _RawHdl_co:
        """Accessor to assert that you are working in write mode."""
        if self.closed:
            raise SOMAError(f"{self} is closed")
        if self.mode == "w":
            return self._handle
        raise SOMAError(f"cannot write to {self}; it is open for reading")

    def close(self) -> None:
        if self.closed:
            return
        self.metadata._write()
        self._handle.close()
        self.closed = True

    def _check_open(self) -> None:
        if self.closed:
            raise SOMAError(f"{self!r} is closed")

    def __repr__(self) -> str:
        closed_str = " (closed)" if self.closed else ""
        return f"<{type(self).__name__} {self.mode} on {self.uri!r}{closed_str}>"

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()


AnyWrapper = Wrapper[RawHandle]
"""Non-instantiable type representing any Handle."""


@attrs.define(frozen=True)
class GroupEntry:
    uri: str
    wrapper_type: Type[AnyWrapper]

    @classmethod
    def from_soma_group_entry(cls, obj: Tuple[str, str]) -> "GroupEntry":
        uri, type = obj[0], obj[1]
        if type == "SOMAArray":
            return GroupEntry(uri, SOMAArrayWrapper)
        if type == "SOMAGroup":
            return GroupEntry(uri, SOMAGroupWrapper)
        raise SOMAError(f"internal error: unknown object type {uri}")


class SOMAGroupWrapper(Wrapper[_SOMAObjectType]):
    """Base class for Pybind11 SOMAGroupWrapper handles."""

    _WRAPPED_TYPE: Type[_SOMAObjectType]

    clib_type = "SOMAGroup"

    @classmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> clib.SOMAGroup:
        open_mode = clib.OpenMode.read if mode == "r" else clib.OpenMode.write
        return cls._WRAPPED_TYPE.open(
            uri,
            mode=open_mode,
            context=context.native_context,
            timestamp=(0, timestamp),
        )

    def _do_initial_reads(self, group: clib.SOMAGroup) -> None:
        super()._do_initial_reads(group)

        self.initial_contents = {
            name: GroupEntry.from_soma_group_entry(entry)
            for name, entry in group.members().items()
        }

    @property
    def meta(self) -> "MetadataWrapper":
        return self.metadata

    def members(self) -> Dict[str, Tuple[str, str]]:
        return cast(Dict[str, Tuple[str, str]], self._handle.members())


class CollectionWrapper(SOMAGroupWrapper[clib.SOMACollection]):
    """Wrapper around a Pybind11 CollectionWrapper handle."""

    _WRAPPED_TYPE = clib.SOMACollection


class ExperimentWrapper(SOMAGroupWrapper[clib.SOMAExperiment]):
    """Wrapper around a Pybind11 ExperimentWrapper handle."""

    _WRAPPED_TYPE = clib.SOMAExperiment


class MeasurementWrapper(SOMAGroupWrapper[clib.SOMAMeasurement]):
    """Wrapper around a Pybind11 MeasurementWrapper handle."""

    _WRAPPED_TYPE = clib.SOMAMeasurement


class MultiscaleImageWrapper(SOMAGroupWrapper[clib.SOMAMultiscaleImage]):
    """Wrapper around a Pybind11 MultiscaleImage handle."""

    _WRAPPED_TYPE = clib.SOMAMultiscaleImage


class SceneWrapper(SOMAGroupWrapper[clib.SOMAScene]):
    """Wrapper around a Pybind11 SceneWrapper handle."""

    _WRAPPED_TYPE = clib.SOMAScene


class SOMAArrayWrapper(Wrapper[_SOMAObjectType]):
    """Base class for Pybind11 SOMAArrayWrapper handles."""

    _WRAPPED_TYPE: Type[_SOMAObjectType]

    clib_type = "SOMAArray"

    @classmethod
    def _opener(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: int,
    ) -> clib.SOMAArray:
        open_mode = clib.OpenMode.read if mode == "r" else clib.OpenMode.write

        return cls._WRAPPED_TYPE.open(
            uri,
            mode=open_mode,
            context=context.native_context,
            timestamp=(0, timestamp),
        )

    def _do_initial_reads(self, reader: RawHandle) -> None:
        """Final setup step before returning the Handle.

        This is passed a raw TileDB object opened in read mode, since writers
        will need to retrieve data from the backing store on setup.
        """
        # non–attrs-managed field
        self.metadata = MetadataWrapper(self, dict(reader.meta))

    @property
    def schema(self) -> pa.Schema:
        return self._handle.schema

    def schema_config_options(self) -> clib.PlatformSchemaConfig:
        """Returns a class containing the TileDB platform configuration options that
        can be read from an array schema.
        """
        return self._handle.schema_config_options()

    def config_options_from_schema(self) -> clib.PlatformConfig:
        return self._handle.config_options_from_schema()

    @property
    def meta(self) -> "MetadataWrapper":
        return self.metadata

    @property
    def ndim(self) -> int:
        return len(self._handle.dimension_names)

    @property
    def domain(self) -> Tuple[Tuple[object, object], ...]:
        from ._util import _cast_domainish

        return _cast_domainish(self._handle.domain())

    @property
    def maxdomain(self) -> Tuple[Tuple[object, object], ...]:
        from ._util import _cast_domainish

        return _cast_domainish(self._handle.maxdomain())

    def non_empty_domain(self) -> Tuple[Tuple[object, object], ...]:
        from ._util import _cast_domainish

        return _cast_domainish(self._handle.non_empty_domain())

    @property
    def attr_names(self) -> Tuple[str, ...]:
        return tuple(
            f.name for f in self.schema if f.name not in self._handle.dimension_names
        )

    @property
    def dim_names(self) -> Tuple[str, ...]:
        return tuple(self._handle.dimension_names)

    @property
    def shape(self) -> Tuple[int, ...]:
        """Not implemented for DataFrame."""
        return cast(Tuple[int, ...], tuple(self._handle.shape))

    @property
    def maxshape(self) -> Tuple[int, ...]:
        """Not implemented for DataFrame."""
        return cast(Tuple[int, ...], tuple(self._handle.maxshape))

    @property
    def maybe_soma_joinid_shape(self) -> int | None:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    @property
    def maybe_soma_joinid_maxshape(self) -> int | None:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    @property
    def tiledbsoma_has_upgraded_shape(self) -> bool:
        """Not implemented for DataFrame."""
        raise NotImplementedError

    @property
    def tiledbsoma_has_upgraded_domain(self) -> bool:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    def resize(self, newshape: Sequence[Union[int, None]]) -> None:
        """Not implemented for DataFrame."""
        raise NotImplementedError

    def tiledbsoma_can_resize(
        self, newshape: Sequence[Union[int, None]]
    ) -> StatusAndReason:
        """Not implemented for DataFrame."""
        raise NotImplementedError

    def tiledbsoma_upgrade_shape(self, newshape: Sequence[Union[int, None]]) -> None:
        """Not implemented for DataFrame."""
        raise NotImplementedError

    def tiledbsoma_can_upgrade_shape(
        self, newshape: Sequence[Union[int, None]]
    ) -> StatusAndReason:
        """Not implemented for DataFrame."""
        raise NotImplementedError

    def resize_soma_joinid_shape(
        self, newshape: int, function_name_for_messages: str
    ) -> None:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    def can_resize_soma_joinid_shape(
        self, newshape: int, function_name_for_messages: str
    ) -> StatusAndReason:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    def upgrade_soma_joinid_shape(
        self, newshape: int, function_name_for_messages: str
    ) -> None:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    def can_upgrade_soma_joinid_shape(
        self, newshape: int, function_name_for_messages: str
    ) -> StatusAndReason:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    def upgrade_domain(
        self, newdomain: Domain, function_name_for_messages: str
    ) -> None:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    def can_upgrade_domain(
        self, newdomain: Domain, function_name_for_messages: str
    ) -> StatusAndReason:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    def change_domain(self, newdomain: Domain, function_name_for_messages: str) -> None:
        """Only implemented for DataFrame."""
        raise NotImplementedError

    def can_change_domain(
        self, newdomain: Domain, function_name_for_messages: str
    ) -> StatusAndReason:
        """Only implemented for DataFrame."""
        raise NotImplementedError


class DataFrameWrapper(SOMAArrayWrapper[clib.SOMADataFrame]):
    """Wrapper around a Pybind11 SOMADataFrame handle."""

    _WRAPPED_TYPE = clib.SOMADataFrame

    @property
    def count(self) -> int:
        return int(self._handle.count)

    def write(self, values: pa.RecordBatch) -> None:
        self._handle.write(values)

    @property
    def maybe_soma_joinid_shape(self) -> int | None:
        """Wrapper-class internals"""
        return cast(Union[int, None], self._handle.maybe_soma_joinid_shape)

    @property
    def maybe_soma_joinid_maxshape(self) -> int | None:
        """Wrapper-class internals"""
        return cast(Union[int, None], self._handle.maybe_soma_joinid_maxshape)

    @property
    def tiledbsoma_has_upgraded_domain(self) -> bool:
        """Wrapper-class internals"""
        return cast(bool, self._handle.tiledbsoma_has_upgraded_domain)

    def resize_soma_joinid_shape(
        self, newshape: int, function_name_for_messages: str
    ) -> None:
        """Wrapper-class internals"""
        self._handle.resize_soma_joinid_shape(newshape, function_name_for_messages)

    def can_resize_soma_joinid_shape(
        self, newshape: int, function_name_for_messages: str
    ) -> StatusAndReason:
        """Wrapper-class internals"""
        return cast(
            StatusAndReason,
            self._handle.can_resize_soma_joinid_shape(
                newshape, function_name_for_messages
            ),
        )

    def upgrade_soma_joinid_shape(
        self, newshape: int, function_name_for_messages: str
    ) -> None:
        """Wrapper-class internals"""
        self._handle.upgrade_soma_joinid_shape(newshape, function_name_for_messages)

    def can_upgrade_soma_joinid_shape(
        self, newshape: int, function_name_for_messages: str
    ) -> StatusAndReason:
        """Wrapper-class internals"""
        return cast(
            StatusAndReason,
            self._handle.can_upgrade_soma_joinid_shape(
                newshape, function_name_for_messages
            ),
        )

    def upgrade_domain(
        self, newdomain: Domain, function_name_for_messages: str
    ) -> None:
        """Wrapper-class internals"""
        self._handle.upgrade_domain(newdomain, function_name_for_messages)

    def can_upgrade_domain(
        self, newdomain: Domain, function_name_for_messages: str
    ) -> StatusAndReason:
        """Wrapper-class internals"""
        return cast(
            StatusAndReason,
            self._handle.can_upgrade_domain(newdomain, function_name_for_messages),
        )

    def change_domain(self, newdomain: Domain, function_name_for_messages: str) -> None:
        """Wrapper-class internals"""
        self._handle.change_domain(newdomain, function_name_for_messages)

    def can_change_domain(
        self, newdomain: Domain, function_name_for_messages: str
    ) -> StatusAndReason:
        """Wrapper-class internals"""
        return cast(
            StatusAndReason,
            self._handle.can_change_domain(newdomain, function_name_for_messages),
        )


class PointCloudDataFrameWrapper(SOMAArrayWrapper[clib.SOMAPointCloudDataFrame]):
    """Wrapper around a Pybind11 SOMAPointCloudDataFrame handle."""

    _WRAPPED_TYPE = clib.SOMAPointCloudDataFrame

    @property
    def count(self) -> int:
        return int(self._handle.count)

    def write(self, values: pa.RecordBatch) -> None:
        self._handle.write(values)


class GeometryDataFrameWrapper(SOMAArrayWrapper[clib.SOMAGeometryDataFrame]):
    """Wrapper around a Pybind11 SOMAGeometryDataFrame handle."""

    _WRAPPED_TYPE = clib.SOMAGeometryDataFrame

    @property
    def count(self) -> int:
        return int(self._handle.count)

    def write(self, values: pa.RecordBatch) -> None:
        self._handle.write(values)


class DenseNDArrayWrapper(SOMAArrayWrapper[clib.SOMADenseNDArray]):
    """Wrapper around a Pybind11 DenseNDArrayWrapper handle."""

    _WRAPPED_TYPE = clib.SOMADenseNDArray

    @property
    def tiledbsoma_has_upgraded_shape(self) -> bool:
        """Wrapper-class internals"""
        return cast(bool, self._handle.tiledbsoma_has_upgraded_shape)

    def resize(self, newshape: Sequence[Union[int, None]]) -> None:
        """Wrapper-class internals"""
        self._handle.resize(newshape)

    def tiledbsoma_can_resize(
        self, newshape: Sequence[Union[int, None]]
    ) -> StatusAndReason:
        """Wrapper-class internals"""
        return cast(StatusAndReason, self._handle.tiledbsoma_can_resize(newshape))

    def tiledbsoma_upgrade_shape(self, newshape: Sequence[Union[int, None]]) -> None:
        """Wrapper-class internals"""
        self._handle.tiledbsoma_upgrade_shape(newshape)

    def tiledbsoma_can_upgrade_shape(
        self, newshape: Sequence[Union[int, None]]
    ) -> StatusAndReason:
        """Wrapper-class internals"""
        return cast(
            StatusAndReason, self._handle.tiledbsoma_can_upgrade_shape(newshape)
        )


class SparseNDArrayWrapper(SOMAArrayWrapper[clib.SOMASparseNDArray]):
    """Wrapper around a Pybind11 SparseNDArrayWrapper handle."""

    _WRAPPED_TYPE = clib.SOMASparseNDArray

    @property
    def nnz(self) -> int:
        return int(self._handle.nnz())

    @property
    def tiledbsoma_has_upgraded_shape(self) -> bool:
        """Wrapper-class internals"""
        return cast(bool, self._handle.tiledbsoma_has_upgraded_shape)

    def resize(self, newshape: Sequence[Union[int, None]]) -> None:
        """Wrapper-class internals"""
        self._handle.resize(newshape)

    def tiledbsoma_can_resize(
        self, newshape: Sequence[Union[int, None]]
    ) -> StatusAndReason:
        """Wrapper-class internals"""
        return cast(StatusAndReason, self._handle.can_resize(newshape))

    def tiledbsoma_upgrade_shape(self, newshape: Sequence[Union[int, None]]) -> None:
        """Wrapper-class internals"""
        self._handle.tiledbsoma_upgrade_shape(newshape)

    def tiledbsoma_can_upgrade_shape(
        self, newshape: Sequence[Union[int, None]]
    ) -> StatusAndReason:
        """Wrapper-class internals"""
        return cast(
            StatusAndReason, self._handle.tiledbsoma_can_upgrade_shape(newshape)
        )


class _DictMod(enum.Enum):
    """State machine to keep track of modifications to a dictionary.

    This whole thing is a hack to allow users to treat the metadata dict
    like an actual dictionary because tiledb currently does not support multiple
    modifications to the same key (e.g., add-then-delete a metadata entry has
    undesired results) [sc-25089].
    """

    # Initially-absent keys are either added or not (added then removed).
    ABSENT = enum.auto()
    """The key is not present in the dict. Initial state."""
    ADDED = enum.auto()
    """The key was originally ABSENT but has been added."""

    # Initially-present keys can be either updated or deleted.
    PRESENT = enum.auto()
    """The key is in the dict and is unchanged. Initial state."""
    UPDATED = enum.auto()
    """The key was originally PRESENT but has been changed."""
    DELETED = enum.auto()
    """The key was originally PRESENT but has been deleted."""

    @classmethod
    def start_state(cls, dct: Mapping[Any, Any], key: Any) -> "_DictMod":
        """Returns the starting state for a DictMod given the key of dct."""
        return cls.PRESENT if key in dct else cls.ABSENT

    def next_state(self, action: Literal["set", "del"]) -> "_DictMod":
        """Determines the next state of an entry given the action."""
        return {
            _DictMod.ABSENT: {
                "set": _DictMod.ADDED,
            },
            _DictMod.ADDED: {
                "set": _DictMod.ADDED,
                "del": _DictMod.ABSENT,
            },
            _DictMod.PRESENT: {
                "set": _DictMod.UPDATED,
                "del": _DictMod.DELETED,
            },
            _DictMod.UPDATED: {
                "set": _DictMod.UPDATED,
                "del": _DictMod.DELETED,
            },
            _DictMod.DELETED: {
                "set": _DictMod.UPDATED,
            },
        }[self][action]


@attrs.define(frozen=True)
class MetadataWrapper(MutableMapping[str, Any]):
    """A wrapper storing the metadata of some TileDB object.

    Because the view of metadata does not change after open time, we immediately
    cache all of it and use that to handle all reads. Writes are then proxied
    through to the backing store and the cache is updated to match.
    """

    owner: Wrapper[RawHandle]
    cache: Dict[str, Any]
    _mods: Dict[str, "_DictMod"] = attrs.field(init=False, factory=dict)
    """Tracks the modifications we have made to cache entries."""

    def __len__(self) -> int:
        self.owner._check_open()
        return len(self.cache)

    def __iter__(self) -> Iterator[str]:
        self.owner._check_open()
        return iter(self.cache)

    def __getitem__(self, key: str) -> Any:
        self.owner._check_open()
        return self.cache[key]

    def __setitem__(self, key: str, value: Any) -> None:
        self.owner.writer  # Ensures we're open in write mode.
        state = self._current_state(key)
        _check_metadata_type(key, value)
        self.cache[key] = value
        self._mods[key] = state.next_state("set")

    def __delitem__(self, key: str) -> None:
        self.owner.writer  # Ensures we're open in write mode.
        state = self._current_state(key)
        del self.cache[key]
        self._mods[key] = state.next_state("del")

    def _current_state(self, key: str) -> _DictMod:
        return self._mods.get(key, _DictMod.start_state(self.cache, key))

    def _write(self) -> None:
        """Writes out metadata changes, if there were any."""
        if not self._mods:
            # There were no changes (e.g., it's a read handle).  Do nothing.
            return
        # Only try to get the writer if there are changes to be made.

        errors: list[Exception] = []
        for key, mod in self._mods.items():
            try:
                if mod in (_DictMod.ADDED, _DictMod.UPDATED):
                    set_metadata = self.owner._handle.set_metadata
                    val = self.cache[key]
                    if isinstance(val, str):
                        set_metadata(key, np.array([val.encode("UTF-8")], "S"))
                    elif isinstance(val, bytes):
                        set_metadata(key, np.array([val], "V"))
                    else:
                        set_metadata(key, np.array([val]))
                if mod is _DictMod.DELETED:
                    self.owner._handle.delete_metadata(key)
            except Exception as e:
                errors.append(e)

        # Temporary hack: When we flush writes, note that the cache
        # is back in sync with disk.
        self._mods.clear()

        if errors:
            raise SOMAError(
                f"[MetadataWrapper][_write] {len(errors)} errors occured while writing metadata to disk."
            )

    def __repr__(self) -> str:
        prefix = f"{type(self).__name__}({self.owner})"
        if self.owner.closed:
            return f"<{prefix}>"
        return f"<{prefix} {self.cache}>"


def _check_metadata_type(key: str, obj: Metadatum) -> None:
    """Pre-checks that a metadata entry can be stored in an array.

    These checks are reproduced from the TileDB Python metadata-setting methods,
    but are slightly more restrictive than what TileDB allows in general:
    TileDB allows (some) arrays as metadata values, but the SOMA spec does not
    allow arrays of any kind.

    We have to pre-check since we don't write metadata changes until closing.
    """
    if not isinstance(key, str):
        raise TypeError(f"metadata keys must be strings, not {type(key)}")
    if isinstance(obj, METADATA_TYPES):
        return
    raise TypeError(f"cannot store {type(obj)} instance as metadata")
