# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

"""Abstractions to more easily manage read and write access to TileDB data.

``open``, ``ArrayWrapper.open``, ``GroupWrapper.open`` are the important parts.
"""

import abc
from typing import (
    Any,
    Dict,
    Generic,
    Optional,
    Tuple,
    Type,
    TypeVar,
    Union,
    cast,
)

import attrs
from somacore import options
from typing_extensions import Self

from . import pytiledbsoma as clib
from ._exception import DoesNotExistError, SOMAError, is_does_not_exist_error
from ._metadata_wrapper import MetadataWrapper
from ._types import OpenTimestamp
from .options._soma_tiledb_context import SOMATileDBContext

CLibHandle = Union[
    clib.SOMAArray,
    clib.SOMADataFrame,
    clib.SOMAPointCloudDataFrame,
    clib.SOMASparseNDArray,
    clib.SOMADenseNDArray,
    clib.SOMAGroup,
    clib.SOMACollection,
    clib.SOMAMeasurement,
    clib.SOMAExperiment,
    clib.SOMAScene,
    clib.SOMAMultiscaleImage,
]
_CLibHandle_co = TypeVar("_CLibHandle_co", bound=CLibHandle, covariant=True)
"""A handle to a pybind11-managed libtiledbsoma object. Covariant because Handles are immutable enough."""


def open(
    uri: str,
    mode: options.OpenMode,
    context: SOMATileDBContext,
    timestamp: Optional[OpenTimestamp],
    clib_type: Optional[str] = None,
) -> "Wrapper[CLibHandle]":
    """Determine whether the URI is an array or group, and open it."""
    open_mode = clib.OpenMode.read if mode == "r" else clib.OpenMode.write

    timestamp_ms = context._open_timestamp_ms(timestamp)

    soma_object = clib.SOMAObject.open(
        uri=uri,
        mode=open_mode,
        context=context.native_context,
        timestamp=(0, timestamp_ms),
        clib_type=clib_type,
    )

    if not soma_object:
        raise DoesNotExistError(f"{uri!r} does not exist")

    _type_to_class = {
        "somadataframe": DataFrameWrapper,
        "somapointclouddataframe": PointCloudDataFrameWrapper,
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
class Wrapper(Generic[_CLibHandle_co], metaclass=abc.ABCMeta):
    """Wrapper for TileDB handles to manage lifecycle and metadata.

    Callers may read and use (non-underscored) members but should never set
    attributes on instances.
    """

    uri: str
    mode: options.OpenMode
    context: SOMATileDBContext
    timestamp_ms: int
    _handle: _CLibHandle_co
    closed: bool = attrs.field(default=False, init=False)
    clib_type: Optional[str] = None

    @classmethod
    def open(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
        timestamp: Optional[OpenTimestamp],
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
    ) -> _CLibHandle_co:
        """Opens and returns a TileDB object specific to this type."""
        raise NotImplementedError()

    def reopen(
        self, mode: options.OpenMode, timestamp: Optional[OpenTimestamp]
    ) -> clib.SOMAObject:
        if mode not in ("r", "w"):
            raise ValueError(
                f"Invalid mode '{mode}' passed. " "Valid modes are 'r' and 'w'."
            )
        ts = self.context._open_timestamp_ms(timestamp)
        return self._handle.reopen(
            clib.OpenMode.read if mode == "r" else clib.OpenMode.write, (0, ts)
        )

    # Covariant types should normally not be in parameters, but this is for
    # internal use only so it's OK.
    def _do_initial_reads(self, reader: _CLibHandle_co) -> None:  # type: ignore[misc]
        """Final setup step before returning the Handle.

        This is passed a raw TileDB object opened in read mode, since writers
        will need to retrieve data from the backing store on setup.
        """
        # non–attrs-managed field
        self.metadata = MetadataWrapper(self, dict(reader.meta))

    @property
    def reader(self) -> _CLibHandle_co:
        """Accessor to assert that you are working in read mode."""
        if self.closed:
            raise SOMAError(f"{self} is closed")
        if self.mode == "r":
            return self._handle
        raise SOMAError(f"cannot read from {self}; it is open for writing")

    @property
    def writer(self) -> _CLibHandle_co:
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


AnyWrapper = Wrapper[CLibHandle]
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


_ClibGroupType = TypeVar("_ClibGroupType", bound=clib.SOMAGroup)


class SOMAGroupWrapper(Wrapper[_ClibGroupType]):
    """Base class for Pybind11 SOMAGroupWrapper handles."""

    _GROUP_WRAPPED_TYPE: Type[_ClibGroupType]

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
        return cls._GROUP_WRAPPED_TYPE.open(
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

    _GROUP_WRAPPED_TYPE = clib.SOMACollection


class ExperimentWrapper(SOMAGroupWrapper[clib.SOMAExperiment]):
    """Wrapper around a Pybind11 ExperimentWrapper handle."""

    _GROUP_WRAPPED_TYPE = clib.SOMAExperiment


class MeasurementWrapper(SOMAGroupWrapper[clib.SOMAMeasurement]):
    """Wrapper around a Pybind11 MeasurementWrapper handle."""

    _GROUP_WRAPPED_TYPE = clib.SOMAMeasurement


class MultiscaleImageWrapper(SOMAGroupWrapper[clib.SOMAMultiscaleImage]):
    """Wrapper around a Pybind11 MultiscaleImage handle."""

    _GROUP_WRAPPED_TYPE = clib.SOMAMultiscaleImage


class SceneWrapper(SOMAGroupWrapper[clib.SOMAScene]):
    """Wrapper around a Pybind11 SceneWrapper handle."""

    _GROUP_WRAPPED_TYPE = clib.SOMAScene


_CLibArrayType = TypeVar("_CLibArrayType", bound=clib.SOMAArray)


class SOMAArrayWrapper(Wrapper[_CLibArrayType]):
    """Base class for Pybind11 SOMAArrayWrapper handles."""

    _ARRAY_WRAPPED_TYPE: Type[_CLibArrayType]

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

        return cls._ARRAY_WRAPPED_TYPE.open(
            uri,
            mode=open_mode,
            context=context.native_context,
            column_names=[],
            result_order=clib.ResultOrder.automatic,
            timestamp=(0, timestamp),
        )

    def _do_initial_reads(self, reader: CLibHandle) -> None:
        """Final setup step before returning the Handle.

        This is passed a raw TileDB object opened in read mode, since writers
        will need to retrieve data from the backing store on setup.
        """
        # non–attrs-managed field
        self.metadata = MetadataWrapper(self, dict(reader.meta))

    @property
    def meta(self) -> "MetadataWrapper":
        return self.metadata


class DataFrameWrapper(SOMAArrayWrapper[clib.SOMADataFrame]):
    """Wrapper around a Pybind11 SOMADataFrame handle."""

    _ARRAY_WRAPPED_TYPE = clib.SOMADataFrame


class PointCloudDataFrameWrapper(SOMAArrayWrapper[clib.SOMAPointCloudDataFrame]):
    """Wrapper around a Pybind11 SOMAPointCloudDataFrame handle."""

    _ARRAY_WRAPPED_TYPE = clib.SOMAPointCloudDataFrame


class DenseNDArrayWrapper(SOMAArrayWrapper[clib.SOMADenseNDArray]):
    """Wrapper around a Pybind11 DenseNDArrayWrapper handle."""

    _ARRAY_WRAPPED_TYPE = clib.SOMADenseNDArray


class SparseNDArrayWrapper(SOMAArrayWrapper[clib.SOMASparseNDArray]):
    """Wrapper around a Pybind11 SparseNDArrayWrapper handle."""

    _ARRAY_WRAPPED_TYPE = clib.SOMASparseNDArray
