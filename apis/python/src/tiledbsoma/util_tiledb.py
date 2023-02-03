from typing import Any, Callable, Generic, MutableMapping, Optional, TypeVar, cast

import attrs
import numpy as np
import pandas._typing as pdt
import scipy.sparse as sp
import tiledb
from pandas.api.types import infer_dtype, is_categorical_dtype
from somacore import options

# This package's pybind11 code
import tiledbsoma.libtiledbsoma as clib

from .exception import DoesNotExistError, SOMAError
from .options import SOMATileDBContext
from .types import NPNDArray, PDSeries, TDBHandle

SOMA_OBJECT_TYPE_METADATA_KEY = "soma_object_type"
SOMA_ENCODING_VERSION_METADATA_KEY = "soma_encoding_version"
SOMA_ENCODING_VERSION = "1"


_DT = TypeVar("_DT", bound=pdt.Dtype)
_MT = TypeVar("_MT", NPNDArray, sp.spmatrix, PDSeries)
_str_to_type = {"boolean": bool, "string": str, "bytes": bytes}


def _to_tiledb_supported_dtype(dtype: _DT) -> _DT:
    """A handful of types are cast into the TileDB type system."""
    # TileDB has no float16 -- cast up to float32
    return cast(_DT, np.dtype("float32")) if dtype == np.dtype("float16") else dtype


def to_tiledb_supported_array_type(x: _MT) -> _MT:
    """
    Converts datatypes unrepresentable by TileDB into datatypes it can represent.
    E.g., categorical strings -> string.

    See also [https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html](https://docs.scipy.org/doc/numpy-1.10.1/reference/arrays.dtypes.html).

    Preferentially converts to the underlying primitive type, as TileDB does not support
    most complex types. NOTE: this does not support ``datetime64`` conversion.

    Categoricals are a special case. If the underlying categorical type is a primitive,
    convert to that. If the array contains NA/NaN (i.e. not in the category, code == -1),
    raise error unless it is a float or string.
    """
    if isinstance(x, (np.ndarray, sp.spmatrix)) or not is_categorical_dtype(x):
        target_dtype = _to_tiledb_supported_dtype(x.dtype)
        return x if target_dtype == x.dtype else x.astype(target_dtype)

    categories = x.cat.categories
    cat_dtype = categories.dtype
    if cat_dtype.kind in ("f", "u", "i"):
        if x.hasnans and cat_dtype.kind == "i":
            raise ValueError(
                "Categorical array contains NaN -- unable to convert to TileDB array."
            )
        target_dtype = _to_tiledb_supported_dtype(cat_dtype)
    else:
        # Into the weirdness. See if Pandas can help with edge cases.
        inferred = infer_dtype(categories)
        if x.hasnans and inferred in ("boolean", "bytes"):
            raise ValueError(
                "Categorical array contains NaN -- unable to convert to TileDB array."
            )
        target_dtype = np.dtype(_str_to_type.get(inferred, object))

    return x.astype(target_dtype)


def is_does_not_exist_error(e: tiledb.TileDBError) -> bool:
    """ "
    Given a TileDBError, return true if it indicates the object does not exist.

    Example
    -------

    try:
        with tiledb.open(uri):
            ...
    except tiledb.TileDBError as e:
        if is_does_not_exist_error(e):
            ...
        raise e
    """
    stre = str(e)
    # Local-disk/S3 does-not-exist exceptions say 'Group does not exist'; TileDB Cloud
    # does-not-exist exceptions are worded less clearly.
    if "does not exist" in stre or "HTTP code 401" in stre:
        return True

    return False


def is_duplicate_group_key_error(e: tiledb.TileDBError) -> bool:
    """
    Given a TileDBError, return try if it indicates a duplicate member
    add request in a tiledb.Group.
    """
    stre = str(e)
    if "member already exists in group" in stre:
        return True

    return False


_TDBO = TypeVar("_TDBO", bound=TDBHandle)
_TDBO_co = TypeVar("_TDBO_co", bound=TDBHandle, covariant=True)
_Self = TypeVar("_Self")


# TODO: We should only need to keep around one of the read or write handles
# after we do the initial open (where we read metadata and group contents).
# For now, we're leaving this as-is.
@attrs.define()
class ReadWriteHandle(Generic[_TDBO_co]):
    """Paired read/write handles to a TileDB object.

    You can (and should) access the members, but do not reassign them.
    """

    uri: str
    context: SOMATileDBContext
    reader: _TDBO_co
    maybe_writer: Optional[_TDBO_co]
    closed: bool = False

    @classmethod
    def open_array(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
    ) -> "ReadWriteHandle[tiledb.Array]":
        return cls._open(tiledb.open, uri, mode, context)

    @classmethod
    def open_group(
        cls,
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
    ) -> "ReadWriteHandle[tiledb.Group]":
        return cls._open(tiledb.Group, uri, mode, context)

    @classmethod
    def _open(
        cls,
        tiledb_cls: Callable[..., _TDBO],
        uri: str,
        mode: options.OpenMode,
        context: SOMATileDBContext,
    ) -> "ReadWriteHandle[_TDBO]":
        try:
            rd = tiledb_cls(uri, "r", ctx=context.tiledb_ctx)
            try:
                wr = (
                    tiledb_cls(uri, "w", ctx=context.tiledb_ctx)
                    if mode == "w"
                    else None
                )
            except Exception:
                rd.close()
                raise
        except tiledb.TileDBError as tdbe:
            if is_does_not_exist_error(tdbe):
                raise DoesNotExistError(
                    f"{tiledb_cls.__name__} {uri!r} does not exist"
                ) from tdbe
            raise
        return ReadWriteHandle(uri, context, reader=rd, maybe_writer=wr)

    @property
    def writer(self) -> _TDBO_co:
        if self.maybe_writer is None:
            raise SOMAError(f"cannot write to {self.uri!r} opened in read-only mode")
        return self.maybe_writer

    @property
    def mode(self) -> options.OpenMode:
        return "r" if self.maybe_writer is None else "w"

    @property
    def metadata(self) -> MutableMapping[str, Any]:
        if self.maybe_writer is not None:
            return self.maybe_writer.meta  # type: ignore[no-any-return]
        return self.reader.meta  # type: ignore[no-any-return]

    def close(self) -> None:
        if self.closed:
            return
        self.closed = True
        self.reader.close()
        if self.maybe_writer is not None:
            self.maybe_writer.close()

    def _flush_hack(self) -> None:
        """Flushes writes. Use rarely, if ever."""
        if self.maybe_writer is None:
            return
        self.maybe_writer.close()
        if isinstance(self.maybe_writer, tiledb.Array):
            self.maybe_writer = tiledb.open(self.uri, "w", ctx=self.context.tiledb_ctx)
            self.reader.reopen()
        else:
            self.maybe_writer = tiledb.Group(self.uri, "w", ctx=self.context.tiledb_ctx)
            self.reader = tiledb.Group(self.uri, "r", ctx=self.context.tiledb_ctx)

    def __enter__(self: _Self) -> _Self:
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()


def stats_enable() -> None:
    """
    Enable TileDB internal statistics.
    """
    clib.stats_enable()


def stats_disable() -> None:
    """
    Disable TileDB internal statistics.
    """
    clib.stats_disable()


def stats_reset() -> None:
    """
    Reset all TileDB internal statistics to 0.
    """
    clib.stats_reset()


def stats_dump() -> None:
    """
    Print TileDB internal statistics.
    """
    clib.stats_dump()
