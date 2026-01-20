# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import pathlib
import sys
from collections.abc import Iterator
from contextlib import contextmanager

import anndata as ad
import pyarrow as pa
from anndata._core import file_backing

from tiledbsoma import pytiledbsoma as clib
from tiledbsoma._exception import SOMAError
from tiledbsoma._soma_context import SOMAContext
from tiledbsoma._types import Path
from tiledbsoma.options import SOMATileDBContext

from ._caching_reader import CachingReader

_pa_type_to_str_fmt = {
    pa.string(): "U",
    pa.binary(): "Z",
    pa.large_string(): "U",
    pa.large_binary(): "Z",
    pa.int8(): "c",
    pa.uint8(): "C",
    pa.int16(): "s",
    pa.uint16(): "S",
    pa.int32(): "i",
    pa.uint32(): "I",
    pa.int64(): "l",
    pa.uint64(): "L",
    pa.float32(): "f",
    pa.float64(): "g",
    pa.bool_(): "b",
}


@contextmanager
def read_h5ad(
    input_path: Path | str, *, mode: str | None = "r", ctx: SOMAContext | SOMATileDBContext | None = None
) -> Iterator[ad.AnnData]:
    """This lets us ingest H5AD with "r" (backed mode) from S3 URIs."""
    if ctx is None:
        if not SOMAContext.has_default():
            SOMAContext.set_default()
        ctx = SOMAContext.get_default()
    elif isinstance(ctx, SOMATileDBContext):
        ctx = ctx._to_soma_context()
    input_handle = CachingReader(
        clib.SOMAFileHandle(str(input_path), ctx._handle),
        memory_budget=64 * 1024**2,
        cache_block_size=8 * 1024**2,
    )
    try:
        adata = ad.read_h5ad(_FSPathWrapper(input_handle, input_path), mode)
        yield adata
    finally:
        # This prevents a race condition with the REPL cleanup in ipython. See sc-65863
        if "adata" in locals() and adata.file:
            adata.file.close()
        input_handle.close()


# This trick lets us ingest H5AD with "r" (backed mode) from S3 URIs.  While h5ad
# supports any file-like object, AnnData specifically wants only an `os.PathLike`
# object. The only thing it does with the PathLike is to use it to get the filename.
class _FSPathWrapper(pathlib.Path):
    """Tricks anndata into thinking a file-like object is an ``os.PathLike``.

    While h5ad supports any file-like object, anndata specifically wants
    an ``os.PathLike object``, which it uses *exclusively* to get the "filename"
    of the opened file.

    We need to provide ``__fspath__`` as a real class method, so simply
    setting ``some_file_obj.__fspath__ = lambda: "some/path"`` won't work,
    so here we just proxy all attributes except ``__fspath__``.
    """

    if sys.version_info >= (3, 12):

        def __init__(self, obj: object, path: Path) -> None:
            super().__init__(path)
            self._obj = obj
            self._path = path

    else:

        def __new__(cls, _obj: object, path: Path) -> _FSPathWrapper:
            return super().__new__(cls, path)

        # ``pathlib.Path`` construction references this attribute (``PosixFlavour`` or ``WindowsFlavour``)
        _flavour = pathlib.Path().__class__._flavour  # type: ignore[attr-defined]

        def __init__(self, obj: object, path: Path) -> None:
            super().__init__()
            self._obj = obj
            self._path = path

    def __fspath__(self) -> str:
        return self._path if isinstance(self._path, str) else str(self._path)

    def __getattr__(self, name: str) -> object:
        return getattr(self._obj, name)


def _monkey_patch_anndata() -> None:
    """Monkey patch the AnnData backed file manager to allow our path-like
    wrapper class in addition to a str|Path.

    As this is a global change whenever tiledbsoma is imported, take
    care to preserve original AnnData setter behavior in cases unrelated
    to the above use.
    """
    original_setter = file_backing.AnnDataFileManager.filename.fset
    original_getter = file_backing.AnnDataFileManager.filename.fget

    def filename(
        self: file_backing.AnnDataFileManager,
        filename: Path | _FSPathWrapper | None,
    ) -> None:
        if isinstance(filename, _FSPathWrapper):
            self._filename = filename
        else:
            original_setter(self, filename)

    file_backing.AnnDataFileManager.filename = property(original_getter, filename)


_monkey_patch_anndata()


def get_arrow_str_format(pa_type: pa.DataType) -> str:
    try:
        return _pa_type_to_str_fmt[pa_type]
    except KeyError:
        raise SOMAError(f"Could not convert {pa_type} to Arrow string format") from None


def _set_and_get_context(context: SOMAContext | SOMATileDBContext | None) -> SOMAContext:
    """Get a SOMAContext from input context parameter.

    If no context is provided, get and possibly set the default context.
    """
    if context is None:
        if not SOMAContext.has_default():
            SOMAContext.set_default()
        return SOMAContext.get_default()
    if isinstance(context, SOMATileDBContext):
        return context._to_soma_context()
    return context
