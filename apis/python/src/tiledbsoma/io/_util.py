# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

import pathlib
import sys
from contextlib import contextmanager
from typing import (
    ContextManager,
    Iterator,
    Union,
)
from unittest import mock

import anndata as ad
import pyarrow as pa
from anndata._core import file_backing

from .. import pytiledbsoma as clib
from .._exception import SOMAError
from .._types import Path
from ..options import SOMATileDBContext

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
    input_path: Path, *, mode: str = "r", ctx: SOMATileDBContext | None = None
) -> Iterator[ad.AnnData]:
    """
    This lets us ingest H5AD with "r" (backed mode) from S3 URIs.
    """
    ctx = ctx or SOMATileDBContext()
    vfs = clib.SOMAVFS(ctx.native_context)
    input_handle = clib.SOMAVFSFilebuf(vfs).open(str(input_path))
    try:
        with _hack_patch_anndata():
            anndata = ad.read_h5ad(_FSPathWrapper(input_handle, input_path), mode)
            yield anndata
    finally:
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

        def __new__(cls, _obj: object, path: Path) -> "_FSPathWrapper":
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


# @typeguard_ignore
def _hack_patch_anndata() -> ContextManager[object]:
    """Part Two of the ``_FSPathWrapper`` trick."""

    @file_backing.AnnDataFileManager.filename.setter  # type: ignore[misc]
    def filename(
        self: file_backing.AnnDataFileManager, filename: Union[Path, _FSPathWrapper]
    ) -> None:
        self._filename = filename

    return mock.patch.object(file_backing.AnnDataFileManager, "filename", filename)


def get_arrow_str_format(pa_type: pa.DataType) -> str:
    try:
        return _pa_type_to_str_fmt[pa_type]
    except KeyError:
        raise SOMAError(f"Could not convert {pa_type} to Arrow string format")
