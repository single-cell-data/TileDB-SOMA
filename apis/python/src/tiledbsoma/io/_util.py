# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.
import pathlib
from contextlib import contextmanager
from typing import (
    ContextManager,
    Iterator,
    Optional,
    Union,
)
from unittest import mock

import anndata as ad
from anndata._core import file_backing

import tiledb

from .._types import Path


@contextmanager
def read_h5ad(
    input_path: Path, *, mode: str = "r", ctx: Optional[tiledb.Ctx] = None
) -> Iterator[ad.AnnData]:
    """
    This lets us ingest H5AD with "r" (backed mode) from S3 URIs.
    """
    input_handle = tiledb.VFS(ctx=ctx).open(input_path)
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
