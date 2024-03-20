# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

from typing import (
    ContextManager,
    Optional,
    Tuple,
)
from unittest import mock

import anndata as ad
import tiledb
from anndata._core import file_backing

from .._types import Path


def read_h5ad(
    input_path: Path, *, mode: str = "r", ctx: Optional[tiledb.Ctx] = None
) -> Tuple[tiledb.vfs.FileIO, ad.AnnData]:
    """
    This lets us ingest H5AD with "r" (backed mode) from S3 URIs.  The caller must close the
    returned handle after processing the returned backed AnnData object.
    """

    # Ideally we'd do a with-open-as. However, since we're returning an AnnData
    # object in backed mode, everything the caller does will have to be done
    # with a still-open handle. Therefore, we need to return the backed
    # AnnData object _and_ the open handle.
    input_handle = tiledb.VFS(ctx=ctx).open(input_path)
    with _hack_patch_anndata():
        anndata = ad.read_h5ad(_FSPathWrapper(input_handle, input_path), mode)
    return (input_handle, anndata)


# This trick lets us ingest H5AD with "r" (backed mode) from S3 URIs.  While h5ad
# supports any file-like object, AnnData specifically wants only an `os.PathLike`
# object. The only thing it does with the PathLike is to use it to get the filename.
class _FSPathWrapper:
    """Tricks anndata into thinking a file-like object is an ``os.PathLike``.

    While h5ad supports any file-like object, anndata specifically wants
    an ``os.PathLike object``, which it uses *exclusively* to get the "filename"
    of the opened file.

    We need to provide ``__fspath__`` as a real class method, so simply
    setting ``some_file_obj.__fspath__ = lambda: "some/path"`` won't work,
    so here we just proxy all attributes except ``__fspath__``.
    """

    def __init__(self, obj: object, path: Path) -> None:
        self._obj = obj
        self._path = path

    def __fspath__(self) -> Path:
        return self._path

    def __getattr__(self, name: str) -> object:
        return getattr(self._obj, name)


def _hack_patch_anndata() -> ContextManager[object]:
    """Part Two of the ``_FSPathWrapper`` trick."""

    @file_backing.AnnDataFileManager.filename.setter  # type: ignore
    def filename(self: file_backing.AnnDataFileManager, filename: Path) -> None:
        self._filename = filename

    return mock.patch.object(file_backing.AnnDataFileManager, "filename", filename)
