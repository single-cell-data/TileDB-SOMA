# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
#

from __future__ import annotations

from pathlib import Path
from typing import Any, Tuple, Union

import h5py
import numpy as np
import pyarrow as pa
from typing_extensions import Self

from ..._exception import SOMAError


def _str_to_int(value: str) -> int:
    if not value.isdigit():
        raise ValueError("{value} is not an integer.")
    return int(value)


def _version_less_than(version: str, upper_bound: Tuple[int, int, int]) -> bool:
    split_version = version.split(".")
    try:
        major = _str_to_int(split_version[0])
        minor = _str_to_int(split_version[1])
        patch = _str_to_int(split_version[2])
    except ValueError as err:
        raise ValueError(f"Unable to parse version {version}.") from err
    return (
        major < upper_bound[0]
        or (major == upper_bound[0] and minor < upper_bound[1])
        or (
            major == upper_bound[0]
            and minor == upper_bound[1]
            and patch < upper_bound[2]
        )
    )


class SpaceRangerMatrixReader:
    """Reader for the SpaceRanger HDF5 matrix group.

    The matrix is stored with features (var) on the rows and spot-barcodes (obs) on the columns.

    The following information should be provided in the provided HDF5 file:

    HDF5 Metadata (called attributes):

    Matrix: Expression matrix stored in CSC format
        /matrix/data
        /matrix/indices
        /matrix/indptr

    Observation Data:
        /matrix/barcodes: h5 string

    Variable Data:
        /matrix/features/_all_tag_keys
        /matrix/features/feature_type  #
        /matrix/features/id  # Unique idenitier for variable
        /matrix/features/name  # Human-readable name for variable
        # TODO: Finish this

    """

    def __init__(self, input_path: Union[str, Path]):
        # File management.
        self._path = input_path
        self._dataset: h5py.File | None = None

        # Metadatata.
        self._version: tuple[int, int, int] | None = None
        self._nobs: int | None = None
        self._nvar: int | None = None

        # X matrix.
        self._data: pa.Array | None = None
        self._indices: pa.Array | None = None  # Row indices for features.
        self._indptr: pa.Array | None = None  # Compressed column ptrs for barcodes.

        # Obs values.
        self._barcodes: pa.Array | None = None

        # Var values.
        self._var_names: pa.Array | None = None  # features/name
        self._feature_types: pa.Array | None = None  # features/feature_type
        self._gene_ids: pa.Array | None = None  # features/id
        self._genome: pa.Array | None = None  # features/genome

    def open(self) -> None:
        if self._dataset is None:
            self._dataset = h5py.File(self._path, "r")
        self._nvar, self._nobs = self._dataset["matrix"]["shape"]

    def __enter__(self) -> Self:
        self.open()
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def _read_dataset(
        self,
        group: h5py.Group,
        name: str,
        desired_type: np.typing.DTypeLike | None = None,
    ) -> np.typing.NDArray:
        if desired_type is None:
            return group[name]
        return group[name].astype(desired_type)

    def _read_version(self) -> tuple[int, int, int]:
        if self._dataset is None:
            raise RuntimeError("Open dataset to read version.")
        try:
            version = self._dataset.attrs["software_version"]
        except KeyError as ke:
            raise SOMAError(
                f"Unable to read software version from gene expression file "
                f"{self._path}."
            ) from ke
        if not isinstance(version, str):
            raise SOMAError(
                f"Unexpected type {type(version)!r} for software version in gene "
                f"expression file {self._path}. Expected a string."
            )
        version_prefix = "spaceranger-"
        if version.startswith(version_prefix):
            version = version[len(version_prefix) :].split(".")
        else:
            version = version.split(".")
        if len(version) not in {3, 4}:
            raise SOMAError(
                f"Unexpected value {version} for software version in gene expresion "
                f"file {self._path}."
            )
        try:
            major = _str_to_int(version[0])
            minor = _str_to_int(version[1])
            patch = _str_to_int(version[2])
        except ValueError:
            raise SOMAError(
                f"Unexpected value {version} for software version in gene expresion "
                f"file {self._version}."
            )
        self._version = (major, minor, patch)
        return self._version

    def close(self) -> None:
        if self._dataset is not None:
            self._dataset.close()

    @property
    def nobs(self) -> int:
        if self._nobs is None:
            raise RuntimeError("Open reader to read nobs.")
        return self._nobs

    @property
    def nvar(self) -> int:
        if self._nvar is None:
            raise RuntimeError("Open reader to read nvar.")
        return self._nvar

    @property
    def version(self) -> tuple[int, int, int]:
        if self._version is None:
            return self._read_version()
        return self._version
