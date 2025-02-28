# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from pathlib import Path
from typing import Any, Optional, Tuple, Union

import h5py
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

    def __init__(self, input_path: Union[str, Path]):
        self._path = input_path
        self._version: Optional[Tuple[int, int, int]] = None
        self._dataset: Optional[h5py.File] = None
        self._nobs: Optional[int] = None
        self._nvars: Optional[int] = None

    def open(self) -> None:
        if self._dataset is None:
            self._dataset = h5py.File(self._path, "r")

    def __enter__(self) -> Self:
        self.open()
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def _read_version(self, dataset: h5py.File) -> Tuple[int, int, int]:
        try:
            version = dataset.attrs["software_version"]
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
        if not version.startswith(version_prefix):
            raise SOMAError(
                f"Unexpected value {version} for software version in gene expresion "
                f"file {self._path}."
            )
        version = version[len(version_prefix) :].split(".")
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
    def version(self) -> Tuple[int, int, int]:
        if self._version is not None:
            return self._version
        if self._dataset is None:
            with h5py.File(self._path) as dataset:
                return self._read_version(dataset)
        return self._read_version(self._dataset)
