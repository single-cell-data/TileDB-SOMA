# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from pathlib import Path
from typing import Tuple, Union

import h5py

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


def _read_visium_software_version(
    gene_expression_path: Union[str, Path]
) -> Tuple[int, int, int]:
    with h5py.File(gene_expression_path) as dataset:
        try:
            version = dataset.attrs["software_version"]
        except KeyError as ke:
            raise SOMAError(
                f"Unable to read software version from gene expression file "
                f"{gene_expression_path}."
            ) from ke
    if not isinstance(version, str):
        raise SOMAError(
            f"Unexpected type {type(version)!r} for software version in gene "
            f"expression file {gene_expression_path}. Expected a string."
        )
    version_prefix = "spaceranger-"
    if not version.startswith(version_prefix):
        raise SOMAError(
            f"Unexpected value {version} for software version in gene expresion "
            f"file {gene_expression_path}."
        )
    version = version[len(version_prefix) :].split(".")
    if len(version) not in {3, 4}:
        raise SOMAError(
            f"Unexpected value {version} for software version in gene expresion "
            f"file {gene_expression_path}."
        )
    try:
        major = _str_to_int(version[0])
        minor = _str_to_int(version[1])
        patch = _str_to_int(version[2])
    except ValueError:
        raise SOMAError(
            f"Unexpected value {version} for software version in gene expresion "
            f"file {gene_expression_path}."
        )
    return (major, minor, patch)
