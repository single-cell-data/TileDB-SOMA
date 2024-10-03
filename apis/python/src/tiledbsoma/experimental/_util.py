# Copyright (c) 2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2024 TileDB, Inc
#
# Licensed under the MIT License.
from pathlib import Path
from typing import Tuple, Union

import h5py

from .._exception import SOMAError


def _str_to_int(value: str) -> int:
    if not value.isdigit():
        raise ValueError("{value} is not an integer.")
    return int(value)


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
