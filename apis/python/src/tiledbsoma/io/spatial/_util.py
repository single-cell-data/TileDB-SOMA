# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.

from __future__ import annotations

from pathlib import Path
from typing import Any

import h5py
import numpy as np
import numpy.typing as npt
import pyarrow as pa
import pyarrow.compute as pacomp
from typing_extensions import Self

from ..._exception import SOMAError


def _expand_compressed_index_pointers(indptr: npt.NDArray, nval: int) -> npt.NDArray:
    """Expand CSR or CSC index pointers into COO indices.

    Args:
        indptr: Array of row/column pointers to be expanded.
        nval: Number of total non-zero values in the sparse matrix.
    """
    indices = np.empty(nval, dtype=np.int64)
    for index in range(len(indptr) - 1):
        indices[indptr[index] : indptr[index + 1]] = index
    return indices


def _unique_ptr(indptr: npt.NDArray) -> npt.NDArray:
    """Returns the unit indices for CSR or CSC points.

    Args:
        indptr: Array of row/column pointers.
    """
    return np.array(
        [
            index
            for index, (prev_nval, total_nval) in enumerate(
                zip(indptr[:-1], indptr[1:])
            )
            if total_nval - prev_nval > 0
        ],
        dtype=np.int64,
    )


def _str_to_int(value: str) -> int:
    if not value.isdigit():
        raise ValueError("{value} is not an integer.")
    return int(value)


def _version_less_than(version: str, upper_bound: tuple[int, int, int]) -> bool:
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
    gene_expression_path: str | Path,
) -> tuple[int, int, int]:
    with TenXCountMatrixReader(gene_expression_path) as reader:
        version = reader.software_version
    return version


class TenXCountMatrixReader:
    """Reader for the SpaceRanger, CellRanger, and XeniumRanger HDF5 count matrix
    group.

    The matrix is stored with features (var) on the rows and spot-barcodes (obs) on the
    columns. The data is not loaded until the user either explicitly loads all data with
    ``load`` or calls for a particular piece of data.

    The following information should be provided in the provided HDF5 file:

    HDF5 Metadata (called attributes):

    Matrix: Expression matrix stored in CSC format
      * /matrix/data: (int64) Count of feature per barcode stored.
      * /matrix/indices: (int64) Indices for the features.
      * /matrix/indptr: (int64) Compressed indices for the barcodes.
      * /matrix/shape: (int32) Length 2 array with the shape of the matrix.

    Observation Data:
      * /matrix/barcodes: Barcodes for obs data.

    Variable Data:
      * /matrix/features/feature_type: The type of feature reference the feature
        belongs to (for example "Gene Expression" or "Antibody Capture".
      * /matrix/features/id: Unique idenitier for feature.
      * /matrix/features/name: Human-readable name for variable. Used as
        `var_id` in scanpy.
      * /matrix/feature/genome: The genome reference for each feature.

    """

    def __init__(self, input_path: str | Path):
        # File management.
        self._path = input_path
        self._root: h5py.File | None = None

        # Metadatata.
        self._software_name: str | None = None
        self._software_version: tuple[int, int, int] | None = None
        self._nobs: int | None = None
        self._nvar: int | None = None

        # X matrix.
        self._data: npt.NDArray | None = None  # matrix/data
        self._feature_indices: npt.NDArray | None = None  # matrix/indices
        self._barcode_indptr: npt.NDArray | None = None  # matrix/indptr
        self._barcode_indices: npt.NDArray | None = None

        # Obs values.
        self._barcodes: npt.NDArray | None = None  # matrix/barcodes

        # Var values.
        self._var_name: npt.NDArray | None = None  # features/name
        self._gene_id: npt.NDArray | None = None  # features/id (unique identifier)
        self._feature_type: npt.NDArray | None = None  # features/feature_type
        self._genome: npt.NDArray | None = None  # features/genome

    def __enter__(self) -> Self:
        self.open()
        return self

    def __exit__(self, *_: Any) -> None:
        self.close()

    def _read_software_version(self) -> tuple[int, int, int]:
        if self._root is None:
            raise RuntimeError(f"[internal] File '{self._path}' is not open reading.")
        try:
            raw_version = self._root.attrs["software_version"]
        except KeyError as ke:
            raise SOMAError(
                f"Unable to read software version from gene expression file "
                f"{self._path}."
            ) from ke
        if not isinstance(raw_version, str):
            raise SOMAError(
                f"Unexpected type {type(raw_version)!r} for software version in gene "
                f"expression file {self._path}. Expected a string."
            )
        version: str | list[str] = raw_version.split("-")
        if len(version) == 1:
            version = version[0]
        elif len(version) == 2:
            self._software = version[0]
            version = version[1]
        else:
            raise SOMAError(
                f"Unexpected value {raw_version} for 'software_version' in gene "
                f"expression file {self._path}."
            )

        version = version.split(".")
        if len(version) not in {3, 4}:
            raise SOMAError(
                f"Unexpected value {raw_version} for 'software_version' in gene "
                f"expression file {self._path}. Expected a version in the form "
                f"'software_name-major.minor.patch'."
            )
        try:
            major = _str_to_int(version[0])
            minor = _str_to_int(version[1])
            patch = _str_to_int(version[2])
        except ValueError as err:
            raise SOMAError(
                f"Unexpected value {raw_version} for 'software_version' in gene "
                f"expression file {self._path}. Expected a version in the form "
                f"'software_name-major.minor.patch'."
            ) from err
        return (major, minor, patch)

    def close(self) -> None:
        if self._root is not None:
            self._root.close()
            self._root = None

    @property
    def data(self) -> pa.Array:
        if self._data is None:
            self._data = self.matrix_group["data"][()]
        return pa.array(self._data)

    @property
    def feature_group(self) -> h5py.Group:
        """Returns the /matrix group in the HDF5 file."""
        if self._root is None:
            raise RuntimeError(f"[internal] File '{self._path}' is not open reading.")
        return self._root["matrix"]["features"]

    @property
    def feature_type(self) -> pa.Array:
        if self._feature_type is None:
            self._feature_type = self.feature_group["feature_type"][()].astype("str")
        return pa.array(self._feature_type)

    @property
    def gene_id(self) -> pa.Array:
        if self._gene_id is None:
            self._gene_id = self.feature_group["id"][()].astype(str)
        return pa.array(self._gene_id)

    @property
    def genome(self) -> pa.Array:
        if self._genome is None:
            self._genome = self.feature_group["genome"][()].astype(str)
        return pa.array(self._genome)

    @property
    def matrix_group(self) -> h5py.Group:
        """Returns the /matrix group in the HDF5 file."""
        if self._root is None:
            raise RuntimeError(f"[internal] File '{self._path}' is not open reading.")
        return self._root["matrix"]

    @property
    def obs_id(self) -> pa.Array:
        if self._barcodes is None:
            self._barcodes = self.matrix_group["barcodes"][()].astype(str)
        return pa.array(self._barcodes)

    @property
    def obs_indices(self) -> pa.Array:
        if self._barcode_indices is None:
            self._barcode_indptr = self.matrix_group["indptr"][()]
            nvalues = self.matrix_group["data"].size
            self._barcode_indices = _expand_compressed_index_pointers(
                self._barcode_indptr, nvalues
            )
        return pa.array(self._barcode_indices)

    def open(self) -> None:
        if self._root is None:
            self._root = h5py.File(self._path, "r")

    def load(self) -> None:
        """Load all data used for SOMA.

        Note: The file version is not loaded.
        """
        # Groups
        if self._root is None:
            raise RuntimeError("Cannot load data. Open reader first.")
        matrix_group = self._root["matrix"]
        feature_group = matrix_group["features"]

        # Shape (convert from np.int32 -> int).
        shape = matrix_group["shape"]
        self._nvar = int(shape[0])
        self._nobs = int(shape[1])

        # X matrix
        self._data = matrix_group["data"][()]
        self._feature_indices = matrix_group["indices"][()]
        self._barcode_indptr = matrix_group["indptr"][()]
        self._barcode_indices = _expand_compressed_index_pointers(
            self._barcode_indptr, self._data.size
        )

        # obs data
        self._barcodes = matrix_group["barcodes"][()].astype(str)

        # var data
        self._var_name = feature_group["name"][()].astype(str)
        self._gene_id = feature_group["id"][()].astype(str)
        self._feature_type = feature_group["feature_type"][()].astype(str)
        self._genome = feature_group["genome"][()].astype(str)

    @property
    def nobs(self) -> int:
        if self._nobs is None:
            self._nobs = int(self.matrix_group["shape"][1])
        return self._nobs

    @property
    def nvar(self) -> int:
        if self._nvar is None:
            self._nvar = int(self.matrix_group["shape"][0])
        return self._nvar

    @property
    def software_version(self) -> tuple[int, int, int]:
        if self._software_version is None:
            self._software_version = self._read_software_version()
        return self._software_version

    def unique_obs_indices(self) -> pa.Array:
        """Returns the unique obs indices that have non-zero values in the X matrix."""
        if self._barcode_indptr is None:
            self._barcode_indptr = self.matrix_group["indptr"][()]
        indices = _unique_ptr(self._barcode_indptr)
        return pa.array(indices)

    def unique_var_indices(self) -> pa.Array:
        """Returns the unique var indices that have non-zero values in the X matrix."""
        var_indices = self.var_indices
        if not isinstance(self.var_indices, pa.Array):
            # This check is required for PyArrow 11.0.0. Otherwise, we get the
            # error:
            # E   pyarrow.lib.ArrowInvalid: Could not convert
            # <pyarrow.Int64Scalar: 0> with type pyarrow.lib.Int64Scalar: did
            # not recognize Python value type when inferring an Arrow data type
            var_indices = pa.array(var_indices)
        return pacomp.unique(var_indices)

    @property
    def var_id(self) -> pa.Array:
        if self._var_name is None:
            self._var_name = self.feature_group["name"][()].astype(str)
        return pa.array(self._var_name)

    @property
    def var_indices(self) -> pa.Array:
        if self._feature_indices is None:
            self._feature_indices = self.matrix_group["indices"][()]
        return pa.array(self._feature_indices)
