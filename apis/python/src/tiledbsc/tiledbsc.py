from __future__ import annotations
from collections.abc import Sequence
from typing import Optional

from dataclasses import dataclass
import numpy as np
import os

import tiledb


def create(uri: str, groups=None, anndata=None):
    """
    Create a TileDB SC Dataset
    """
    pass


def open(uri: str, mode: str = "r", timestamp: Optional[int] = None):
    """
    Open a TileDB SC Dataset or Group
    """
    pass


def close():
    """
    Close a TileDB SC Dataset or Group
    """
    pass


class Group:
    def __init__(
        self,
        obs: str = None,
        obs_index: str = None,
        var: str = None,
        var_index: str = None,
        X: dict = None,
        aux=None,
        meta=None,
    ):
        self._obs = obs
        self._obs_index = obs_index
        self._var = var
        self._var_index = var_index
        self._X = X
        self._aux = aux
        self._meta = meta

    @property
    def obs(self):
        return self._obs

    @property
    def var(self):
        return self._var

    @property
    def X(self) -> dict:
        return self._X

    @property
    def groups(self) -> dict:
        return self._X

    def _X_schema(self, data_dtype=np.float32):
        # TODO TBD layers need multiple attributes
        return tiledb.ArraySchema(
            domain=tiledb.Domain(
                tiledb.Dim(name=self._obs_index, dtype="ascii"),
                tiledb.Dim(name=self._var_index, dtype="ascii"),
            ),
            attrs=[tiledb.Attr(name="data", dtype=data_dtype)],
            sparse=True,
        )

    def _obs_schema(self, data_dtype=np.float32):
        return tiledb.ArraySchema(
            domain=tiledb.Domain(
                tiledb.Dim(name="dim", dtype="ascii"),
            ),
            attrs=[tiledb.Attr(name=self._obs_index, dtype=data_dtype)],
            sparse=True,
        )

    def _var_schema(self, data_dtype=np.float32):
        return tiledb.ArraySchema(
            domain=tiledb.Domain(
                tiledb.Dim(name="dim", dtype="ascii"),
            ),
            attrs=[tiledb.Attr(name=self._var_index, dtype=data_dtype)],
            sparse=True,
        )

    def _create(self, base: str, group_name: str):
        group_uri = os.path.join(base, group_name)
        tiledb.group_create(group_uri)

        tiledb.Array.create(os.path.join(group_uri, self._obs), self._obs_schema())
        tiledb.Array.create(os.path.join(group_uri, self._var), self._var_schema())

        # TODO incorrect. need to add multiple attributes per array as necessary
        for array_name in self._X:
            X_uri = os.path.join(group_uri, array_name)
            tiledb.Array.create(X_uri, self._X_schema())

    def query(self) -> Query:
        pass

    def __getitem__(self):
        pass

    def to_array(self):
        pass

    def to_dict(self):
        pass

    def to_df(self):
        pass

    def to_anndata(self):
        pass


class Dataset:
    def __init__(self, uri, groups):
        tiledb.group_create(uri)

        for group_name in groups:
            groups[group_name]._create(uri, group_name)
