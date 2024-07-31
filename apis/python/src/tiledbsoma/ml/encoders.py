# Copyright (c) 2021-2024 The Chan Zuckerberg Initiative Foundation
# Copyright (c) 2021-2024 TileDB, Inc.
#
# Licensed under the MIT License.

from __future__ import annotations

import abc
import functools
from typing import List

import numpy.typing as npt
import pandas as pd
from sklearn.preprocessing import LabelEncoder as SklearnLabelEncoder


class Encoder(abc.ABC):
    """Base class for ``obs`` encoders.

    To define a custom encoder, five methods must be implemented:

    - ``fit``: defines how the encoder will be fitted to the data.
    - ``transform``: defines how the encoder will be applied to the data
      in order to create an ``obs`` tensor.
    - ``inverse_transform``: defines how to decode the encoded values back
      to the original values.
    - ``name``: The name of the encoder. This will be used as the key in the
      dictionary of encoders. Each encoder passed to a :class:`.pytorch.ExperimentDataPipe` must have a unique name.
    - ``columns``: List of columns in ``obs`` that the encoder will be applied to.

    See the implementation of :class:`LabelEncoder` for an example.
    """

    @abc.abstractmethod
    def fit(self, obs: pd.DataFrame) -> None:
        """Fit the encoder with obs."""
        pass

    @abc.abstractmethod
    def transform(self, df: pd.DataFrame) -> npt.ArrayLike:
        """Transform the obs :class:`pandas.DataFrame` into a :class:`pandas.DataFrame` of encoded values."""
        pass

    @abc.abstractmethod
    def inverse_transform(self, encoded_values: npt.ArrayLike) -> npt.ArrayLike:
        """Inverse transform the encoded values back to the original values."""
        pass

    @property
    @abc.abstractmethod
    def name(self) -> str:
        """Name of the encoder."""
        pass

    @property
    @abc.abstractmethod
    def columns(self) -> List[str]:
        """Columns in ``obs`` that the encoder will be applied to."""
        pass


class LabelEncoder(Encoder):
    """Default encoder based on :class:`sklearn.preprocessing.LabelEncoder`."""

    def __init__(self, col: str) -> None:
        self._encoder = SklearnLabelEncoder()
        self.col = col

    def fit(self, obs: pd.DataFrame) -> None:
        """Fit the encoder with ``obs``."""
        self._encoder.fit(obs[self.col].unique())

    def transform(self, df: pd.DataFrame) -> npt.ArrayLike:
        """Transform the obs :class:`pandas.DataFrame` into a :class:`numpy.typing.ArrayLike` of encoded values."""
        return self._encoder.transform(df[self.col])  # type: ignore

    def inverse_transform(self, encoded_values: npt.ArrayLike) -> npt.ArrayLike:
        """Inverse transform the encoded values back to the original values."""
        return self._encoder.inverse_transform(encoded_values)  # type: ignore

    @property
    def name(self) -> str:
        """Name of the encoder."""
        return self.col

    @property
    def columns(self) -> List[str]:
        """Columns in ``obs`` that the encoder will be applied to."""
        return [self.col]

    @property
    def classes_(self):  # type: ignore
        """Classes of the encoder."""
        return self._encoder.classes_


class BatchEncoder(Encoder):
    """An encoder that concatenates and encodes several ``obs`` columns."""

    def __init__(self, cols: List[str], name: str = "batch"):
        self.cols = cols
        from sklearn.preprocessing import LabelEncoder

        self._name = name
        self._encoder = LabelEncoder()

    def _join_cols(self, df: pd.DataFrame):  # type: ignore
        return functools.reduce(
            lambda a, b: a + b, [df[c].astype(str) for c in self.cols]
        )

    def transform(self, df: pd.DataFrame) -> npt.ArrayLike:
        """Transform the obs :class:`pandas.DataFrame` into a :class:`pandas.DataFrame` of encoded values."""
        arr = self._join_cols(df)
        return self._encoder.transform(arr)  # type: ignore

    def inverse_transform(self, encoded_values: npt.ArrayLike) -> npt.ArrayLike:
        """Inverse transform the encoded values back to the original values."""
        return self._encoder.inverse_transform(encoded_values)  # type: ignore

    def fit(self, obs: pd.DataFrame) -> None:
        """Fit the encoder with ``obs``."""
        arr = self._join_cols(obs)
        self._encoder.fit(arr.unique())

    @property
    def columns(self) -> List[str]:
        """Columns in ``obs`` that the encoder will be applied to."""
        return self.cols

    @property
    def name(self) -> str:
        """Name of the encoder."""
        return self._name

    @property
    def classes_(self):  # type: ignore
        """Classes of the encoder."""
        return self._encoder.classes_
