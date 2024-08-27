import abc
import json
from typing import Any, Dict, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt
import pyarrow as pa
from somacore import coordinates
from typing_extensions import Self


class Axis(coordinates.Axis):
    """A description of an axis of a coordinate system

    TODO: Note if this class remains more or less as is the base class in somacore
    can be implemented as a ``dataclasses.dataclass``.

    Lifecycle: experimental
    """

    @classmethod
    def from_json(cls, data: str) -> Self:
        """Create from a json blob.

        Args:
           data: json blob to deserialize.

        Lifecycle: experimental
        """
        kwargs = json.loads(data)
        if not isinstance(kwargs, dict):
            raise ValueError()
        return cls(**kwargs)

    def to_dict(self) -> Dict[str, Any]:
        """TODO: Add docstring"""
        kwargs: Dict[str, Any] = {"name": self.name}
        if self.units is not None:
            kwargs["units"] = self.units
        if self.scale is not None:
            kwargs["scale"] = self.scale
        return kwargs

    def to_json(self) -> str:
        """TODO: Add docstring"""
        return json.dumps(self.to_dict())


class CoordinateSpace(coordinates.CoordinateSpace):
    """A coordinate system for spatial data."""

    @classmethod
    def from_json(cls, data: str) -> Self:
        """TODO: Add docstring"""
        # TODO: Needs good, comprehensive error handling.
        raw = json.loads(data)
        return cls(tuple(Axis(**axis) for axis in raw))

    def __init__(self, axes: Sequence[Axis]):
        """TODO: Add docstring"""
        # TODO: Needs good, comprehensive error handling.
        if len(tuple(axes)) == 0:
            raise ValueError("Coordinate space must have at least one axis.")
        self._axes = tuple(axes)

    def __len__(self) -> int:
        return len(self._axes)

    def __getitem__(self, index: int) -> Axis:
        return self._axes[index]

    def __repr__(self) -> str:
        output = f"<{type(self).__name__}\n"
        for axis in self._axes:
            output += f"  {axis},\n"
        return output + ">"

    @property
    def axes(self) -> Tuple[Axis, ...]:
        """TODO: Add docstring"""
        return self._axes

    @property
    def axis_names(self) -> Tuple[str, ...]:
        return tuple(axis.name for axis in self._axes)

    def rank(self) -> int:
        return len(self)

    def to_json(self) -> str:
        """TODO: Add docstring"""
        return json.dumps(tuple(axis.to_dict() for axis in self._axes))


class CoordinateTransform(coordinates.CoordinateTransform):
    def __init__(
        self,
        input_axes: Union[str, Sequence[str]],
        output_axes: Union[str, Sequence[str]],
    ):
        self._input_axes = (
            (input_axes,) if isinstance(input_axes, str) else tuple(input_axes)
        )
        self._output_axes = (
            (output_axes,) if isinstance(output_axes, str) else tuple(output_axes)
        )

    @abc.abstractmethod
    def __mul__(self, other: Any) -> "CoordinateTransform":
        raise NotImplementedError()

    @abc.abstractmethod
    def __rmul__(self, other: Any) -> "CoordinateTransform":
        raise NotImplementedError()

    @abc.abstractmethod
    def apply(self, data: Union[pa.Tensor, pa.Table]) -> Union[pa.Tensor, pa.Table]:
        raise NotImplementedError()

    @property
    def input_axes(self) -> Tuple[str, ...]:
        return self._input_axes

    @property
    def input_rank(self) -> int:
        return len(self._input_axes)

    @property
    def output_axes(self) -> Tuple[str, ...]:
        return self._output_axes

    @property
    def output_rank(self) -> int:
        return len(self._output_axes)

    def to_dict(self) -> Dict[str, Any]:
        return {"input_axes": self.input_axes, "output_axes": self.output_axes}


class IdentityCoordinateTransform(CoordinateTransform):
    """TODO: Add docstring"""

    def __init__(
        self,
        input_axes: Union[str, Sequence[str]],
        output_axes: Union[str, Sequence[str]],
    ):
        super().__init__(input_axes, output_axes)
        if self.input_rank != self.output_rank:
            raise ValueError("Incompatible rank of input and output axes")

    def __mul__(self, other: Any) -> CoordinateTransform:
        if np.isscalar(other):
            raise NotImplementedError(
                "Support for multiplying by scalars is not yet implemented."
            )
        if isinstance(other, CoordinateTransform):
            if isinstance(other, IdentityCoordinateTransform):
                if self.output_axes != other.input_axes:
                    raise ValueError("Axis mismatch between transformations.")
                return IdentityCoordinateTransform(self.input_axes, other.output_axes)
            return other.__rmul__(self)
        if isinstance(other, np.ndarray):
            raise NotImplementedError(
                "Support for multiplying by numpy arrays is not yet implemented."
            )
        raise TypeError(
            f"Cannot multiply a CoordinateTransform by type {type(other)!r}."
        )

    def __rmul__(self, other: Any) -> CoordinateTransform:
        if np.isscalar(other):
            return self.__mul__(other)
        if isinstance(other, CoordinateTransform):
            if isinstance(other, IdentityCoordinateTransform):
                if other.output_axes != self.input_axes:
                    raise ValueError("Axis mismatch between transformations.")
                return IdentityCoordinateTransform(other.input_axes, self.output_axes)
            return other.__mul__(self)
        if isinstance(other, np.ndarray):
            raise NotImplementedError(
                "Support for multiplying by numpy arrays is not yet implemented."
            )
        raise TypeError(
            f"Cannot multiply a CoordinateTransform by type {type(other)!r}."
        )

    def apply(self, data: Union[pa.Tensor, pa.Table]) -> Union[pa.Tensor, pa.Table]:
        # TODO: Check valid rank
        return data

    def to_dict(self) -> Dict[str, Any]:
        return super().to_dict()


class AffineCoordinateTransform(CoordinateTransform):
    """TODO: Add docstring"""

    def __init__(
        self,
        input_axes: Union[str, Sequence[str]],
        output_axes: Union[str, Sequence[str]],
        matrix: npt.ArrayLike,
    ):
        super().__init__(input_axes, output_axes)

        # Check the rank of the input/output axes match.
        if self.input_rank != self.output_rank:
            raise ValueError(
                "The input axes and output axes must be the same length for an "
                "affine transform."
            )
        rank = self.input_rank

        # Create and validate the augmented matrix.
        self._matrix: npt.NDArray[np.float64] = np.array(matrix, dtype=np.float64)
        if self._matrix.shape == (rank + 1, rank + 1):
            if not (
                self._matrix[-1, -1] == 1.0
                and np.array_equal(self._matrix[:-1, -1], np.zeros((rank,)))
            ):
                raise ValueError(
                    "Input matrix has augmented matrix shape, but is not a valid "
                    "augmented matrix."
                )
        elif self._matrix.shape == (rank, rank + 1):
            self._matrix = np.vstack(
                (
                    self._matrix,
                    np.hstack((np.zeros((rank,)), np.array([1]))),
                )
            )
        else:
            raise ValueError(
                f"Unexpected shape {self._matrix.shape} for the input affine matrix."
            )

    def __mul__(self, other: Any) -> CoordinateTransform:
        if np.isscalar(other):
            return AffineCoordinateTransform(
                self.input_axes,
                self.output_axes,
                other * self.augmented_matrix,  # type: ignore
            )
        if isinstance(other, CoordinateTransform):
            if self.output_axes != other.input_axes:
                raise ValueError("Axis mismatch between transformations.")
            if isinstance(other, IdentityCoordinateTransform):
                return AffineCoordinateTransform(
                    self.input_axes, other.output_axes, self._matrix
                )
            if isinstance(other, AffineCoordinateTransform):
                return AffineCoordinateTransform(
                    self.input_axes, other.output_axes, self._matrix @ other._matrix
                )
        if isinstance(other, np.ndarray):
            raise NotImplementedError(
                "Support for multiplying by numpy arrays is not yet implemented."
            )
        raise TypeError(
            f"Cannot multiply a CoordinateTransform by type {type(other)!r}."
        )

    def __rmul__(self, other: Any) -> CoordinateTransform:
        if np.isscalar(other):
            return self.__mul__(other)
        if isinstance(other, CoordinateTransform):
            if other.output_axes != self.input_axes:
                raise ValueError("Axis mismatch between transformations.")
            if isinstance(other, IdentityCoordinateTransform):
                return AffineCoordinateTransform(
                    other.input_axes, self.output_axes, self._matrix
                )
            if isinstance(other, AffineCoordinateTransform):
                return AffineCoordinateTransform(
                    other.input_axes, self.output_axes, other._matrix @ self._matrix
                )
        if isinstance(other, np.ndarray):
            raise NotImplementedError(
                "Support for multiplying by numpy arrays is not yet implemented."
            )
        raise TypeError(
            f"Cannot multiply a CoordinateTransform by type {type(other)!r}."
        )

    def apply(self, data: Union[pa.Tensor, pa.Table]) -> Union[pa.Tensor, pa.Table]:
        """TODO: Add docstring"""
        raise NotImplementedError()

    @property
    def augmented_matrix(self) -> npt.NDArray[np.float64]:
        """Returns the augmented affine matrix for the transformation."""
        return self._matrix

    def to_dict(self) -> Dict[str, Any]:
        kwargs = super().to_dict()
        kwargs["matrix"] = self._matrix.tolist()
        return kwargs


def transform_from_json(data: str) -> CoordinateTransform:
    """TODO: Add docstring"""
    raw = json.loads(data)
    try:
        transform_type = raw.pop("transform_type")
    except KeyError:
        raise KeyError()  # TODO Add error message
    try:
        kwargs = raw.pop("transform")
    except KeyError:
        raise KeyError()  # TODO Add error message
    if transform_type == "IdentityCoordinateTransform":
        return IdentityCoordinateTransform(**kwargs)
    elif transform_type == "AffineCoordinateTransform":
        return AffineCoordinateTransform(**kwargs)
    else:
        raise KeyError("Unrecognized transform type key 'transform_type'")


def transform_to_json(transform: CoordinateTransform) -> str:
    kwargs = transform.to_dict()
    transform_type = type(transform).__name__
    return json.dumps({"transform_type": transform_type, "transform": kwargs})
