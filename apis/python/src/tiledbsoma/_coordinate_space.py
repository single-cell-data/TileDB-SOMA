# Copyright (c) TileDB, Inc. and The Chan Zuckerberg Initiative Foundation
#
# Licensed under the MIT License.
"""Definitions of types related to coordinate systems."""

from __future__ import annotations

import abc
import collections.abc
import itertools
from collections.abc import Iterable, Sequence
from typing import Union

import attrs
import numpy as np
import numpy.typing as npt
from typing_extensions import Self

from ._types import str_or_seq_length, to_string_tuple


@attrs.define(frozen=True)
class Axis:
    """A description of an axis of a coordinate system.

    Lifecycle: experimental
    """

    name: str = attrs.field()
    """Name of the axis."""
    unit: Union[str, None] = None  # noqa: UP007
    """Optional string name for the units of the axis."""

    @name.validator
    def check(self, _, value: str) -> None:  # type: ignore[no-untyped-def]  # noqa: ANN001
        if value.startswith("soma_"):
            raise ValueError(f"Invalid axis name '{value}'. Cannot start with 'soma_'.")


@attrs.define(frozen=True)
class CoordinateSpace(collections.abc.Sequence[Axis]):
    """A coordinate space for spatial data.

    Args:
        axes: The axes of the coordinate system in order.

    Lifecycle: experimental
    """

    axes: tuple[Axis, ...] = attrs.field(converter=tuple)

    @classmethod
    def from_axis_names(cls, axis_names: Sequence[str,]) -> Self:
        return cls(tuple(Axis(name) for name in axis_names))  # type: ignore[misc]

    @axes.validator
    def _validate(self, _, axes: tuple[Axis, ...]) -> None:  # type: ignore[no-untyped-def]  # noqa: ANN001
        if not axes:
            raise ValueError("The coordinate space must have at least one axis.")
        if len({axis.name for axis in self.axes}) != len(axes):
            raise ValueError("The names for the axes must be unique.")

    def __len__(self) -> int:
        return len(self.axes)

    def __getitem__(self, index: int) -> Axis:  # type: ignore[override]
        return self.axes[index]

    @property
    def axis_names(self) -> tuple[str, ...]:
        """The names of the axes in order.

        Lifecycle: experimental
        """
        return tuple(axis.name for axis in self.axes)


class CoordinateTransform(metaclass=abc.ABCMeta):
    """A coordinate transformation from one coordinate space to another.

    Args:
        input_axes: The names of the axes for the input coordinate space.
        output_axes: The names of the axes for the output coordinate space.

    CoordinateTransform classes are composable using the ``@`` (__matmul__) operator.

    Lifecycle: experimental
    """

    def __init__(
        self,
        input_axes: str | Sequence[str],
        output_axes: str | Sequence[str],
    ) -> None:
        self._input_axes = to_string_tuple(input_axes)
        self._output_axes = to_string_tuple(output_axes)

    def _check_matmul_inner_axes(self, other: CoordinateTransform) -> None:
        """Throws a ``ValueError`` if ``self @ other`` has mismatched axes."""
        if self.input_axes != other.output_axes:
            raise ValueError(f"Input axes of {type(self).__name__} must match output axes of {type(other).__name__}.")

    def _check_rmatmul_inner_axes(self, other: CoordinateTransform) -> None:
        """Throws a ``ValueError `` if ``other @ self`` has mismatched axes."""
        if self.output_axes != other.input_axes:
            raise ValueError(f"Input axes of {type(other).__name__} must match output axes of {type(self).__name__}.")

    @abc.abstractmethod
    def _contents_lines(self) -> Iterable[str]:
        return
        yield

    def _my_repr(self) -> Iterable[str]:
        yield f"{type(self).__name__}"
        yield f"  input axes: {self._input_axes}"
        yield f"  output axes: {self._output_axes}"

    def __repr__(self) -> str:
        content = self._contents_lines
        lines = self._my_repr() if content is None else itertools.chain(self._my_repr(), self._contents_lines())
        return "<" + "\n".join(lines) + ">"

    @abc.abstractmethod
    def __matmul__(self, other: object) -> CoordinateTransform:
        raise NotImplementedError

    @abc.abstractmethod
    def inverse_transform(self) -> CoordinateTransform:
        """Returns the inverse coordinate transform.

        Lifecycle: experimental
        """
        raise NotImplementedError

    @property
    def input_axes(self) -> tuple[str, ...]:
        """The names of the axes of the input coordinate space.

        Lifecycle: experimental
        """
        return self._input_axes

    @property
    def output_axes(self) -> tuple[str, ...]:
        """The names of the axes of the output coordinate space.

        Lifecycle: experimental
        """
        return self._output_axes


class AffineTransform(CoordinateTransform):
    """An affine coordinate trasformation from one coordinate space to another.

    An affine transform is a combination of a linear transformation and a translation.

    Args:
        input_axes: The names of the axes for the input coordinate space.
        output_axes: The names of the axes for the output coordinate space.
        matrix: Matrix (perhaps augmented) that represents the affine transformation.
            Can be provided as just the linear transform (if no translation), the
            full augmented matrix, or the augmented matrix without the final row.

    Lifecycle: experimental
    """

    def __init__(
        self,
        input_axes: str | Sequence[str],
        output_axes: str | Sequence[str],
        matrix: npt.ArrayLike | Sequence[int | float],
    ) -> None:
        super().__init__(input_axes, output_axes)

        # Check the rank of the input/output axes match.
        if len(self.input_axes) != len(self.output_axes):
            raise ValueError("The input axes and output axes must be the same length for an affine transform.")
        rank = len(self.input_axes)

        # Create and validate the augmented matrix.
        self._matrix: npt.NDArray[np.float64] = np.array(matrix, dtype=np.float64)
        if self._matrix.shape == (rank + 1, rank + 1):
            if not (self._matrix[-1, -1] == 1.0 and np.array_equal(self._matrix[-1, :-1], np.zeros((rank,)))):
                raise ValueError(
                    f"Input matrix {self._matrix} has augmented matrix shape, but is not a valid augmented matrix."
                )
        elif self._matrix.shape == (rank, rank + 1):
            self._matrix = np.vstack((
                self._matrix,
                np.hstack((np.zeros((rank,)), np.array([1]))),
            ))
        elif self._matrix.shape == (rank, rank):
            self._matrix = np.vstack((
                np.hstack((self._matrix, np.zeros((rank, 1)))),
                np.hstack((np.zeros((rank,)), np.array([1]))),
            ))
        else:
            raise ValueError(f"Unexpected shape {self._matrix.shape} for the input affine matrix.")

    def _contents_lines(self) -> Iterable[str]:
        yield "  augmented matrix:"
        yield "    " + str(self._matrix).replace("\n", "\n    ")

    def __matmul__(self, other: object) -> CoordinateTransform:
        if not isinstance(other, CoordinateTransform):
            raise NotImplementedError(f"Matrix multiply is not implemented with type {type(other)!r}.")
        self._check_matmul_inner_axes(other)
        if isinstance(other, IdentityTransform):
            return AffineTransform(other.input_axes, self.output_axes, self._matrix)
        if isinstance(other, AffineTransform):
            return AffineTransform(
                other.input_axes,
                self.output_axes,
                self.augmented_matrix @ other.augmented_matrix,
            )
        raise NotImplementedError(f"Cannot multiply a CoordinateTransform by type {type(other)!r}.")

    @property
    def augmented_matrix(self) -> npt.NDArray[np.float64]:
        """Returns the augmented affine matrix for the transformation.

        Lifecycle: experimental
        """
        return self._matrix

    def inverse_transform(self) -> AffineTransform:
        """Returns the inverse coordinate transform.

        Lifecycle: experimental
        """
        rank = len(self.output_axes)
        inv_a = np.linalg.inv(self._matrix[:-1, :-1])
        b2 = -inv_a @ self._matrix[:-1, -1].reshape((rank, 1))
        inv_augmented: npt.NDArray[np.float64] = np.vstack((
            np.hstack((inv_a, b2)),
            np.hstack((np.zeros(rank), np.array([1]))),
        ))
        return AffineTransform(self.output_axes, self.input_axes, inv_augmented)


class ScaleTransform(AffineTransform):
    """A scale coordinate transformation from one coordinate space to another.

    Args:
        input_axes: The names of the axes for the input coordinate space.
        output_axes: The names of the axes for the output coordinate space.
        scale_factors: The scale factors for the transformation. There must be one
            value per axis.

    Lifecycle: experimental
    """

    def __init__(
        self,
        input_axes: str | Sequence[str],
        output_axes: str | Sequence[str],
        scale_factors: npt.ArrayLike | Sequence[int | float],
    ) -> None:
        rank = str_or_seq_length(input_axes)
        self._scale_factors: npt.NDArray[np.float64] = np.array(scale_factors, dtype=np.float64)
        if self._scale_factors.size != rank:
            raise ValueError(
                f"Scale factors have unexpected shape={self._scale_factors.shape} for a transform with rank={rank}."
            )
        self._scale_factors = self._scale_factors.reshape((rank,))

        super().__init__(input_axes, output_axes, np.diag(self._scale_factors))

    def _contents_lines(self) -> Iterable[str]:
        yield f"  scales: {self._scale_factors}"

    def __matmul__(self, other: object) -> CoordinateTransform:
        if not isinstance(other, CoordinateTransform):
            raise NotImplementedError(f"Matrix multiply is not implemented with type {type(other)!r}.")
        self._check_matmul_inner_axes(other)
        if isinstance(other, ScaleTransform):
            return ScaleTransform(
                other.input_axes,
                self.output_axes,
                self.scale_factors * other.scale_factors,
            )
        return super().__matmul__(other)

    def inverse_transform(self) -> ScaleTransform:
        """Returns the inverse coordinate transform.

        Lifecycle: experimental
        """
        return ScaleTransform(self.output_axes, self.input_axes, 1.0 / self._scale_factors)

    @property
    def scale_factors(self) -> npt.NDArray[np.float64]:
        """Returns the scale factors as an one-dimensional numpy array.

        Lifecycle: experimental
        """
        return self._scale_factors


class UniformScaleTransform(ScaleTransform):
    """A scale coordinate transformation from one coordinate space to another.

    Args:
        input_axes: The names of the axes for the input coordinate space.
        output_axes: The names of the axes for the output coordinate space.
        scale: The scale factor for all axes.

    Lifecycle: experimental
    """

    def __init__(
        self,
        input_axes: str | Sequence[str],
        output_axes: str | Sequence[str],
        scale: int | float | np.float64,
    ) -> None:
        self._scale = float(scale)
        rank = str_or_seq_length(input_axes)
        super().__init__(input_axes, output_axes, rank * [self._scale])

    def _contents_lines(self) -> Iterable[str]:
        yield f"  scale: {self._scale}"

    def __matmul__(self, other: object) -> CoordinateTransform:
        if not isinstance(other, CoordinateTransform):
            raise NotImplementedError(f"Matrix multiply is not implemented with type {type(other)!r}.")
        if isinstance(other, UniformScaleTransform):
            self._check_matmul_inner_axes(other)
            return UniformScaleTransform(other.input_axes, self.output_axes, self.scale * other.scale)
        return super().__matmul__(other)

    def inverse_transform(self) -> UniformScaleTransform:
        """Returns the inverse coordinate transform.

        Lifecycle: experimental
        """
        return UniformScaleTransform(self.output_axes, self.input_axes, 1.0 / self._scale)

    @property
    def scale(self) -> float:
        """Returns the scale factor for the uniform scale transform.

        Lifecycle: experimental
        """
        return self._scale


class IdentityTransform(UniformScaleTransform):
    """The identify transform from one coordinate space to another.

    This transform only changes the name of the axes.

    Args:
        input_axes: The names of the axes for the input coordinate space.
        output_axes: The names of the axes for the output coordinate space.

    Lifecycle: experimental
    """

    def __init__(
        self,
        input_axes: str | Sequence[str],
        output_axes: str | Sequence[str],
    ) -> None:
        super().__init__(input_axes, output_axes, 1)

    def _contents_lines(self) -> Iterable[str]:
        return
        yield

    def __matmul__(self, other: object) -> CoordinateTransform:
        if not isinstance(other, CoordinateTransform):
            raise NotImplementedError(f"Matrix multiply is not implemented with type {type(other)!r}.")
        if isinstance(other, IdentityTransform):
            self._check_matmul_inner_axes(other)
            return IdentityTransform(other.input_axes, self.output_axes)
        return super().__matmul__(other)

    def inverse_transform(self) -> IdentityTransform:
        """Returns the inverse coordinate transform.

        Lifecycle: experimental
        """
        return IdentityTransform(self.output_axes, self.input_axes)
