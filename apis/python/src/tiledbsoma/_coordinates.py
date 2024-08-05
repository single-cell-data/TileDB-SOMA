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
        self._axes = tuple(axes)

    def __len__(self) -> int:
        return len(self._axes)

    def __getitem__(self, index: int) -> Axis:
        return self._axes[index]

    @property
    def axes(self) -> Tuple[Axis, ...]:
        """TODO: Add docstring"""
        return self._axes

    def to_json(self) -> str:
        """TODO: Add docstring"""
        return json.dumps(tuple(axis.to_dict() for axis in self._axes))


class BaseCoordinateTransform(coordinates.CoordinateTransform):
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


class IdentityCoordinateTransform(BaseCoordinateTransform):
    """TODO: Add docstring"""

    def __init__(
        self,
        input_axes: Union[str, Sequence[str]],
        output_axes: Union[str, Sequence[str]],
    ):
        super().__init__(input_axes, output_axes)
        if self.input_rank != self.output_rank:
            raise ValueError("Incompatible rank of input and output axes")

    def apply(self, data: Union[pa.Tensor, pa.Table]) -> Union[pa.Tensor, pa.Table]:
        # TODO: Check valid rank
        return data

    def to_dict(self) -> Dict[str, Any]:
        return super().to_dict()


class AffineCoordinateTransform(BaseCoordinateTransform):
    """TODO: Add docstring"""

    def __init__(
        self,
        input_axes: Union[str, Sequence[str]],
        output_axes: Union[str, Sequence[str]],
        affine_matrix: npt.NDArray[np.float64],
    ):
        super().__init__(input_axes, output_axes)
        # TODO: Check the rank of the affine matrix
        self._affine_matrix = affine_matrix

    def apply(self, data: Union[pa.Tensor, pa.Table]) -> Union[pa.Tensor, pa.Table]:
        """TODO: Add docstring"""
        raise NotImplementedError()

    def to_dict(self) -> Dict[str, Any]:
        kwargs = super().to_dict()
        kwargs["affine_matrix_data"] = self._affine_matrix.tolist()
        kwargs["affine_matrix_shape"] = self._affine_matrix.shape
        return kwargs


def transform_from_json(data: str) -> coordinates.CoordinateTransform:
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
        affine_matrix_data = kwargs.pop("affine_matrix_data")
        affine_matrix_shape = kwargs.pop("affine_matrix_shape")
        return AffineCoordinateTransform(
            affine_matrix=np.array(affine_matrix_data).reshape(affine_matrix_shape),
            **kwargs,
        )
    else:
        raise KeyError("Unrecognized transform type key 'transform_type'")


def transform_to_json(transform: BaseCoordinateTransform) -> str:
    kwargs = transform.to_dict()
    transform_type = type(transform).__name__
    return json.dumps({"transform_type": transform_type, "transform": kwargs})
