import abc
import json
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt
from somacore import coordinates
from typing_extensions import Self


class Axis(coordinates.Axis):
    """A description of an axis of a coordinate system

    TODO: Note if this class remains more or less as is the base class in somacore
    can be implemented as a ``dataclasses.dataclass``.

    Lifecycle: experimental
    """

    @classmethod
    def from_dict(cls, data: Mapping[str, str]) -> Self:
        """TODO: Docstring"""
        return cls(
            axis_name=data["name"],
            axis_type=data.get("type"),
            axis_unit=data.get("unit"),
        )

    @classmethod
    def from_json(cls, data: str) -> Self:
        """Create from a json blob.

        Args:
           data: json blob to deserialize.

        Lifecycle: experimental
        """
        kwargs = json.loads(data)
        if not isinstance(kwargs, dict):
            # TODO: Needs good, comprehensive error handling for expected metedata.
            raise ValueError()
        return cls.from_dict(kwargs)

    def __init__(
        self,
        *,
        axis_name: str,
        axis_type: Optional[str] = None,
        axis_unit: Optional[str] = None
    ):
        # TODO: Needs good, comprehensive error handling for valid input.
        self._name: str = axis_name
        self._type: Optional[str] = axis_type
        self._unit: Optional[str] = axis_unit

    def __eq__(self, other: Any) -> bool:
        """Equality test for two axes."""
        if not isinstance(other, Axis):
            return False
        return (
            self.name == other.name
            and self.type == other.type
            and self.unit == other.unit
        )

    @property
    def name(self) -> str:
        """TODO: Add docstring"""
        return self._name

    @property
    def type(self) -> Optional[str]:
        """TODO: Add docstring"""
        return self._type

    @property
    def unit(self) -> Optional[str]:
        """TODO: Add docstring"""
        return self._unit

    def to_dict(self) -> Dict[str, str]:
        """TODO: Add docstring"""
        kwargs = {"name": self._name}
        if self._type is not None:
            kwargs["type"] = self._type
        if self._unit is not None:
            kwargs["unit"] = self._unit
        return kwargs

    def to_json(self) -> str:
        """TODO: Add docstring"""
        return json.dumps(self.to_dict())


class CoordinateSystem(coordinates.CoordinateSystem):
    """A coordinate system for spatial data."""

    @classmethod
    def from_json(cls, data: str) -> Self:
        """TODO: Add docstring"""
        # TODO: Needs good, comprehensive error handling.
        raw = json.loads(data)
        return cls(tuple(Axis.from_dict(axis) for axis in raw))

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


class CoordinateTransform(coordinates.CoordinateTransform, metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def to_json(self) -> str:
        """TODO: Add docstring"""
        raise NotImplementedError()

    @abc.abstractmethod
    def to_dict(self) -> Dict[str, Any]:
        """TODO: Add docstring"""


class IdentityTransform(CoordinateTransform):
    """TODO: Add docstring"""

    def to_dict(self) -> Dict[str, Any]:
        """TODO: Add docstring"""
        return {"type": "identity"}

    def to_numpy(self) -> npt.NDArray[np.float64]:
        """TODO: Add docstring"""
        raise NotImplementedError()

    def to_json(self) -> str:
        """TODO: Add docstring"""
        return json.dumps(self.to_dict())


class ScaleTransform(CoordinateTransform):
    """TODO: Add docstring"""

    def __init__(self, scale: npt.ArrayLike):
        self._scale: npt.NDArray[np.float64] = np.array(scale, dtype=np.float64)
        self._nvalues = len(self._scale)

    @property
    def scale(self) -> npt.NDArray[np.float64]:
        """TODO: Add docstring"""
        return self._scale

    @property
    def shape(self) -> Tuple[int, int]:
        """TODO: Add docstring"""
        return (self._nvalues, self._nvalues)

    def to_dict(self) -> Dict[str, Any]:
        """TODO: Add docstring"""
        return {"type": "scale", "scale": self._scale.tolist()}

    def to_json(self) -> str:
        """TODO: Add docstring"""
        return json.dumps(self.to_dict())

    def to_numpy(self) -> npt.NDArray[np.float64]:
        """TODO: Add docstring"""
        raise NotImplementedError()


class TranslateTransform(CoordinateTransform):
    """TODO: Add docstring"""

    def __init__(self, translate: npt.ArrayLike):
        self._translate: npt.NDArray[np.float64] = np.array(translate, dtype=np.float64)
        self._nvalues = len(self._translate)

    def shape(self) -> Tuple[int, int]:
        """TODO: Add docstring"""
        return (self._nvalues, self._nvalues)

    def to_dict(self) -> Dict[str, Any]:
        """TODO: Add docstring"""
        return {"type": "translate", "translate": self._translate.tolist()}

    def to_json(self) -> str:
        """TODO: Add docstring"""
        return json.dumps(self.to_dict())

    def to_numpy(self) -> npt.NDArray[np.float64]:
        """TODO: Add docstring"""
        raise NotImplementedError()


class CompositeTransform(coordinates.CoordinateTransform):
    """TODO: Add docstring"""

    @classmethod
    def from_json(cls, data: str) -> Self:
        """TODO: Add docstring"""
        raw = json.loads(data)
        transforms: List[CoordinateTransform] = []
        for transform in raw:
            transform_type = transform["type"]
            if transform_type == "identity":
                transforms.append(IdentityTransform())
            elif transform_type == "scale":
                transforms.append(ScaleTransform(transform["scale"]))
            elif transform_type == "translate":
                transforms.append(TranslateTransform(transform["translate"]))
            elif transform_type == "affine":
                raise NotImplementedError()
            else:
                raise ValueError()
        return cls(transforms)

    def __init__(self, transforms: Sequence[CoordinateTransform]):
        self._transforms = tuple(transforms)

    def __len__(self) -> int:
        return len(self._transforms)

    def __getitem__(self, key: int) -> CoordinateTransform:
        return self._transforms[key]

    def to_json(self) -> str:
        """TODO: Add docstring"""
        tmp = [transform.to_dict() for transform in self._transforms]
        return json.dumps(tmp)

    def to_numpy(self) -> npt.NDArray[np.float64]:
        """TODO: Add docstring"""
        raise NotImplementedError()
