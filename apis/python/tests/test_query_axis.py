from typing import Any

import numpy as np
import pytest
from pytest import mark

import tiledbsoma
from tiledbsoma._core_options import SparseDFCoord


@mark.parametrize(
    ["coords", "want"],
    [
        ((), ()),
        ((slice(1, 10),), (slice(1, 10),)),
        ([0, 1, 2], (0, 1, 2)),
        ([1, 1.5, 2], (1, 1.5, 2)),
        ((slice(None), [0, 88, 1001]), (slice(None), (0, 88, 1001))),
        ((slice(2.5, 3.5),), (slice(2.5, 3.5),)),
        (
            (slice(np.datetime64(946684802, "s"), np.datetime64(946684803, "s")),),
            (slice(np.datetime64(946684802, "s"), np.datetime64(946684803, "s")),),
        ),
        (("string-coord", [b"lo", b"hi"]), ("string-coord", (b"lo", b"hi"))),
        ((slice(4, 5), True, None), (slice(4, 5), True, None)),
    ],
)
def test_canonicalization(coords: Any, want: tuple[SparseDFCoord, ...]) -> None:
    axq = tiledbsoma.AxisQuery(coords=coords)
    assert want == axq.coords


def test_canonicalization_nparray() -> None:
    axq = tiledbsoma.AxisQuery(coords=(1, np.array([1, 2, 3])))

    one, arr = axq.coords
    assert one == 1
    assert (np.array([1, 2, 3]) == arr).all()


@mark.parametrize(
    ["coords"],
    [
        ("forbid bare strings",),
        (b"forbid bare byteses",),
    ],
)
def test_canonicalization_bad(coords) -> None:
    with pytest.raises(TypeError):
        tiledbsoma.AxisQuery(coords=coords)
