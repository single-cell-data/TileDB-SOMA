import tiledb, numpy as np
from tiledb import *

# Reference/label array

s1 = ArraySchema(
    domain=Domain(
        *[
            Dim(name="z_name", tile=None, dtype="ascii"),
        ]
    ),
    attrs=[
        Attr(name="z_idx", dtype="uint32", var=False, nullable=False),
    ],
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
    sparse=True,
    allows_duplicates=True,
    coords_filters=FilterList([ZstdFilter(level=-1)]),
)

tiledb.Array.create("ref", s1)

ref_data = [
    "the",
    "quick",
    "brown",
    "fox",
    "jumped",
    "over",
    "the",
    "lazy",
    "spotted",
    "dog",
]
with tiledb.open("ref", "w") as A:
    A[ref_data] = np.arange(10)

# Target array

s2 = ArraySchema(
    domain=Domain(
        *[
            Dim(name="z_idx", domain=(0, 99), tile=10, dtype="uint32"),
        ]
    ),
    attrs=[
        Attr(name="data", dtype="int32", var=False, nullable=False),
    ],
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
    sparse=True,
    allows_duplicates=False,
    coords_filters=FilterList([ZstdFilter(level=-1)]),
)

tiledb.Array.create("tgt", s2)

with tiledb.open("tgt", "w") as A:
    A[np.arange(10)] = np.arange(30, 40)
