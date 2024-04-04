import numpy as np

import tiledb

# Reference/label array

s1 = tiledb.ArraySchema(
    domain=tiledb.Domain(
        tiledb.Dim(name="z_name", tile=None, dtype="ascii"),
    ),
    attrs=[
        tiledb.Attr(name="z_idx", dtype="uint32", var=False, nullable=False),
    ],
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
    sparse=True,
    allows_duplicates=True,
    coords_filters=tiledb.FilterList([tiledb.ZstdFilter(level=-1)]),
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

s2 = tiledb.ArraySchema(
    domain=tiledb.Domain(
        tiledb.Dim(name="z_idx", domain=(0, 99), tile=10, dtype="uint32"),
    ),
    attrs=[
        tiledb.Attr(name="data", dtype="int32", var=False, nullable=False),
    ],
    cell_order="row-major",
    tile_order="row-major",
    capacity=10000,
    sparse=True,
    allows_duplicates=False,
    coords_filters=tiledb.FilterList([tiledb.ZstdFilter(level=-1)]),
)

tiledb.Array.create("tgt", s2)

with tiledb.open("tgt", "w") as A:
    A[np.arange(10)] = np.arange(30, 40)
