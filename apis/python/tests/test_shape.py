from __future__ import annotations

import pyarrow as pa
import pytest

import tiledbsoma

from tests._util import maybe_raises


@pytest.mark.parametrize(
    "element_dtype",
    [
        pa.float64(),
        pa.float32(),
        pa.int64(),
        pa.uint16(),
    ],
)
@pytest.mark.parametrize(
    "shape_exc",
    [
        # Note: non-None exceptions are coming on https://github.com/single-cell-data/TileDB-SOMA/issues/2407
        [(100,), None],
        [(100, 200), None],
        [(100, 200, 300), None],
    ],
)
def test_sparse_nd_array_basics(
    tmp_path,
    element_dtype,
    shape_exc,
):
    uri = tmp_path.as_posix()
    arg_shape, arg_create_exc = shape_exc
    ndim = len(arg_shape)

    # Create the array
    with maybe_raises(arg_create_exc):
        snda = tiledbsoma.SparseNDArray.create(
            uri,
            type=element_dtype,
            shape=arg_shape,
        )
    if arg_create_exc is not None:
        return

    assert tiledbsoma.SparseNDArray.exists(uri)

    # Test the various accessors
    with tiledbsoma.SparseNDArray.open(uri) as snda:

        assert snda.shape == arg_shape

        # More to come on https://github.com/single-cell-data/TileDB-SOMA/issues/2407
        assert not snda.has_upgraded_shape

        # Before current-domain support: shape is maxshape.
        #
        # With current-domain support: We expect the maxshape to be set to a big
        # signed int32. (There are details on the exact value of that number,
        # involving R compatibility, and leaving room for a single tile
        # capacity, etc ...  we could check for some magic value but it suffices
        # to check that it's over 2 billion.)
        assert snda.shape == snda.maxshape
        # for e in snda.maxshape:
        #    assert e > 2_000_000_000

        # No data have been written for this test case
        assert snda.non_empty_domain() == tuple([(0, 0)] * ndim)

    # soma_dim_0: (0,1)
    # soma_dim_1: (2,3)
    # soma_dim_2: (4,5)
    coords = []
    dim_names = []
    for i in range(ndim):
        dim_names.append(f"soma_dim_{i}")
        coords.append((2 * i, 2 * i + 1))
    coords = tuple(coords)

    # Write some data
    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        dikt = {"soma_data": [4, 5]}
        for i in range(ndim):
            dikt[dim_names[i]] = coords[i]
        table = pa.Table.from_pydict(dikt)
        snda.write(table)

    # Test the various accessors
    with tiledbsoma.SparseNDArray.open(uri) as snda:
        assert snda.shape == arg_shape
        # This will change with current-domain support
        assert snda.shape == snda.maxshape
        # for e in snda.maxshape:
        #    assert e > 2_000_000_000
        assert snda.non_empty_domain() == coords

    # Test reads out of bounds
    with tiledbsoma.SparseNDArray.open(uri) as snda:
        # https://github.com/single-cell-data/TileDB-SOMA/issues/2407
        with pytest.raises(tiledbsoma.SOMAError):
            coords = tuple(arg_shape[i] + 10 for i in range(ndim))
            snda.read(coords).tables().concat()

    # Test writes out of bounds
    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        with pytest.raises(tiledbsoma.SOMAError):
            dikt = {"soma_data": [30]}
            dikt = {name: [shape + 20] for name, shape in zip(dim_names, arg_shape)}
            table = pa.Table.from_pydict(dikt)
            snda.write(table)

    with tiledbsoma.SparseNDArray.open(uri) as snda:
        assert snda.shape == arg_shape
        assert snda.shape == snda.maxshape


## Pending 2.27 timeframe for dense support for current domain, including resize
## TODO: mark these with a linked GitHub tracking issue
def test_dense_nd_array_basics(tmp_path):
    uri = tmp_path.as_posix()
    shape = (100, 200)
    tiledbsoma.DenseNDArray.create(uri, type=pa.float64(), shape=shape)

    with tiledbsoma.DenseNDArray.open(uri) as dnda:
        assert dnda.shape == (100, 200)
        assert dnda.maxshape == (100, 200)

        assert dnda.non_empty_domain() == ((0, 0), (0, 0))


@pytest.mark.parametrize(
    "soma_joinid_domain",
    [
        # TODO: https://github.com/single-cell-data/TileDB-SOMA/issues/2407
        # None,
        (0, 1),
        (0, 3),
        (0, 100),
    ],
)
@pytest.mark.parametrize(
    "index_column_names",
    [
        ["soma_joinid"],
        ["soma_joinid", "myint"],
        ["soma_joinid", "mystring"],
        ["mystring", "myint"],
    ],
)
def test_dataframe_basics(tmp_path, soma_joinid_domain, index_column_names):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("mystring", pa.string()),
            ("myint", pa.int16()),
            ("myfloat", pa.float32()),
        ]
    )

    data = pa.Table.from_pydict(
        {
            "soma_joinid": [0, 1, 2, 3],
            "mystring": ["a", "b", "a", "b"],
            "myint": [20, 30, 40, 50],
            "myfloat": [1.0, 2.5, 4.0, 5.5],
        }
    )

    domain_slots = {
        "soma_joinid": soma_joinid_domain,
        "mystring": None,
        "myint": (-1000, 1000),
        "myfloat": (-999.5, 999.5),
    }

    domain = tuple([domain_slots[name] for name in index_column_names])

    soma_joinid_coords = data["soma_joinid"]
    oob_write = any(
        e.as_py() < soma_joinid_domain[0] or e.as_py() > soma_joinid_domain[1]
        for e in soma_joinid_coords
    )
    oob_write = oob_write and "soma_joinid" in index_column_names

    with tiledbsoma.DataFrame.create(
        uri,
        schema=schema,
        index_column_names=index_column_names,
        domain=domain,
    ) as sdf:
        if oob_write:
            with pytest.raises(tiledbsoma.SOMAError):
                sdf.write(data)
        else:
            sdf.write(data)

    with tiledbsoma.DataFrame.open(uri) as sdf:
        has_sjid_dim = "soma_joinid" in index_column_names
        if has_sjid_dim:
            assert sdf._maybe_soma_joinid_shape == 1 + soma_joinid_domain[1]
            assert sdf._maybe_soma_joinid_maxshape == 1 + soma_joinid_domain[1]
        else:
            assert sdf._maybe_soma_joinid_shape is None
            assert sdf._maybe_soma_joinid_maxshape is None

        assert len(sdf.non_empty_domain()) == len(index_column_names)