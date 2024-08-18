from __future__ import annotations

import pyarrow as pa
import pytest

import tiledbsoma

from tests._util import maybe_raises

# ================================================================
# SHORT LIST:
#
# k snda creation: with shape and maxshape
#
# k snda accessor: shape
#   o try fallback on old data (check in to repo)
# k snda accessor: maxshape
# k snda accessor: used_shape
#   o deprecation notice ...
# k snda accessor: non_empty_domain
#
# k snda bounds-checking on reads
# k snda bounds-checking on writes
#
# k snda mutator: resize
#   o raise NotImplementedError for old arrays
#
# TODO: non-2D SNDA cases
#
# ----------------------------------------------------------------
# * sdf creation: with domain and ... ? also shape or maxshape? implicit?
#
# * sdf accessor: shape
#   o try fallback on old data (check in to repo)
# k sdf accessor: maxshape
# k sdf accessor: used_shape -- does not exist anyway
#   k deprecation notice -- b/c it does not exist anyway
# k sdf accessor: non_empty_domain
# k sdf accessor: domain
#
# k sdf bounds-checking on reads
# k sdf bounds-checking on writes
#
# * sdf mutator: resize
#   o raise NotImplementedError for old arrays
# * sdf mutator: tiledbsoma_upgrade_shape
#   o no-op for new arrays -- ?
#
# k all: partials w/ extra dims
#
# ----------------------------------------------------------------
# * both: raise IndexError rather than SOMAError for OOB accesses??
#
# * both mutator: tiledbsoma_upgrade_shape
#   o array.schema.version to see if needed
#   o no-op for new arrays -- ?
#   o use core storage-version-update logic
#   o fail if outside domain
#
# ----------------------------------------------------------------
# * experiment mutator: tiledbsoma.io.resize
#   o do-it-all w/ new nobs/nvar -- ?
# ================================================================


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
        [(100,), None],
        [(100,), None],
        [(100, 200), None],
        [(100, 200, 300), None],
        [(100, 200), None],
        [(100, 200), None],
        [(100, 200), None],
        [(100, 200), None],
        [(100, 200), None],
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

        # TODO: need a saved-off array in UT-data land

        # We expect XXX to be set to a big signed int32. (There are details on the exact value of
        # that number, involving R compatibility, and leaving room for a single tile capacity, etc
        # ...  we could check for some magic value but it suffices to check that it's over 2
        # billion.)
        for e in snda.maxshape:
            assert e > 2_000_000_000

        # TODO: used_shape
        #   o as-is
        #   o deprecation notice ...

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
        for e in snda.maxshape:
            assert e > 2_000_000_000
        assert snda.non_empty_domain() == coords

    # Test reads out of bounds
    with tiledbsoma.SparseNDArray.open(uri) as snda:
        with pytest.raises(ValueError):
            coords = tuple([arg_shape[i] + 10 for i in range(ndim)])
            snda.read(coords).tables().concat()

    # Test writes out of bounds
    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        with pytest.raises(tiledbsoma.SOMAError):
            dikt = {"soma_data": [30]}
            for i in range(ndim):
                dikt[dim_names[i]] = [arg_shape[i] + 20]
            table = pa.Table.from_pydict(dikt)
            snda.write(table)

    with tiledbsoma.SparseNDArray.open(uri) as snda:
        assert snda.shape == arg_shape

    # Test resize down
    new_shape = tuple([arg_shape[i] - 50 for i in range(ndim)])
    # TODO: why union with tiledb.cc.TileDBError -- needed in sandbox
    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        with pytest.raises(ValueError):
            snda.resize(new_shape)

    with tiledbsoma.SparseNDArray.open(uri) as snda:
        assert snda.shape == arg_shape

    # Test resize too big
    new_shape = tuple([4_000_000_000 for i in range(ndim)])
    with pytest.raises(ValueError):
        with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
            snda.resize(new_shape)
    with tiledbsoma.SparseNDArray.open(uri) as snda:
        assert snda.shape == arg_shape

    # Test reasonable resize
    new_shape = tuple([arg_shape[i] + 50 for i in range(ndim)])
    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        snda.resize(new_shape)

        dikt = {}
        for i in range(ndim):
            dikt[dim_names[i]] = [arg_shape[i] + 20]
        dikt["soma_data"] = pa.array([34.5], type=element_dtype)
        table = pa.Table.from_pydict(dikt)

        # Re-test writes out of old bounds, within new bounds
        with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
            # Implicitly checking there's no raise
            snda.write(table)

        # Re-test reads out of old bounds, within new bounds
        with tiledbsoma.SparseNDArray.open(uri) as snda:
            assert snda.shape == new_shape

            coords = tuple([(arg_shape[i] + 20,) for i in range(ndim)])
            # Implicitly checking there's no raise
            readback = snda.read(coords).tables().concat()
            assert readback == table


# Pending 2.26 timeframe for dense support
# TODO: mark these with a linked GitHub tracking issue
def test_dense_nd_array_basics(tmp_path):
    uri = tmp_path.as_posix()
    shape = (100, 200)
    tiledbsoma.DenseNDArray.create(uri, type=pa.float64(), shape=shape)

    with tiledbsoma.DenseNDArray.open(uri) as dnda:
        assert dnda.shape == (100, 200)

        with pytest.raises(NotImplementedError):
            assert dnda.maxshape == (100, 200)

        with pytest.raises(NotImplementedError):
            dnda.resize((300, 400))

        assert dnda.non_empty_domain() == ((0, 0), (0, 0))

        assert dnda.shape == (100, 200)


@pytest.mark.parametrize(
    "domain0",
    [
        None,
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
    ],
)
def test_dataframe_basics(tmp_path, domain0, index_column_names):
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

    if domain0 is None:
        shape0 = None
    else:
        shape0 = domain0[1] + 1
    domain = [None] * len(index_column_names)
    domain[0] = domain0
    domain = tuple(domain)

    with tiledbsoma.DataFrame.create(
        uri,
        schema=schema,
        index_column_names=index_column_names,
        domain=domain,
    ) as sdf:
        if shape0 is not None and len(data) > shape0:
            with pytest.raises(tiledbsoma.SOMAError):
                sdf.write(data)
        else:
            sdf.write(data)

    with tiledbsoma.DataFrame.open(uri) as sdf:
        assert len(sdf.shape) == 1
        if shape0 is not None:
            assert sdf.shape[0] == shape0
        assert len(sdf.maxshape) == 1
        assert sdf.maxshape[0] > 2_000_000_000  # XXX COMMENT
        assert len(sdf.non_empty_domain()) == len(index_column_names)

    # XXX guard against ... not just the type-checker ...
    if domain0 is not None:
        new_size = 10000
        with tiledbsoma.DataFrame.open(uri, "r") as sdf:
            # Must be open for write.
            # XXX TO DO fix this
            # with pytest.raises(tiledbsoma.SOMAError):
            with pytest.raises(ValueError):
                sdf.resize([new_size])
        with tiledbsoma.DataFrame.open(uri, "w") as sdf:
            sdf.resize([new_size])

        with tiledbsoma.DataFrame.open(uri, "w") as sdf:
            sdf.write(data)

        with tiledbsoma.DataFrame.open(uri) as sdf:
            assert len(sdf.shape) == 1
            assert sdf.shape[0] == new_size
            assert len(sdf.maxshape) == 1
            assert sdf.maxshape[0] > 2_000_000_000  # XXX COMMENT
            assert len(sdf.non_empty_domain()) == len(index_column_names)

    # XXX MORE
    # XXX have here too a saved-off old array
    # XXX new OOB test
