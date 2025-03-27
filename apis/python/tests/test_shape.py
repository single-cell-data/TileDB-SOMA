from __future__ import annotations

import io
import tarfile

import pyarrow as pa
import pytest

import tiledbsoma
import tiledbsoma.io

from ._util import TESTDATA


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
    "arg_shape",
    [
        (100,),
        (100, 200),
        (100, 200, 300),
    ],
)
def test_sparse_nd_array_basics(
    tmp_path,
    element_dtype,
    arg_shape,
):
    uri = tmp_path.as_posix()
    ndim = len(arg_shape)

    # Create the array
    tiledbsoma.SparseNDArray.create(
        uri,
        type=element_dtype,
        shape=arg_shape,
    )

    assert tiledbsoma.SparseNDArray.exists(uri)

    # Test the various accessors
    with tiledbsoma.SparseNDArray.open(uri) as snda:

        assert snda.shape == arg_shape

        assert snda.tiledbsoma_has_upgraded_shape

        # Before current-domain support: shape is maxshape.
        #
        # With current-domain support: We expect the maxshape to be set to a big
        # signed int32. (There are details on the exact value of that number,
        # involving R compatibility, and leaving room for a single tile
        # capacity, etc ...  we could check for some magic value but it suffices
        # to check that it's over 2 billion.)

        for e in snda.maxshape:
            assert e > 2_000_000_000

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
        for e in snda.maxshape:
            assert e > 2_000_000_000
        assert snda.non_empty_domain() == coords

    # Test reads out of bounds
    with tiledbsoma.SparseNDArray.open(uri) as snda:
        coords = tuple(arg_shape[i] + 10 for i in range(ndim))
        with pytest.raises(tiledbsoma.SOMAError):
            snda.read(coords).tables().concat()

    # Test writes out of bounds
    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        with pytest.raises(tiledbsoma.SOMAError):
            dikt = {name: [shape + 20] for name, shape in zip(dim_names, arg_shape)}
            dikt["soma_data"] = [30]
            table = pa.Table.from_pydict(dikt)
            snda.write(table)

    with tiledbsoma.SparseNDArray.open(uri) as snda:
        assert snda.shape == arg_shape

    with tiledbsoma.SparseNDArray.open(uri) as snda:
        ok, msg = snda.tiledbsoma_upgrade_shape(arg_shape, check_only=True)
        assert not ok
        assert (
            msg
            == "tiledbsoma_can_upgrade_shape: array already has a shape: please use resize"
        )

    # Test resize down
    new_shape = tuple([arg_shape[i] - 50 for i in range(ndim)])
    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        (ok, msg) = snda.resize(new_shape, check_only=True)
        assert not ok
        assert (
            msg
            == "[can_resize] index-column name 'soma_dim_0': new upper 49 < old upper 99 (downsize is unsupported)"
        )
        # TODO: check draft spec
        # with pytest.raises(ValueError):
        with pytest.raises(tiledbsoma.SOMAError):
            snda.resize(new_shape)

    with tiledbsoma.SparseNDArray.open(uri) as snda:
        assert snda.shape == arg_shape

    # Test writes out of bounds
    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        with pytest.raises(tiledbsoma.SOMAError):
            dikt = {name: [shape + 20] for name, shape in zip(dim_names, arg_shape)}
            dikt["soma_data"] = [30]
            table = pa.Table.from_pydict(dikt)
            snda.write(table)

    # Test resize
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

    with tiledbsoma.SparseNDArray.open(uri) as snda:
        assert snda.shape == new_shape

        (ok, msg) = snda.resize(new_shape, check_only=True)
        assert ok
        assert msg == ""

        too_small = tuple(e - 1 for e in new_shape)
        (ok, msg) = snda.resize(too_small, check_only=True)
        assert not ok
        assert (
            msg
            == "[can_resize] index-column name 'soma_dim_0': new upper 148 < old upper 149 (downsize is unsupported)"
        )

    with tiledbsoma.SparseNDArray.open(uri, "w") as snda:
        (ok, msg) = snda.resize(new_shape, check_only=True)


def test_dense_nd_array_basics(tmp_path):
    uri = tmp_path.as_posix()
    shape = (100, 200)
    tiledbsoma.DenseNDArray.create(uri, type=pa.float64(), shape=shape)

    with tiledbsoma.DenseNDArray.open(uri) as dnda:
        assert dnda.shape == (100, 200)
        assert len(dnda.maxshape)
        assert dnda.maxshape[0] > 2**62
        assert dnda.maxshape[1] > 2**62

        assert dnda.non_empty_domain() == ((0, 0), (0, 0))

    with tiledbsoma.DenseNDArray.open(uri, "w") as dnda:
        dnda.resize((300, 400))

    with tiledbsoma.DenseNDArray.open(uri) as dnda:
        assert dnda.non_empty_domain() == ((0, 0), (0, 0))
        assert dnda.shape == (300, 400)

    with tiledbsoma.DenseNDArray.open(uri) as dnda:
        ok, msg = dnda.tiledbsoma_upgrade_shape((600, 700), check_only=True)
        assert not ok
        assert (
            msg
            == "tiledbsoma_can_upgrade_shape: array already has a shape: please use resize"
        )


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

    data_dict = {
        "soma_joinid": [0, 1, 2, 3],
        "mystring": ["a", "b", "a", "b"],
        "myint": [20, 30, 40, 50],
        "myfloat": [1.0, 2.5, 4.0, 5.5],
    }

    data = pa.Table.from_pydict(data_dict)

    domain_slots = {
        "soma_joinid": soma_joinid_domain,
        "mystring": None,
        "myint": (-1000, 1000),
        "myfloat": (-999.5, 999.5),
    }

    has_soma_joinid_dim = "soma_joinid" in index_column_names

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
        else:
            assert sdf._maybe_soma_joinid_shape is None

        assert len(sdf.non_empty_domain()) == len(index_column_names)

        # This may be None if soma_joinid is not an index column
        shape_at_create = sdf._maybe_soma_joinid_shape

    # Test resize down
    new_shape = 0
    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        ok, msg = sdf.tiledbsoma_resize_soma_joinid_shape(new_shape, check_only=True)
        if has_soma_joinid_dim:
            # TODO: check draft spec
            # with pytest.raises(ValueError):
            assert not ok
            assert (
                "tiledbsoma_resize_soma_joinid_shape: new soma_joinid shape 0 < existing shape"
                in msg
            )
            with pytest.raises(tiledbsoma.SOMAError):
                sdf.tiledbsoma_resize_soma_joinid_shape(new_shape)
        else:
            assert ok
            assert msg == ""
            sdf.tiledbsoma_resize_soma_joinid_shape(new_shape)

    with tiledbsoma.DataFrame.open(uri) as sdf:
        assert sdf._maybe_soma_joinid_shape == shape_at_create

    # Test writes out of bounds, before resize
    offset = shape_at_create if has_soma_joinid_dim else 100
    data_dict["soma_joinid"] = [e + offset for e in data_dict["soma_joinid"]]
    data = pa.Table.from_pydict(data_dict)

    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        if has_soma_joinid_dim:
            with pytest.raises(tiledbsoma.SOMAError):
                sdf.write(data)
        else:
            sdf.write(data)

    # Test resize
    new_shape = 0 if shape_at_create is None else shape_at_create + 100
    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        sdf.tiledbsoma_resize_soma_joinid_shape(new_shape)

    # Test writes out of old bounds, within new bounds, after resize
    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        sdf.write(data)


def test_domain_mods(tmp_path):
    uri = tmp_path.as_posix()

    schema = pa.schema(
        [
            ("soma_joinid", pa.int64()),
            ("mystring", pa.string()),
            ("myint", pa.int16()),
            ("myfloat", pa.float32()),
            ("mybool", pa.bool_()),  # not supported as an index type
        ]
    )
    index_column_names = ["soma_joinid", "mystring", "myint", "myfloat"]

    domain_for_create = [
        [0, 3],
        None,
        [20, 50],
        [0.0, 6.0],
    ]

    data_dict = {
        "soma_joinid": [0, 1, 2, 3],
        "mystring": ["a", "b", "a", "b"],
        "myint": [20, 30, 40, 50],
        "myfloat": [1.0, 2.5, 4.0, 5.5],
        "mybool": [True, False, True, True],
    }

    data = pa.Table.from_pydict(data_dict)

    with tiledbsoma.DataFrame.create(
        uri,
        schema=schema,
        index_column_names=index_column_names,
        domain=domain_for_create,
    ) as sdf:
        sdf.write(data)

    # Check "expand" to same
    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        newdomain = [[0, 3], None, [20, 50], [0.0, 6.0]]
        ok, msg = sdf.change_domain(newdomain, check_only=True)
        assert ok
        assert msg == ""

    # Shrink
    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        newdomain = [[0, 2], None, [20, 50], [0.0, 6.0]]
        ok, msg = sdf.change_domain(newdomain, check_only=True)
        assert not ok
        assert "downsize is unsupported" in msg

    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        newdomain = [[0, 3], None, [20, 40], [0.0, 6.0]]
        ok, msg = sdf.change_domain(newdomain, check_only=True)
        assert not ok
        assert "downsize is unsupported" in msg

    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        newdomain = [[0, 3], None, [20, 50], [1.0, 6.0]]
        ok, msg = sdf.change_domain(newdomain, check_only=True)
        assert not ok
        assert "downsize is unsupported" in msg

    # String domain cannot be specified
    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        newdomain = [
            [0, 3],
            ["a", "z"],
            [20, 50],
            [0.0, 6.0],
        ]
        ok, msg = sdf.change_domain(newdomain, check_only=True)
        assert not ok
        assert "domain cannot be set for string index columns" in msg

    # All clear
    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        newdomain = [[0, 9], None, [0, 100], [-10.0, 10.0]]
        ok, msg = sdf.change_domain(newdomain, check_only=True)
        assert ok
        assert msg == ""
        sdf.change_domain(newdomain)

    # Check for success
    with tiledbsoma.DataFrame.open(uri, "r") as sdf:
        dom = sdf.domain
        assert dom[0] == (0, 9)
        assert dom[1] == ("", "")
        assert dom[2] == (0, 100)
        assert dom[3] == (-10.0, 10.0)


@pytest.mark.parametrize("has_shapes", [False, True])
def test_canned_experiments(tmp_path, has_shapes):
    uri = tmp_path.as_posix()

    if not has_shapes:
        tgz = TESTDATA / "pbmc-exp-without-shapes.tgz"
    else:
        tgz = TESTDATA / "pbmc-exp-with-shapes.tgz"

    with tarfile.open(tgz) as handle:
        handle.extractall(uri)

    def _assert_huge_domainish(d):
        assert len(d) == 1
        assert len(d[0]) == 2
        assert d[0][0] == 0
        # Exact number depends on tile extent, and is unimportant in any case
        assert d[0][1] > 2**62

    def _check_dataframe(
        sdf, has_shapes, expected_count, *, count_must_match: bool = True
    ):
        if count_must_match:
            # OK match case: 2000 populated rows and shape is 2000.
            # OK mismatch case: 2000 populated rows and a reshape to 3000 has been done.
            assert sdf.count == expected_count
        if not has_shapes:
            _assert_huge_domainish(sdf.domain)
        else:
            assert sdf.domain == ((0, expected_count - 1),)
        _assert_huge_domainish(sdf.maxdomain)
        assert sdf.tiledbsoma_has_upgraded_domain == has_shapes

    def _assert_huge_shape(d):
        assert len(d) == 2
        # Exact number depends on tile extent, and is unimportant in any case
        assert d[0] > 2**62
        assert d[1] > 2**62

    def _check_ndarray(ndarray, has_shapes, expected_shape):
        if not has_shapes:
            _assert_huge_shape(ndarray.shape)
        else:
            assert ndarray.shape == expected_shape
        _assert_huge_shape(ndarray.maxshape)

    with tiledbsoma.Experiment.open(uri) as exp:

        _check_dataframe(exp.obs, has_shapes, 2638)

        assert "raw" in exp.ms
        assert "data" in exp.ms["raw"].X

        _check_dataframe(exp.ms["raw"].var, has_shapes, 13714)
        _check_ndarray(exp.ms["raw"].X["data"], has_shapes, (2638, 13714))

        _check_dataframe(exp.ms["RNA"].var, has_shapes, 1838)
        _check_ndarray(exp.ms["RNA"].X["data"], has_shapes, (2638, 1838))
        _check_ndarray(exp.ms["RNA"].obsm["X_pca"], has_shapes, (2638, 50))
        _check_ndarray(exp.ms["RNA"].obsp["connectivities"], has_shapes, (2638, 2638))
        _check_ndarray(exp.ms["RNA"].varm["PCs"], has_shapes, (1838, 50))

    # Check tiledbsoma.io.show_experiment_shapes. This is mainly for interactive
    # use, so in a unit test, we want to check basics:
    # * Check that many lines were printed
    # * Make sure exceptions occurred
    # * Do some spot-checks on wording

    handle = io.StringIO()
    tiledbsoma.io.show_experiment_shapes(uri, output_handle=handle)
    handle.seek(0)
    lines = handle.readlines()
    handle.close()
    # Exact line count doesn't matter: make sure it's a lot.
    assert len(lines) > 50
    body = "\n".join(lines)
    assert "[SparseNDArray] ms/RNA/obsp/distances" in body
    assert "ms/RNA/obsm/X_draw_graph_fr" in body

    # Check upgrade_domain for dataframes
    with tiledbsoma.Experiment.open(uri, "w") as exp:

        ok, msg = exp.obs.tiledbsoma_upgrade_domain([[10, 4]], check_only=True)
        if has_shapes:
            assert not ok
            assert "dataframe already has a domain" in msg
        else:
            assert not ok
            assert "new lower 10 > new upper 4" in msg

        ok, msg = exp.obs.tiledbsoma_upgrade_domain([[0, 1]], check_only=True)
        if has_shapes:
            assert not ok
            assert "dataframe already has a domain" in msg
        else:
            assert ok
            assert msg == ""

        with pytest.raises(ValueError):
            exp.obs.tiledbsoma_upgrade_domain([[0, 1, 2]], check_only=True)

    # Check dry run of tiledbsoma.io.upgrade_experiment_shapes
    handle = io.StringIO()
    upgradeable = tiledbsoma.io.upgrade_experiment_shapes(
        uri, check_only=True, output_handle=handle
    )
    handle.seek(0)
    lines = handle.readlines()
    handle.close()
    # Exact line count doesn't matter: make sure it's a lot.
    assert len(lines) > 50
    body = "\n".join(lines)
    assert "[SparseNDArray] ms/RNA/obsp/distances" in body
    assert "ms/RNA/obsm/X_draw_graph_fr" in body
    assert upgradeable != has_shapes

    # Check dry run of tiledbsoma.io.resize_experiment -- plenty of room
    handle = io.StringIO()
    resizeable = tiledbsoma.io.resize_experiment(
        uri,
        nobs=100000,
        nvars={"RNA": 100000, "raw": 200000},
        check_only=True,
        output_handle=handle,
    )
    handle.seek(0)
    lines = handle.readlines()
    handle.close()
    # Exact line count doesn't matter: make sure it's a lot.
    assert len(lines) > 50
    body = "\n".join(lines)
    assert "[SparseNDArray] ms/RNA/obsp/distances" in body
    assert "ms/RNA/obsm/X_draw_graph_fr" in body
    assert resizeable == has_shapes

    # Check dry run of tiledbsoma.io.resize_experiment -- no change
    handle = io.StringIO()
    resizeable = tiledbsoma.io.resize_experiment(
        uri,
        nobs=2638,
        nvars={"RNA": 1838, "raw": 13714},
        check_only=True,
        output_handle=handle,
    )
    assert resizeable == has_shapes

    # Check dry run of tiledbsoma.io.resize_experiment -- downsize
    handle = io.StringIO()
    resizeable = tiledbsoma.io.resize_experiment(
        uri,
        nobs=2638,
        nvars={"RNA": 1838, "raw": 13713},
        check_only=True,
        output_handle=handle,
    )
    assert not resizeable
    handle.seek(0)
    lines = handle.readlines()
    handle.close()
    body = "\n".join(lines)
    if not has_shapes:
        assert (
            "Not OK: can_resize: array currently has no shape: please upgrade the array"
            in body
        )
    else:
        assert (
            "Not OK: [can_resize] index-column name 'soma_dim_1': new upper 13712 < old upper 13713 (downsize is unsupported)"
            in body
        )

    # Check real run of tiledbsoma.io.upgrade_experiment_shapes
    handle = io.StringIO()
    upgraded = tiledbsoma.io.upgrade_experiment_shapes(uri, output_handle=handle)
    handle.seek(0)
    lines = handle.readlines()
    handle.close()
    # Experiment-level upgrade is idempotent.
    assert upgraded

    # Check post-upgrade shapes
    with tiledbsoma.Experiment.open(uri) as exp:

        _check_dataframe(exp.obs, True, 2638)

        assert "raw" in exp.ms
        assert "data" in exp.ms["raw"].X

        _check_dataframe(exp.ms["raw"].var, True, 13714)
        _check_ndarray(exp.ms["raw"].X["data"], True, (2638, 13714))

        _check_dataframe(exp.ms["RNA"].var, True, 1838)
        _check_ndarray(exp.ms["RNA"].X["data"], True, (2638, 1838))
        _check_ndarray(exp.ms["RNA"].obsm["X_pca"], True, (2638, 50))
        _check_ndarray(exp.ms["RNA"].obsp["connectivities"], True, (2638, 2638))
        _check_ndarray(exp.ms["RNA"].varm["PCs"], True, (1838, 50))

    # Check real same-size resize
    handle = io.StringIO()
    resized = tiledbsoma.io.resize_experiment(
        uri,
        nobs=2638,
        nvars={"RNA": 1838, "raw": 13714},
        output_handle=handle,
    )
    assert resized

    # Check real down-size resize
    with pytest.raises(tiledbsoma.SOMAError):
        tiledbsoma.io.resize_experiment(
            uri,
            nobs=2637,
            nvars={"RNA": 1838, "raw": 13714},
            output_handle=handle,
        )

    # Check real up-size resize
    handle = io.StringIO()
    resized = tiledbsoma.io.resize_experiment(
        uri,
        nobs=2639,
        nvars={"RNA": 1839, "raw": 13720},
        output_handle=handle,
    )
    assert resized

    # Check new shapes
    with tiledbsoma.Experiment.open(uri) as exp:
        _check_dataframe(exp.obs, True, 2639, count_must_match=False)

        assert "raw" in exp.ms
        assert "data" in exp.ms["raw"].X

        _check_dataframe(exp.ms["raw"].var, True, 13720, count_must_match=False)
        _check_ndarray(exp.ms["raw"].X["data"], True, (2639, 13720))

        _check_dataframe(exp.ms["RNA"].var, True, 1839, count_must_match=False)
        _check_ndarray(exp.ms["RNA"].X["data"], True, (2639, 1839))
        _check_ndarray(exp.ms["RNA"].obsm["X_pca"], True, (2639, 50))
        _check_ndarray(exp.ms["RNA"].obsp["connectivities"], True, (2639, 2639))
        _check_ndarray(exp.ms["RNA"].varm["PCs"], True, (1839, 50))


def test_canned_nonstandard_dataframe_upgrade(tmp_path):
    uri = tmp_path.as_posix()

    tgz = TESTDATA / "nonstandard-dataframe-without-shapes.tgz"

    with tarfile.open(tgz) as handle:
        handle.extractall(uri)

    # ----------------------------------------------------------------
    # As of tiledbsoma 1.15 we no longer write dataframes/arrays without
    # core current domain (soma domain/shape) so it is crucial that in order
    # to test upgrade-shape, we test saved-off data written before 1.15.0.
    #
    # Here is a dataframe created with tiledbsoma 1.14.5:
    # ----------------------------------------------------------------
    # import tiledbsoma
    # import pyarrow as pa
    #
    # schema = pa.schema([
    #     ("soma_joinid", pa.int64()),
    #     ("mystring", pa.string()),
    #     ("myint", pa.int32()),
    #     ("myfloat", pa.float32()),
    # ])
    #
    # data = pa.Table.from_pydict({
    #     "soma_joinid": [10, 20],
    #     "mystring": ["hello", "world"],
    #     "myint": [1330, 1440],
    #     "myfloat": [4.5, 5.5],
    # })
    #
    # with tiledbsoma.DataFrame.create(
    #     uri="data-sdf-dom-multi-py",
    #     schema=schema,
    #     index_column_names=["soma_joinid", "myint", "mystring"],
    #     domain=None,
    # ) as sdf:
    #     sdf.write(data)
    # ----------------------------------------------------------------

    with tiledbsoma.DataFrame.open(uri) as sdf:
        assert not sdf.tiledbsoma_has_upgraded_domain
        assert sdf.non_empty_domain() == ((10, 20), (1330, 1440), ("hello", "world"))
        assert sdf.domain == ((0, 2147483646), (-2147483648, 2147481598), ("", ""))
        assert sdf.maxdomain == ((0, 2147483646), (-2147483648, 2147481598), ("", ""))

    with tiledbsoma.DataFrame.open(uri, "w") as sdf:
        ok, msg = sdf.tiledbsoma_upgrade_soma_joinid_shape(1, check_only=True)

        sdf.tiledbsoma_upgrade_soma_joinid_shape(100)

    with tiledbsoma.DataFrame.open(uri) as sdf:
        assert sdf.tiledbsoma_has_upgraded_domain
        assert sdf.non_empty_domain() == ((10, 20), (1330, 1440), ("hello", "world"))
        assert sdf.domain == ((0, 99), (-2147483648, 2147481598), ("", ""))
        assert sdf.maxdomain == ((0, 2147483646), (-2147483648, 2147481598), ("", ""))
