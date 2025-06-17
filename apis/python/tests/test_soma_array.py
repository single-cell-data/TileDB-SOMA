#!/usr/bin/env python

import os

import pyarrow as pa
import pytest

import tiledbsoma.pytiledbsoma as clib

VERBOSE = False

TEST_DIR = os.path.dirname(__file__)
SOMA_URI = f"{TEST_DIR}/../../../data/soco/pbmc3k_processed"

if VERBOSE:
    clib.config_logging("debug")


def test_soma_array_basic_getters():
    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    column_names = [
        "soma_joinid",
        "obs_id",
        "n_genes",
        "percent_mito",
        "n_counts",
        "louvain",
        "is_b_cell",
    ]
    result_order = clib.ResultOrder.colmajor

    sr = clib.SOMAArray(uri=uri)
    assert sr.uri == uri
    assert sr.nnz() == 2638

    mq = clib.ManagedQuery(sr, sr.context(), name)
    mq.set_layout(result_order)
    mq.select_columns(column_names)
    assert mq.result_order == result_order
    assert mq.column_names == column_names


def test_soma_array_obs():
    """Read all values from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == 2638


def test_soma_array_var():
    """Read all values from var array into an arrow table."""

    name = "var"
    uri = os.path.join(SOMA_URI, "ms/RNA", name)
    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == 1838


def test_soma_array_var_x_data():
    """Read all values from x/data array into an arrow table."""

    name = "X/data"
    uri = os.path.join(SOMA_URI, "ms/RNA", name)
    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())

    # iterate read batches until all results have been processed
    total_num_rows = 0
    for arrow_table in iter(mq.next, StopIteration):
        total_num_rows += arrow_table.num_rows

    assert total_num_rows == 4848644


def test_soma_array_dim_points():
    """Read scalar dimension slice from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)

    obs_id_points = list(range(0, 100, 2))

    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_dim_points_int64("soma_joinid", obs_id_points)
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == len(obs_id_points)


def test_soma_array_empty_dim_points():
    """Read scalar dimension slice from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)

    obs_id_points = []

    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_dim_points_int64("soma_joinid", obs_id_points)
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == len(obs_id_points)


def test_soma_array_dim_points_arrow_array():
    """Read scalar dimension slice from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)

    obs_id_points = pa.array([0, 2, 4, 6, 8])

    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_dim_points_arrow("soma_joinid", obs_id_points)
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == len(obs_id_points)


def test_soma_array_dim_ranges():
    """Read range dimension slice from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)

    obs_id_ranges = [
        [1000, 1004],
        [2000, 2004],
    ]

    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_dim_ranges_int64("soma_joinid", obs_id_ranges)
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == 10


def test_soma_array_dim_mixed():
    """Read scalar and range dimension slice from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)

    obs_id_points = list(range(0, 100, 2))

    obs_id_ranges = [
        [1000, 1004],
        [2000, 2004],
    ]

    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_dim_points_int64("soma_joinid", obs_id_points)
    mq.set_dim_ranges_int64("soma_joinid", obs_id_ranges)
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == 60


def test_soma_array_obs_slice_x():
    """Read X/data sliced by obs."""

    # read obs
    # ---------------------------------------------------------------1
    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)

    obs_id_points = list(range(0, 100, 2))

    obs_id_ranges = [
        [1000, 1004],
        [2000, 2004],
    ]

    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_dim_points_int64("soma_joinid", obs_id_points)
    mq.set_dim_ranges_int64("soma_joinid", obs_id_ranges)
    obs = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert obs.num_rows == 60

    # read X/data
    # ---------------------------------------------------------------1
    name = "X/data"
    uri = os.path.join(SOMA_URI, "ms/RNA", name)
    sr = clib.SOMAArray(uri)

    # slice X/data read with obs.soma_joinid column
    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_dim_points_arrow("soma_dim_0", obs.column("soma_joinid"))

    # iterate read batches until all results have been processed
    total_num_rows = 0
    for x_data in iter(mq.next, StopIteration):
        total_num_rows += x_data.num_rows

    assert total_num_rows == 110280


def test_soma_array_column_names():
    """Read specified column names of obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())
    mq.select_columns(["soma_joinid", "louvain"])
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_columns == 2


def test_soma_array_reset():
    """Submit a query with a SOMAArray object, reset the SOMAArray, and submit another query."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = clib.SOMAArray(uri)

    mq = clib.ManagedQuery(sr, sr.context())
    mq.select_columns(["soma_joinid", "louvain"])
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_columns == 2
    assert arrow_table.num_rows == 2638

    # reset and submit new query with open array
    # ---------------------------------------------------------------
    obs_id_points = pa.array([0, 2, 4, 6, 8])
    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_dim_points_arrow("soma_joinid", obs_id_points)
    arrow_table = mq.next()

    # test that all results are present in the arrow table (no incomplete queries)
    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_columns == 7
    assert arrow_table.num_rows == 5


if __name__ == "__main__":
    test_soma_array_obs_slice_x()
