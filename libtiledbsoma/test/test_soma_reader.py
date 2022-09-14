import os
import pandas as pd
import pyarrow as pa
import tiledbsoma.libtiledbsoma as sc

VERBOSE = True

TEST_DIR = os.path.dirname(__file__)
SOMA_URI = f"{TEST_DIR}/../../test/soco/pbmc3k_processed"

if VERBOSE:
    sc.config_logging("debug")


def test_soma_reader_obs():
    """Read all values from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = sc.SOMAReader(uri)
    sr.submit()
    arrow_table = sr.read_next()

    # test that all results are present in the arrow table (no incomplete queries)
    assert sr.results_complete()
    assert arrow_table.num_rows == 2638


def test_soma_reader_var():
    """Read all values from var array into an arrow table."""

    name = "var"
    uri = os.path.join(SOMA_URI, "ms/mRNA", name)
    sr = sc.SOMAReader(uri)
    sr.submit()
    arrow_table = sr.read_next()

    # test that all results are present in the arrow table (no incomplete queries)
    assert sr.results_complete()
    assert arrow_table.num_rows == 1838


def test_soma_reader_var_x_data():
    """Read all values from x/data array into an arrow table."""

    name = "X/data"
    uri = os.path.join(SOMA_URI, "ms/mRNA", name)
    sr = sc.SOMAReader(uri)
    sr.submit()

    # iterate read batches until all results have been processed
    total_num_rows = 0
    while arrow_table := sr.read_next():
        total_num_rows += arrow_table.num_rows

    # test that all results are not present in the arrow table (incomplete queries)
    assert not sr.results_complete()
    assert total_num_rows == 4848644


def test_soma_reader_dim_points():
    """Read scalar dimension slice from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = sc.SOMAReader(uri)

    obs_id_points = list(range(0, 100, 2))

    sr.set_dim_points("soma_rowid", obs_id_points)

    sr.submit()
    arrow_table = sr.read_next()

    # test that all results are present in the arrow table (no incomplete queries)
    assert sr.results_complete()
    assert arrow_table.num_rows == len(obs_id_points)


def test_soma_reader_dim_ranges():
    """Read range dimension slice from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = sc.SOMAReader(uri)

    obs_id_ranges = [
        [1000, 1004],
        [2000, 2004],
    ]

    sr.set_dim_ranges("soma_rowid", obs_id_ranges)

    sr.submit()
    arrow_table = sr.read_next()

    # test that all results are present in the arrow table (no incomplete queries)
    assert sr.results_complete()
    assert arrow_table.num_rows == 10


def test_soma_reader_dim_mixed():
    """Read scalar and range dimension slice from obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = sc.SOMAReader(uri)

    obs_id_points = list(range(0, 100, 2))

    obs_id_ranges = [
        [1000, 1004],
        [2000, 2004],
    ]

    sr.set_dim_points("soma_rowid", obs_id_points)
    sr.set_dim_ranges("soma_rowid", obs_id_ranges)

    sr.submit()
    arrow_table = sr.read_next()

    # test that all results are present in the arrow table (no incomplete queries)
    assert sr.results_complete()
    assert arrow_table.num_rows == 60


def test_soma_reader_column_names():
    """Read specified column names of obs array into an arrow table."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = sc.SOMAReader(uri, column_names=["soma_rowid", "louvain"])

    sr.submit()
    arrow_table = sr.read_next()

    # test that all results are present in the arrow table (no incomplete queries)
    assert sr.results_complete()
    assert arrow_table.num_columns == 2


# TODO: implement value filter
def todo_soma_reader_value_filter():
    """Read obs array into an arrow table with value filter applied."""

    name = "obs"
    uri = os.path.join(SOMA_URI, name)
    sr = sc.SOMAReader(uri)

    sr.set_value_filter("foo > 1")

    sr.submit()
    arrow_table = sr.read_next()

    # test that all results are present in the arrow table (no incomplete queries)
    assert sr.results_complete()


if __name__ == "__main__":
    test_soma_reader_var()
