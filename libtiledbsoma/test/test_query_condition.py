#!/usr/bin/env python

import os

import pytest
import tiledb

import tiledbsoma.libtiledbsoma as clib
from tiledbsoma.exception import SOMAError
from tiledbsoma.query_condition import QueryCondition

VERBOSE = False

TEST_DIR = os.path.dirname(__file__)
SOMA_URI = f"{TEST_DIR}/../../test/soco/pbmc3k_processed"

if VERBOSE:
    clib.config_logging("debug")


def pandas_query(uri, condition):
    sr = clib.SOMAReader(uri)
    sr.submit()
    arrow_table = sr.read_next()
    assert sr.results_complete()

    return arrow_table.to_pandas().query(condition)


def soma_query(uri, condition):
    qc = QueryCondition(condition)
    schema = tiledb.open(uri).schema

    sr = clib.SOMAReader(uri, query_condition=qc, schema=schema)
    sr.submit()
    arrow_table = sr.read_next()
    assert sr.results_complete()

    return arrow_table


@pytest.mark.parametrize(
    "condition",
    [
        # types
        "n_genes > 500",  # int
        'louvain == "NK cells"',  # string
        "percent_mito > 0.02",  # float
        "is_b_cell == True",  # bool
        # compare_op
        "n_genes == 480",
        "n_genes != 480",
        "n_genes <= 500",
        "n_genes >= 500",
        "n_genes > 100",
        "n_genes < 100",
        "0 <= n_genes < 500",
        "louvain in ['CD4 T cells', 'B cells', 'NK cells']",
        "1 > n_genes",
        # unary ops
        "n_genes == +480",
        "n_genes >= -1",
        "n_genes > -(+(-1))",
        # boolean logic
        "percent_mito > 0.02 and n_genes > 700",  # and
        "percent_mito > 0.02 or n_genes > 700",  # or
        '(percent_mito > 0.02 and n_genes > 700) or (percent_mito < 0.015 and louvain == "B cells")',  # and or
        "(percent_mito > 0.02) & (n_genes > 700)",  # bit and
        "(percent_mito > 0.02) | (n_genes > 700)",  # bit or
    ],
)
def test_query_condition(condition):
    """Compare condition to result from Pandas - they should match"""
    uri = os.path.join(SOMA_URI, "obs")
    pandas = pandas_query(uri, condition)
    soma_arrow = soma_query(uri, condition)
    assert len(pandas.index) == soma_arrow.num_rows
    assert (
        (pandas.reset_index(drop=True) == soma_arrow.to_pandas().reset_index(drop=True))
        .all()
        .all()
    )


@pytest.mark.parametrize(
    "condition, pandas_equivalent_condition",
    [
        ["attr('n_genes') == 480", "n_genes == 480"],
        ["n_genes > val(100)", "n_genes > 100"],
        ["n_genes > val(100.0)", "n_genes > 100.0"],
        ["is_b_cell == val(True)", "is_b_cell == True"],
        ["louvain == val('B cells')", "louvain == 'B cells'"],
    ],
)
def test_query_condition_extensions(condition, pandas_equivalent_condition):
    """Test conditions which Pandas does not support directly - extensions"""
    uri = os.path.join(SOMA_URI, "obs")
    pandas = pandas_query(uri, pandas_equivalent_condition)
    soma_arrow = soma_query(uri, condition)
    assert len(pandas.index) == soma_arrow.num_rows
    assert (
        (pandas.reset_index(drop=True) == soma_arrow.to_pandas().reset_index(drop=True))
        .all()
        .all()
    )


def test_query_condition_select_columns():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02"

    qc = QueryCondition(condition)
    schema = tiledb.open(uri).schema

    sr = clib.SOMAReader(
        uri, query_condition=qc, schema=schema, column_names=["n_genes"]
    )
    sr.submit()
    arrow_table = sr.read_next()

    assert sr.results_complete()
    assert arrow_table.num_rows == 1332
    assert arrow_table.num_columns == 2


def test_query_condition_all_columns():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02"

    qc = QueryCondition(condition)
    schema = tiledb.open(uri).schema

    sr = clib.SOMAReader(uri, query_condition=qc, schema=schema)
    sr.submit()
    arrow_table = sr.read_next()

    assert sr.results_complete()
    assert arrow_table.num_rows == 1332
    assert arrow_table.num_columns == 7


def test_query_condition_reset():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02"

    qc = QueryCondition(condition)
    schema = tiledb.open(uri).schema

    sr = clib.SOMAReader(uri, query_condition=qc, schema=schema)
    sr.submit()
    arrow_table = sr.read_next()

    assert sr.results_complete()
    assert arrow_table.num_rows == 1332
    assert arrow_table.num_columns == 7

    # reset and submit new query with open array
    # ---------------------------------------------------------------
    condition = "percent_mito < 0.02"
    qc = QueryCondition(condition)
    sr.reset(column_names=["percent_mito"], query_condition=qc, schema=schema)

    sr.submit()
    arrow_table = sr.read_next()

    assert sr.results_complete()
    assert arrow_table.num_rows == 1306
    assert arrow_table.num_columns == 1


@pytest.mark.parametrize(
    "malformed_condition",
    [
        None,
        0,
        "",
        " ",
        b"",
        b" ",
        3.1415,
        object,
        {},
        ["a < 3"],
        "< a + 1",
        "a > val('1.0'",
        """
        a
        <
        2
        """,
        '"',
        "'",
        "attr(3) > 1attr(b) == 3",
    ],
)
def test_parsing_error_conditions(malformed_condition):
    """Conditions which should not parse."""

    with pytest.raises(SOMAError, match=r"Could not parse"):
        QueryCondition(malformed_condition)


@pytest.mark.parametrize(
    "malformed_condition",
    [
        "n_genes > b'abc'",
        "n_genes > 1+1",
        "percent_mito < n_counts",
        "attr('n_genes', 'foo') > 1",
        "somefunction('foo') == 'bar'",
        "a > ~0",
        "attr('n_g'+'enes') > 1",
        "attr(foo) > 1",
        "attr(1) > 1",
        "attr(True) > 1",
        "n_genes == val()",
        "attr() > 20",
        "n_genes < -val(-1)",
    ],
)
def test_eval_error_conditions(malformed_condition):
    """Conditions which should not evaluate (but WILL parse)"""
    uri = os.path.join(SOMA_URI, "obs")
    schema = tiledb.open(uri).schema

    # TODO: these raise the wrong error - it should be SOMAError. Change the test
    # when https://github.com/single-cell-data/TileDB-SOMA/issues/783 is fixed
    #
    with pytest.raises(RuntimeError):
        qc = QueryCondition(malformed_condition)
        sr = clib.SOMAReader(uri, query_condition=qc, schema=schema)
        sr.submit()
        sr.read_next()


if __name__ == "__main__":
    test_query_condition_select_columns()
