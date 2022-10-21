#!/usr/bin/env python

import os

import tiledb

import tiledbsoma.libtiledbsoma as clib
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


def test_query_condition_int():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "n_genes > 500"

    pandas = pandas_query(uri, condition)

    soma_arrow = soma_query(uri, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_string():
    uri = os.path.join(SOMA_URI, "obs")
    condition = 'louvain == "NK cells"'

    pandas = pandas_query(uri, condition)

    soma_arrow = soma_query(uri, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_float():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02"

    pandas = pandas_query(uri, condition)

    soma_arrow = soma_query(uri, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_bool():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "is_b_cell == True"

    pandas = pandas_query(uri, condition)

    soma_arrow = soma_query(uri, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_and():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02 and n_genes > 700"

    pandas = pandas_query(uri, condition)

    soma_arrow = soma_query(uri, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_and_or():
    uri = os.path.join(SOMA_URI, "obs")
    condition = '(percent_mito > 0.02 and n_genes > 700) or (percent_mito < 0.015 and louvain == "B cells")'

    pandas = pandas_query(uri, condition)

    soma_arrow = soma_query(uri, condition)

    assert len(pandas.index) == soma_arrow.num_rows


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


if __name__ == "__main__":
    test_query_condition_select_columns()
