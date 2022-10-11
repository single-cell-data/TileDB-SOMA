#!/usr/bin/env python

import os
import pandas as pd
import pyarrow as pa
import tiledbsoma.libtiledbsoma as sc
import tiledb
from tiledbsoma.query_condition import QueryCondition


VERBOSE = True

TEST_DIR = os.path.dirname(__file__)
SOMA_URI = f"{TEST_DIR}/../../test/soco/pbmc3k_processed"

if VERBOSE:
    sc.config_logging("debug")


def pandas_query(uri, condition):
    sr = sc.SOMAReader(uri)
    sr.submit()
    arrow_table = sr.read_next()
    assert sr.results_complete()

    return arrow_table.to_pandas().query(condition)


def soma_query(uri, column_names, condition):
    qc = QueryCondition(condition)
    schema = tiledb.open(uri).schema

    sr = sc.SOMAReader(
        uri, column_names=column_names, query_condition=qc, schema=schema
    )
    sr.submit()
    arrow_table = sr.read_next()
    assert sr.results_complete()

    return arrow_table


def test_query_condition_int():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "n_genes > 500"

    pandas = pandas_query(uri, condition)

    column_names = ["n_genes"]
    soma_arrow = soma_query(uri, column_names, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_string():
    uri = os.path.join(SOMA_URI, "obs")
    condition = 'louvain == "NK cells"'

    pandas = pandas_query(uri, condition)

    column_names = ["louvain"]
    soma_arrow = soma_query(uri, column_names, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_float():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02"

    pandas = pandas_query(uri, condition)

    column_names = ["percent_mito"]
    soma_arrow = soma_query(uri, column_names, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_and():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02 and n_genes > 700"

    pandas = pandas_query(uri, condition)

    column_names = ["percent_mito", "n_genes"]
    soma_arrow = soma_query(uri, column_names, condition)

    assert len(pandas.index) == soma_arrow.num_rows


def test_query_condition_and_or():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "(percent_mito > 0.02 and n_genes > 700) or (percent_mito < 0.015 and louvain == 'B cells')"

    pandas = pandas_query(uri, condition)

    column_names = ["percent_mito", "n_genes", "louvain"]
    soma_arrow = soma_query(uri, column_names, condition)

    assert len(pandas.index) == soma_arrow.num_rows


if __name__ == "__main__":
    test_query_condition_and_or()
