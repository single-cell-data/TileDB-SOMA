#!/usr/bin/env python

import os
import pathlib

import pytest

import tiledbsoma.pytiledbsoma as clib
from tiledbsoma import DataFrame
from tiledbsoma._exception import SOMAError
from tiledbsoma._query_condition import QueryCondition

VERBOSE = False

TEST_DIR = pathlib.Path(__file__).parent
SOMA_URI = f"{TEST_DIR}/../../../data/soco/pbmc3k_processed"

if VERBOSE:
    clib.config_logging("debug")


def pandas_query(uri, condition):
    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())
    arrow_table = mq.next()

    with pytest.raises(StopIteration):
        mq.next()

    return arrow_table.to_pandas().query(condition)


def soma_query(uri, condition):
    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_condition(QueryCondition(condition), sr.schema)
    arrow_table = mq.next()

    with pytest.raises(StopIteration):
        mq.next()

    return arrow_table


@pytest.mark.parametrize(
    "soma_pd_condition",
    [
        # types
        ("int == None or int == 1", "(int != int) | (int == 1)"),
        ("int == None and bool == None", "int.isna() & bool.isna()"),
        ("int == None and ord != 'g1'", "int.isna() & ord != 'g1'"),
    ],
)
def test_query_with_combo_none_condition(tmp_path, soma_pd_condition):
    import pandas as pd
    import pyarrow as pa

    import tiledbsoma as soma

    uri = tmp_path.as_posix()
    asch = pa.schema(
        [
            pa.field("int", pa.int32()),
            pa.field("bool", pa.bool_()),
            pa.field("ord", pa.dictionary(pa.int64(), pa.string())),
        ],
        metadata={
            "int": "nullable",
            "bool": "nullable",
            "ord": "nullable",
        },
    )

    pydict = {}
    pydict["soma_joinid"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    pydict["int"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["bool"] = [True, True, True, False, True, False, None, False, None, None]
    pydict["ord"] = pd.Categorical(["g1", "g2", "g3", None, "g2", "g3", "g1", None, "g3", "g1"])
    data = pa.Table.from_pydict(pydict, schema=asch.insert(0, pa.field("soma_joinid", pa.int64())))

    with soma.DataFrame.create(uri, schema=asch, domain=[[0, 9]]) as sdf:
        sdf.write(data)

    soma_arrow = soma_query(uri, soma_pd_condition[0])
    pandas = pandas_query(uri, soma_pd_condition[1])
    assert len(pandas.index) == soma_arrow.num_rows


@pytest.mark.parametrize(
    "soma_pd_condition",
    [
        # types
        ("ord == None", "ord.isna()"),
        ("ord != None", "ord.notna()"),
        ("ord <= None", "ord < ord"),
        ("ord < None", "ord < ord"),
        ("ord >= None", "ord >= ord"),
        ("ord >= None", "ord >= ord"),
    ],
)
def test_query_with_none_condition(tmp_path, soma_pd_condition):
    """Compare condition to result from Pandas - they should match"""

    import pandas as pd
    import pyarrow as pa

    import tiledbsoma as soma

    uri = tmp_path.as_posix()
    asch = pa.schema(
        [
            pa.field("int", pa.int32()),
            pa.field("bool", pa.bool_()),
            pa.field("ord", pa.dictionary(pa.int64(), pa.string())),
        ],
        metadata={
            "int": "nullable",
            "bool": "nullable",
            "ord": "nullable",
        },
    )

    pydict = {}
    pydict["soma_joinid"] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    pydict["int"] = [1, 2, 3, 4, 5, 6, None, 8, None, None]
    pydict["bool"] = [True, True, True, False, True, False, None, False, None, None]
    pydict["ord"] = pd.Categorical(["g1", "g2", "g3", None, "g2", "g3", "g1", None, "g3", "g1"])
    data = pa.Table.from_pydict(pydict, schema=asch.insert(0, pa.field("soma_joinid", pa.int64())))

    with soma.DataFrame.create(uri, schema=asch, domain=[[0, 9]]) as sdf:
        sdf.write(data)

    if any(op in soma_pd_condition[0] for op in ["==", "!="]):
        soma_arrow = soma_query(uri, soma_pd_condition[0])
        pandas = pandas_query(uri, soma_pd_condition[1])
        assert len(pandas.index) == soma_arrow.num_rows
    else:
        with pytest.raises(SOMAError) as e:
            soma_arrow = soma_query(uri, soma_pd_condition[0])
            assert "Null value can only be used with equality operators" in str(e)


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
        # not ops
        "not n_genes == 480",
        "not n_genes > 500",
        "not n_genes < 500",
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
    assert pandas.empty or (
        (pandas.reset_index(drop=True) == soma_arrow.to_pandas().reset_index(drop=True)).all().all()
    )


@pytest.mark.parametrize(
    "condition, pandas_equivalent_condition",
    [
        ["attr('n_genes') == 480", "n_genes == 480"],
        ["n_genes > val(100)", "n_genes > 100"],
        ["n_genes > val(100.0)", "n_genes > 100.0"],
        ["is_b_cell == val(True)", "is_b_cell == True"],
        ["louvain == val('B cells')", "louvain == 'B cells'"],
        [
            "louvain == 'B cells' or louvain == 'NK cells'",
            "louvain == 'B cells' or louvain == 'NK cells'",
        ],
        [
            "attr('louvain') == 'B cells' or attr('louvain') == 'NK cells'",
            "louvain == 'B cells' or louvain == 'NK cells'",
        ],
        [
            "louvain in ['B cells', 'NK cells']",
            "louvain == 'B cells' or louvain == 'NK cells'",
        ],
        [
            "attr('louvain') in ['B cells', 'NK cells']",
            "louvain == 'B cells' or louvain == 'NK cells'",
        ],
    ],
)
def test_query_condition_extensions(condition, pandas_equivalent_condition):
    """Test conditions which Pandas does not support directly - extensions"""
    uri = os.path.join(SOMA_URI, "obs")
    pandas = pandas_query(uri, pandas_equivalent_condition)
    soma_arrow = soma_query(uri, condition)
    assert len(pandas.index) == soma_arrow.num_rows
    assert (pandas.reset_index(drop=True) == soma_arrow.to_pandas().reset_index(drop=True)).all().all()


def test_query_condition_select_columns():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02"

    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())
    mq.select_columns(["n_genes"])
    mq.set_condition(QueryCondition(condition), sr.schema)
    arrow_table = mq.next()

    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == 1332
    assert arrow_table.num_columns == 2


def test_query_condition_all_columns():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02"

    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_condition(QueryCondition(condition), sr.schema)
    arrow_table = mq.next()

    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == 1332
    assert arrow_table.num_columns == 7


def test_query_condition_reset():
    uri = os.path.join(SOMA_URI, "obs")
    condition = "percent_mito > 0.02"

    sr = clib.SOMAArray(uri)
    mq = clib.ManagedQuery(sr, sr.context())
    mq.set_condition(QueryCondition(condition), sr.schema)
    arrow_table = mq.next()

    with pytest.raises(StopIteration):
        mq.next()

    assert arrow_table.num_rows == 1332
    assert arrow_table.num_columns == 7

    # reset and submit new query with open array
    # ---------------------------------------------------------------
    condition = "percent_mito < 0.02"

    mq = clib.ManagedQuery(sr, sr.context())
    mq.select_columns(["percent_mito"])
    mq.set_condition(QueryCondition(condition), sr.schema)
    arrow_table = mq.next()

    with pytest.raises(StopIteration):
        mq.next()

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
        "not > a",
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
        "louvain in []",
    ],
)
def test_eval_error_conditions(malformed_condition):
    """Conditions which should not evaluate (but WILL parse)"""
    uri = os.path.join(SOMA_URI, "obs")
    qc = QueryCondition(malformed_condition)

    with pytest.raises(SOMAError):
        sr = clib.SOMAArray(uri)
        mq = clib.ManagedQuery(sr, sr.context())
        mq.set_condition(qc, sr.schema)

    with pytest.raises(SOMAError):
        # test function directly for codecov
        qc.init_query_condition(sr.schema, [])
        qc.init_query_condition(sr.schema, ["bad_query_attr"])


@pytest.mark.parametrize(
    "expression_and_message",
    [
        ["foo is True", "the `is` operator is not supported"],
        ["foo is not True", "the `is not` operator is not supported"],
        [
            "foo &&& bar",
            "Could not parse the given QueryCondition statement: foo &&& bar",
        ],
        [
            "louvain == leukocyte",
            "Incorrect type for comparison value.* right-hand sides must be constant expressions, not variables",
        ],
        # Minus sign looks like an operator to the Python parser:
        [
            "louvain == leuko-cyte",
            "Unable to parse expression component.* did you mean to quote it as a string",
        ],
        # Test the dot "operator" (valid name in R, not in Python)
        [
            'lou.vain == "leukocyte"',
            "if your attribute name has a dot in it",
        ],
    ],
)
def test_query_condition_syntax_handling(expression_and_message):
    uri = os.path.join(SOMA_URI, "obs")
    expression, message = expression_and_message
    with DataFrame.open(uri) as obs, pytest.raises(SOMAError, match=message):
        obs.read(value_filter=expression).concat()


if __name__ == "__main__":
    test_query_condition_select_columns()
