import re
from contextlib import contextmanager, nullcontext
from pathlib import Path
from typing import Any, List, Optional, Tuple, Type, Union

import numpy as np
import pandas as pd
import pytest
from _pytest._code import ExceptionInfo
from _pytest.logging import LogCaptureFixture
from _pytest.python_api import E, RaisesContext
from anndata import AnnData
from numpy import array_equal
from pandas._testing import assert_frame_equal, assert_series_equal
from scipy.sparse import spmatrix
from typeguard import suppress_type_checks


def assert_uns_equal(uns0, uns1):
    assert (
        uns0.keys() == uns1.keys()
    ), f"extra keys: {uns0.keys() - uns1.keys()}, missing keys: {uns1.keys() - uns0.keys()}"
    for k, v0 in uns0.items():
        try:
            v1 = uns1[k]
            if isinstance(v0, dict):
                assert isinstance(v1, dict)
                assert_uns_equal(v0, v1)
            elif isinstance(v0, list):
                assert isinstance(v1, list)
                assert len(v0) == len(v1)
                for e0, e1 in zip(v0, v1):
                    assert_uns_equal(e0, e1)
            elif isinstance(v0, pd.DataFrame):
                assert_frame_equal(v0, v1)
            elif isinstance(v0, pd.Series):
                assert_series_equal(v0, v1)
            elif isinstance(v0, np.ndarray):
                assert array_equal(v0, v1)
            elif isinstance(v0, (int, float, str, bool)):
                assert v0 == v1, f"{v0} != {v1}"
            else:
                raise ValueError(f"Unsupported type: {type(v0)}")
        except AssertionError:
            raise AssertionError(f"assert_uns_equal: key {k} mismatched")


def assert_array_dicts_equal(d0, d1):
    assert d0.keys() == d1.keys()
    for k in d0.keys():
        assert_array_equal(d0[k], d1[k])


def assert_array_equal(a0, a1):
    assert type(a0) is type(a1)
    if isinstance(a0, np.ndarray):
        assert array_equal(a0, a1)
    elif isinstance(a0, spmatrix):
        assert array_equal(a0.todense(), a1.todense())
    else:
        raise ValueError(f"Unsupported type: {type(a0)}")


def assert_adata_equal(ad0: AnnData, ad1: AnnData):
    assert_frame_equal(ad0.obs, ad1.obs)
    assert_frame_equal(ad0.var, ad1.var)
    assert_uns_equal(ad0.uns, ad1.uns)
    assert_array_equal(ad0.X, ad1.X)
    assert_array_dicts_equal(ad0.obsm, ad1.obsm)
    assert_array_dicts_equal(ad0.varm, ad1.varm)
    assert_array_dicts_equal(ad0.obsp, ad1.obsp)
    assert_array_dicts_equal(ad0.varp, ad1.varp)


def parse_col(col_str: str) -> Tuple[Optional[str], List[str]]:
    """Parse a "column string" of the form ``val1,val2,...`` or ``name=val1,val2,...``."""
    pcs = col_str.split("=")
    if len(pcs) == 1:
        return None, col_str.split(",")
    elif len(pcs) == 2:
        name, vals_str = pcs
        vals = vals_str.split(",")
        return name, vals
    else:
        raise ValueError(f"Invalid column string: {col_str}")


def make_df(index_str: Optional[str] = None, **cols) -> pd.DataFrame:
    """DataFrame construction helper, for tests.

    - index and columns are provided as strings of the form ``name=val1,val2,...``.
    - ``name=`` is optional for the initial (``index_str``) arg.
    """
    cols = dict([(col, parse_col(col_str)[1]) for col, col_str in cols.items()])
    index = None
    index_name = None
    if index_str:
        index_name, index = parse_col(index_str)
    df = pd.DataFrame(cols, index=index)
    df.index.name = index_name
    return df


HERE = Path(__file__).parent
PY_ROOT = HERE.parent
TESTDATA = PY_ROOT / "testdata"


@contextmanager
def raises_no_typeguard(exc: Type[Exception], *args: Any, **kwargs: Any):
    """
    Temporarily suppress typeguard checks in order to verify a runtime exception is raised.

    Otherwise, most errors end up manifesting as ``TypeCheckError``s, during tests (thanks to
    ``typeguard``'s import hook).
    """
    with suppress_type_checks():
        with pytest.raises(exc, *args, **kwargs):
            yield


# Alias for several types that can be used to specify expected exceptions in `maybe_raises`
Err = Union[str, Type[E], Tuple[Type[E], str]]


def maybe_raises(
    expected_exception: Optional[Err],
    *args: Any,
    **kwargs: Any,
) -> Union[RaisesContext[E], ExceptionInfo[E], nullcontext]:
    """
    Wrapper around ``pytest.raises`` that additionally accepts ``None`` (signifying no exception
    should be raised), a string (signifying a message match) or a tuple of (exception, message).

    Useful in test cases that are parameterized to test both valid and invalid inputs.
    """
    if expected_exception is not None:
        if isinstance(expected_exception, type):
            exc = expected_exception
        else:
            if isinstance(expected_exception, str):
                exc, match = Exception, expected_exception
            else:
                exc, match = expected_exception
            if "match" in kwargs:
                raise ValueError(
                    "Cannot specify 'match' in both kwargs and `expected_exception`"
                )
            kwargs["match"] = match
        return pytest.raises(exc, *args, **kwargs)
    else:
        return nullcontext()


def verify_logs(caplog: LogCaptureFixture, expected_logs: Optional[List[str]]) -> None:
    """Verify that expected log messages are present in a pytest "caplog" (captured logs) fixture."""
    if not expected_logs:
        return

    for expected_log in expected_logs:
        found = False
        for record in caplog.records:
            log_msg = record.getMessage()
            if re.search(expected_log, log_msg):
                found = True
                break
        if not found:
            raise AssertionError(f"Expected log message not found: {expected_log}")
