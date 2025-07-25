from __future__ import annotations

import re
from contextlib import contextmanager, nullcontext
from pathlib import Path
from typing import Any, Union

import _pytest
import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from _pytest._code import ExceptionInfo
from _pytest.logging import LogCaptureFixture
from packaging.version import Version
from typing_extensions import TypeVar

import tiledbsoma

if Version(_pytest.__version__) < Version("8.4.0"):
    from _pytest.python_api import RaisesContext

    E = TypeVar("E", bound=BaseException, default=BaseException)
    MaybeRaisesReturn = Union[RaisesContext[E], ExceptionInfo[E], nullcontext]
else:
    from _pytest.raises import RaisesExc

    E = TypeVar("E", bound=BaseException, default=BaseException)
    MaybeRaisesReturn = Union[RaisesExc[E], ExceptionInfo[E], nullcontext]


from anndata import AnnData
from numpy import array_equal
from pandas._testing import assert_frame_equal, assert_series_equal
from scipy.sparse import spmatrix
from somacore import (
    AffineTransform,
    AxisQuery,
    CoordinateTransform,
    IdentityTransform,
    ScaleTransform,
    UniformScaleTransform,
)
from typeguard import suppress_type_checks


def assert_uns_equal(uns0, uns1):
    assert uns0.keys() == uns1.keys(), (
        f"extra keys: {uns0.keys() - uns1.keys()}, missing keys: {uns1.keys() - uns0.keys()}"
    )
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
    for k in d0:
        assert_array_equal(d0[k], d1[k])


def assert_array_equal(a0, a1):
    assert type(a0) is type(a1)
    if isinstance(a0, np.ndarray):
        assert array_equal(a0, a1)
    elif isinstance(a0, spmatrix):
        assert type(a0) is type(a1)
        assert a0.shape == a1.shape
        assert (a0 != a1).nnz == 0
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


def assert_transform_equal(actual: CoordinateTransform, expected: CoordinateTransform) -> None:
    assert actual.input_axes == expected.input_axes
    assert actual.output_axes == expected.output_axes
    if isinstance(expected, IdentityTransform):
        assert isinstance(actual, IdentityTransform)
    elif isinstance(expected, UniformScaleTransform):
        assert isinstance(actual, UniformScaleTransform)
        assert actual.scale == expected.scale
    elif isinstance(expected, ScaleTransform):
        assert isinstance(actual, ScaleTransform)
        np.testing.assert_array_equal(actual.scale_factors, expected.scale_factors)
    elif isinstance(expected, AffineTransform):
        assert isinstance(actual, AffineTransform)
        np.testing.assert_array_equal(actual.augmented_matrix, expected.augmented_matrix)
    else:
        assert False


def parse_col(col_str: str) -> tuple[str | None, list[str]]:
    """Parse a "column string" of the form ``val1,val2,...`` or ``name=val1,val2,...``."""
    pcs = col_str.split("=")
    if len(pcs) == 1:
        return None, col_str.split(",")
    if len(pcs) == 2:
        name, vals_str = pcs
        vals = vals_str.split(",")
        return name, vals
    raise ValueError(f"Invalid column string: {col_str}")


def make_pd_df(index_str: str | None = None, **cols) -> pd.DataFrame:
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
PROJECT_ROOT = PY_ROOT.parent.parent
TESTDATA = PY_ROOT / "testdata"
ROOT_DATA_DIR = PROJECT_ROOT / "data"


@contextmanager
def raises_no_typeguard(exc: type[Exception], *args: Any, **kwargs: Any):
    """
    Temporarily suppress typeguard checks in order to verify a runtime exception is raised.

    Otherwise, most errors end up manifesting as ``TypeCheckError``s, during tests (thanks to
    ``typeguard``'s import hook).
    """
    with suppress_type_checks(), pytest.raises(exc, *args, **kwargs):
        yield


# Alias for several types that can be used to specify expected exceptions in `maybe_raises`
Err = Union[str, type[E], tuple[type[E], str]]


def maybe_raises(expected_exception: Err | None, *args: Any, **kwargs: Any) -> MaybeRaisesReturn:
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
                raise ValueError("Cannot specify 'match' in both kwargs and `expected_exception`")
            kwargs["match"] = match
        return pytest.raises(exc, *args, **kwargs)
    return nullcontext()


def verify_logs(caplog: LogCaptureFixture, expected_logs: list[str] | None) -> None:
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


def filter(value_filter: str) -> AxisQuery:
    """Shorthand for creating an ``AxisQuery`` with a value_filter, in tests."""
    return AxisQuery(value_filter=value_filter)


def create_basic_object(soma_type, uri, **kwargs) -> tiledbsoma.SOMAObject:
    """Create a basic SOMA object of the requested type."""

    if soma_type == "SOMAExperiment":
        return tiledbsoma.Experiment.create(uri, **kwargs)
    if soma_type == "SOMAMeasurement":
        return tiledbsoma.Measurement.create(uri, **kwargs)
    if soma_type == "SOMACollection":
        return tiledbsoma.Collection.create(uri, **kwargs)
    if soma_type == "SOMAScene":
        return tiledbsoma.Scene.create(uri, **kwargs)
    if soma_type == "SOMADataFrame":
        kwargs.setdefault("schema", pa.schema([pa.field("myint", pa.int64())]))
        return tiledbsoma.DataFrame.create(uri, **kwargs)
    if soma_type == "SOMAGeometryDataFrame":
        kwargs.setdefault("schema", pa.schema([("quality", pa.float32())]))
        return tiledbsoma.GeometryDataFrame.create(uri, **kwargs)
    if soma_type == "SOMAPointCloudDataFrame":
        kwargs.setdefault("schema", pa.schema([("x", pa.float64()), ("y", pa.float64())]))
        return tiledbsoma.PointCloudDataFrame.create(uri, **kwargs)
    if soma_type == "SOMADenseNDArray":
        kwargs.setdefault("type", pa.float64())
        kwargs.setdefault("shape", (100, 100))
        return tiledbsoma.DenseNDArray.create(uri, **kwargs)
    if soma_type == "SOMASparseNDArray":
        kwargs.setdefault("type", pa.float64())
        kwargs.setdefault("shape", (100, 100))
        return tiledbsoma.SparseNDArray.create(uri, **kwargs)
    if soma_type == "SOMAMultiscaleImage":
        kwargs.setdefault("type", pa.uint8())
        kwargs.setdefault("level_shape", (3, 64, 128))
        return tiledbsoma.MultiscaleImage.create(uri, **kwargs)

    raise f"Internal error: Unexepcted soma type '{soma_type}'."
