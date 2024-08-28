import logging
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
from _pytest.logging import LogCaptureFixture

import tiledbsoma
from tiledbsoma import Experiment
from tiledbsoma.io._common import UnsDict, UnsMapping
from tiledbsoma.io.update_uns import Strict, _update_uns

from tests._util import Err, assert_uns_equal, make_pd_df, maybe_raises, verify_logs
from tests.parametrize_cases import parametrize_cases
from tests.test_basic_anndata_io import TEST_UNS, make_uns_adata

ValidUpdates = Union[None, str, List[str], Dict[str, "ValidUpdates"]]
Logs = Optional[List[str]]


@dataclass
class Case:
    """A test case for `update_uns`.

    :param id: Human-readable identifier for the test case.
    :param uns_updates: Updates to attempt to apply (to the default ``TEST_UNS`` dict).
    :param strict: How to handle conflicts with existing ``uns`` keys. Default (``True``) performs
        an initial "dry run" pass to verify all updates can be applied, then a second pass actually
        applying the updates. Pass ``"debug"``, ``"info"``, or ``"warn"`` to instead skip
        conflicting keys, and log a message at the corresponding level. See :func:`update_uns` for
        more details.
    :param err: If present, verify that :func:`update_uns` raises an ``Exception`` containing this
        regex.
    :param valid_updates: If present, only these updates are expected to have been applied. If
        ``strict=True`` (the default), all updates must be applied (``valid_updates=None``), or
        no updates (``valid_updates=[]``), based on presence or absence of an expected ``err``.
    :param logs: If present, verify that these log messages are present in the captured logs. Used
        only when ``strict`` is a log level (``"debug"``, ``"info"``, ``"warn"``).
    """

    id: str
    uns_updates: UnsMapping
    strict: Strict = True
    err: Optional[Err] = None
    valid_updates: ValidUpdates = None
    logs: Logs = None


def case(
    id: str,
    err: Optional[Err] = None,
    valid_updates: ValidUpdates = None,
    logs: Logs = None,
    strict: Strict = True,
    **uns_updates,
) -> Case:
    """Helper for construction a :class:`Case`.

    ``err``s are verified to be ``ValueError``s, by default, and "kwargs" become ``uns_updates``.
    """
    if isinstance(err, str):
        err = (ValueError, err)
    return Case(
        id=id,
        uns_updates=uns_updates,
        err=err,
        strict=strict,
        valid_updates=valid_updates,
        logs=logs,
    )


# fmt: off
@parametrize_cases([
    case(
        "Update one scalar",
        int_scalar=11,
    ),
    case(
        "Update scalar, change type",
        int_scalar="aaa",
    ),
    case(
        "Update multiple scalar values, add a new DataFrame",
        int_scalar=11,
        float_scalar=2.2,
        string_scalar="HELLO 2",
        new_df=make_pd_df(a="1,2,3"),
    ),
    case(
        "Update multiple scalar values, including inside a Collection",
        int_scalar=11,
        float_scalar=2.2,
        string_scalar="HELLO 2",
        new_df=make_pd_df(a="1,2,3"),
        strings=dict(
            aaa="AAA",
            nnn=111,
        )
    ),
    case(
        "Overwrite np.array inside collection (raise)",
        strings=dict(
           string_np_ndarray_1d=np.asarray(list("abc")),
        ),
        err=r"ms/RNA/uns/strings\[string_np_ndarray_1d]: already exists \(type DataFrame\), refusing to overwrite with \['a' 'b' 'c']"
    ),
    case(
        "Overwrite np.array inside collection (skip)",
        strings=dict(
            string_np_ndarray_1d=np.asarray(list("abc")),
        ),
        strict="info",
        valid_updates=[],
        logs=[r"ms/RNA/uns/strings\[string_np_ndarray_1d]: already exists \(type DataFrame\), refusing to overwrite with \['a' 'b' 'c']"]
    ),
    case(
        "No partial updates inside collection",
        strings=dict(
            foo="FOO",
            string_np_ndarray_1d=np.asarray(list("abc")),
        ),
        err=r"ms/RNA/uns/strings\[string_np_ndarray_1d]: already exists \(type DataFrame\), refusing to overwrite with \['a' 'b' 'c']"
    ),
    case(
        "Partial update inside collection",
        strings=dict(
            foo="FOO",
            string_np_ndarray_1d=np.asarray(list("abc")),
        ),
        strict="info",
        valid_updates={"strings": "foo"},
        logs=[r"ms/RNA/uns/strings\[string_np_ndarray_1d]: already exists \(type DataFrame\), refusing to overwrite with \['a' 'b' 'c']"]
    ),
    case(
        "Overwrite scalar with DataFrame",
        int_scalar=make_pd_df(a="1,2,3"),
        err=(
                r"ms/RNA/uns\[int_scalar]: already exists \(type int\), refusing to overwrite with    a\n"
                r"0  1\n"
                r"1  2\n"
                r"2  3"
        ),
    ),
    case(
        "Overwrite existing DataFrame (raise)",
        int_scalar=22,
        pd_df_indexed=TEST_UNS["pd_df_indexed"].copy(),
        err=(
            r"ms/RNA/uns\[pd_df_indexed]: already exists \(type DataFrame\), refusing to overwrite with   column_1\n"
            r"0        d\n"
            r"1        e\n"
            r"2        f"
        ),
    ),
    case(
        "Overwrite existing DataFrame (skip)",
        int_scalar=22,
        pd_df_indexed=TEST_UNS["pd_df_indexed"].assign(column_1=list("ghi")),
        strict="info",
        valid_updates=['int_scalar'],  # `int_scalar` is updated, `pd_df_indexed` is not ("info" msg logged below)
        logs=[
            r"ms/RNA/uns\[pd_df_indexed]: already exists \(type DataFrame\), refusing to overwrite with   column_1\n"
            r"0        g\n"
            r"1        h\n"
            r"2        i"
        ],
    ),
])
# fmt: on
def test_update_uns(
    caplog: LogCaptureFixture,
    tmp_path: Path,
    uns_updates: UnsDict,
    strict: Strict,
    err: Optional[Err],
    valid_updates: ValidUpdates,
    logs: Logs,
):
    caplog.set_level(logging.INFO)
    soma_uri, adata = make_uns_adata(tmp_path)

    with Experiment.open(soma_uri, "w") as exp:
        with maybe_raises(err):
            _update_uns(exp, uns_updates, measurement_name="RNA", strict=strict)

    verify_logs(caplog, logs)

    if err:
        # In all cases, an error during `update_uns` should result in no updates being applied (nor
        # should any non-empty `valid_updates` have been provided as part of the test-case "spec").
        assert valid_updates is None or valid_updates == []
        valid_updates = []

    with Experiment.open(soma_uri) as exp:
        adata2 = tiledbsoma.io.to_anndata(exp, measurement_name="RNA")

    expected = deepcopy(TEST_UNS)
    merge_updates(expected, uns_updates, valid_updates)
    assert_uns_equal(adata2.uns, expected)


def merge_updates(
    uns: UnsDict,
    updates: UnsDict,
    valid_updates: ValidUpdates,
) -> None:
    """Apply a subset of ``updates`` to an ``uns`` dict, filtering based on ``valid_updates``.

    - If ``valid_updates`` is ``None``, all updates are applied.
    - If ``valid_updates`` is a ``list`` of strings, ``updates`` must be a ``dict``, and only keys
        from the ``list`` are applied.
    - If ``valid_updates`` is a ``dict``, its keys must be present in ``updates``, and filtering is
        applied recursively, with ``valid_updates``' values applying to corresponding items in
        ``updates``.
    """
    if isinstance(valid_updates, str):
        valid_updates = [valid_updates]
    for k, v in updates.items():
        if isinstance(v, dict):
            if valid_updates is None:
                merge_updates(uns[k], v, None)
            elif isinstance(valid_updates, list):
                if k in valid_updates:
                    merge_updates(uns[k], v, None)
            else:
                assert isinstance(valid_updates, dict)
                if k in valid_updates:
                    merge_updates(uns[k], v, valid_updates[k])
        else:
            if valid_updates is None:
                uns[k] = v
            else:
                assert isinstance(valid_updates, list)
                if k in valid_updates:
                    uns[k] = v
