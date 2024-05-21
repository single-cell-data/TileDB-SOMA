from contextlib import contextmanager, nullcontext
from pathlib import Path
from typing import Any, Tuple, Type, Union

import pytest
from _pytest._code import ExceptionInfo
from _pytest.python_api import E, RaisesContext
from typeguard import suppress_type_checks

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


def maybe_raises(
    expected_exception: Union[None, Type[E], Tuple[Type[E], ...]],
    *args: Any,
    **kwargs: Any
) -> Union[RaisesContext[E], ExceptionInfo[E]]:
    """
    Wrapper around ``pytest.raises`` that accepts None (signifying no exception should be raised).
    This is only necessary since ``pytest.raises`` does not itself accept None, so we are
    decorating.

    Useful in test cases that are parameterized to test both valid and invalid inputs.
    """
    return (
        nullcontext()
        if expected_exception is None
        else pytest.raises(expected_exception, *args, **kwargs)
    )
