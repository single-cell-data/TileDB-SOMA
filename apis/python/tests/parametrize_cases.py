from dataclasses import asdict, fields
from inspect import getfullargspec
from typing import List

import pytest


def parametrize_cases(cases: List):
    """Parametrize a test with a list of test cases (each an instance of a dataclass).

    The cases are expected to have a ``name`` string, which is used as the test-case "ID".
    """

    cls = cases[0].__class__
    for case in cases:
        if case.__class__ is not cls:
            raise ValueError(
                f"Expected all cases to be of type {cls}, but found {case.__class__}"
            )

    def wrapper(fn):
        # Test-case IDs
        ids = [case.name for case in cases]
        # Convert each case to a "values" array; also filter and reorder to match kwargs expected
        # by the wrapped "test_*" function.
        spec = getfullargspec(fn)
        field_names = [f.name for f in fields(cls)]
        names = [arg for arg in spec.args if arg in field_names]
        values = [
            {name: rt_dict[name] for name in names}.values()
            for rt_dict in [asdict(case) for case in cases]
        ]
        # Delegate to PyTest `parametrize`
        return pytest.mark.parametrize(
            names,  # kwarg names
            values,  # arg value lists
            ids=ids,  # test-case names
        )(fn)

    return wrapper
