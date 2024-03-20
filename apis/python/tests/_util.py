from contextlib import contextmanager
from typing import Any, Type

import pytest
from typeguard import suppress_type_checks


@contextmanager
def raises_no_typeguard(exc: Type[Exception], *args: Any, **kwargs: Any):
    """Temporarily suppress typeguard checks in order to verify a runtime exception is raised.

    Otherwise, most errors end up manifesting as ``TypeCheckError``s, during tests (thanks to
    ``typeguard``'s import hook).
    """
    with suppress_type_checks():
        with pytest.raises(exc, *args, **kwargs):
            yield
