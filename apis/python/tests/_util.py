from contextlib import contextmanager
from typing import Type

import pytest
from typeguard import suppress_type_checks


@contextmanager
def raises(exc: Type[Exception]):
    """Temporarily suppress typeguard checks in order to verify a runtime exception is raised.

    Otherwise, most errors end up manifesting as ``TypeCheckError``s, during tests (thanks to
    ``typeguard``'s import hook).
    """
    with suppress_type_checks():
        with pytest.raises(exc):
            yield
